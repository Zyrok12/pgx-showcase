#!/usr/bin/env python3
"""
review_triage_dl.py — learned variant-quality model that TRIAGES the GATK REVIEW tier.

Our tiered filter emits PASS / REVIEW / FAIL. We proved REVIEW holds real variants (including it
lifted GIAB F1), but it's a mixed bag — some true calls, some artifacts. Hard thresholds can't
separate them cleanly. This model scores each REVIEW-tier call -> P(real) and reclassifies it:
rescue the likely-real ones to PASS, drop the likely-artifacts to FAIL, keep the genuinely
ambiguous as REVIEW. It's a lightweight, tabular cousin of DeepVariant — the right DL/ML fit
because the training signal is HUGE and free: every variant in a GIAB run is TP/FP-labelable
against the NIST truth set, and the mapping (pileup/annotation features -> real vs artifact) is
exactly what a calibrated tree ensemble learns well.

Why this matters for the product: it directly serves "don't miss variants" (rescues real REVIEW
calls regression-style filters would drop) WHILE cutting false positives — measurable as held-out
GIAB F1 vs the two dumb baselines (keep-all-REVIEW, drop-all-REVIEW).

Features (GATK VCF annotations): QUAL, QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum, BaseQRankSum,
DP, GQ, allele-balance (from AD), is_indel, indel_len. Validation: split BY CHROMOSOME (no leakage).
Degrades gracefully: needs pysam for real VCFs; `--demo` runs the full train/validate loop on
synthetic features to prove the plumbing.
"""
from __future__ import annotations

import argparse
import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

try:
    import pysam
    _HAVE_PYSAM = True
except Exception:
    _HAVE_PYSAM = False

try:
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.metrics import roc_auc_score, average_precision_score
    _HAVE_SKLEARN = True
except Exception:
    _HAVE_SKLEARN = False

INFO_FIELDS = ["QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum", "BaseQRankSum", "DP"]
FEATURES = ["QUAL"] + INFO_FIELDS + ["GQ", "allele_balance", "is_indel", "indel_len", "n_alt"]


# ----------------------------------------------------------------- VCF feature extraction

def extract_vcf_features(vcf_path: str, sample: Optional[str] = None) -> pd.DataFrame:
    if not _HAVE_PYSAM:
        raise RuntimeError("pysam required to read VCFs")
    vf = pysam.VariantFile(vcf_path)
    samples = list(vf.header.samples)
    smp = sample or (samples[0] if samples else None)
    rows = []
    for rec in vf:
        alt = rec.alts[0] if rec.alts else ""
        ref = rec.ref or ""
        is_indel = int(len(ref) != len(alt))
        f = {"chrom": rec.chrom, "pos": rec.pos, "ref": ref, "alt": alt,
             "filter": ";".join(rec.filter.keys()) if len(rec.filter.keys()) else "PASS",
             "QUAL": rec.qual if rec.qual is not None else np.nan,
             "is_indel": is_indel, "indel_len": abs(len(ref) - len(alt)),
             "n_alt": len(rec.alts) if rec.alts else 1}
        for k in INFO_FIELDS:
            v = rec.info.get(k)
            f[k] = float(v[0]) if isinstance(v, tuple) else (float(v) if v is not None else np.nan)
        gq, ab = np.nan, np.nan
        if smp is not None and smp in rec.samples:
            s = rec.samples[smp]
            gq = float(s["GQ"]) if "GQ" in s and s["GQ"] is not None else np.nan
            ad = s.get("AD")
            if ad and len(ad) >= 2 and sum(ad) > 0:
                ab = ad[1] / sum(ad)         # alt fraction; ~0.5 het / ~1.0 hom for real, skewed for artifacts
        f["GQ"], f["allele_balance"] = gq, ab
        rows.append(f)
    return pd.DataFrame(rows)


def label_against_truth(feat: pd.DataFrame, truth_vcf: str) -> pd.DataFrame:
    """Label each call TP(1)/FP(0) by exact chrom:pos:ref:alt match against a truth VCF."""
    if not _HAVE_PYSAM:
        raise RuntimeError("pysam required")
    truth = set()
    for rec in pysam.VariantFile(truth_vcf):
        for a in (rec.alts or [""]):
            truth.add((rec.chrom, rec.pos, rec.ref, a))
    feat = feat.copy()
    feat["label"] = [int((c, p, r, a) in truth) for c, p, r, a in
                     zip(feat.chrom, feat.pos, feat.ref, feat.alt)]
    return feat


# ----------------------------------------------------------------- model

@dataclass
class TriageConfig:
    rescue_threshold: float = 0.70     # P(real) >= -> rescue REVIEW to PASS
    drop_threshold: float = 0.30       # P(real) <  -> drop REVIEW to FAIL ; between = keep REVIEW
    seed: int = 0


class ReviewTriageModel:
    def __init__(self, cfg: Optional[TriageConfig] = None):
        if not _HAVE_SKLEARN:
            raise RuntimeError("scikit-learn required")
        self.cfg = cfg or TriageConfig()
        self.model = None
        self.feat_cols = FEATURES

    def _prep(self, df):
        X = df.reindex(columns=self.feat_cols).astype(float)
        return np.nan_to_num(X.values, nan=0.0)

    def fit(self, df: pd.DataFrame):
        X, y = self._prep(df), df["label"].astype(int).values
        base = HistGradientBoostingClassifier(max_depth=3, learning_rate=0.05, max_iter=300,
                                              l2_regularization=1.0, random_state=self.cfg.seed)
        self.model = CalibratedClassifierCV(base, method="isotonic", cv=3) \
            if min(np.bincount(y)) >= 5 else base
        self.model.fit(X, y)
        return self

    def score(self, df: pd.DataFrame) -> np.ndarray:
        return self.model.predict_proba(self._prep(df))[:, 1]

    def triage_filter(self, p_real: float, current: str) -> str:
        if "REVIEW" not in current.upper():
            return current
        if p_real >= self.cfg.rescue_threshold:
            return "PASS"
        if p_real < self.cfg.drop_threshold:
            return "FAIL"
        return "REVIEW"


def _f1_at(scores, labels, t):
    keep = scores >= t
    tp = float(labels[keep].sum()); fp = float((1 - labels[keep]).sum()); fn = float(labels[~keep].sum())
    d = 2 * tp + fp + fn
    return (2 * tp / d) if d else 0.0


def _best_threshold(scores, labels):
    """Rescue threshold on REVIEW scores that maximizes F1 (rescue real / drop artifact)."""
    if len(scores) == 0:
        return 0.5
    ts = np.unique(np.round(scores, 3))
    if len(ts) > 80:
        ts = np.quantile(scores, np.linspace(0.02, 0.98, 80))
    return float(max(ts, key=lambda t: _f1_at(scores, labels, t)))


def train_validate_by_chrom(labeled: pd.DataFrame, cfg: Optional[TriageConfig] = None) -> dict:
    """Hold out whole chromosomes (no positional leakage). The rescue threshold is TUNED on the
    TRAINING REVIEW variants each fold and applied to the held-out REVIEW variants — an honest
    deployed estimate (no threshold-peeking on the test set). Money metric: does the tuned model
    beat keep-all-REVIEW on F1?"""
    cfg = cfg or TriageConfig()
    chroms = sorted(labeled["chrom"].unique())
    folds = [chroms[i::3] for i in range(3)] if len(chroms) >= 3 else [chroms]
    aurocs, auprcs, thr_used = [], [], []
    tp_keep = fp_keep = tp_model = fp_model = fn = 0
    for held in folds:
        te = labeled[labeled.chrom.isin(held)]
        tr = labeled[~labeled.chrom.isin(held)]
        if tr.label.nunique() < 2 or len(te) == 0:
            continue
        m = ReviewTriageModel(cfg).fit(tr)
        tr_rev = tr[tr["filter"].str.upper().str.contains("REVIEW")]
        t = (_best_threshold(m.score(tr_rev), tr_rev.label.values.astype(float))
             if len(tr_rev) and tr_rev.label.nunique() == 2 else cfg.rescue_threshold)
        thr_used.append(round(t, 3))
        p = m.score(te)
        if te.label.nunique() == 2:
            aurocs.append(roc_auc_score(te.label, p))
            auprcs.append(average_precision_score(te.label, p))
        rev = te["filter"].str.upper().str.contains("REVIEW").values
        for is_rev, lab, pr in zip(rev, te.label.values, p):
            if not is_rev:
                continue
            tp_keep += lab; fp_keep += (1 - lab)         # keep-all-REVIEW baseline
            if pr >= t:                                  # tuned model keeps as PASS
                tp_model += lab; fp_model += (1 - lab)
            else:                                        # model drops -> if it was real, that's a miss
                fn += lab
    def f1(tp, fp, fn_):
        denom = 2 * tp + fp + fn_
        return (2 * tp / denom) if denom else 0.0
    keep_f1, triaged_f1 = f1(tp_keep, fp_keep, 0), f1(tp_model, fp_model, fn)
    return {"auroc": round(float(np.mean(aurocs)), 4) if aurocs else None,
            "auprc": round(float(np.mean(auprcs)), 4) if auprcs else None,
            "review_keep_all_f1": round(keep_f1, 4),
            "review_model_triaged_f1": round(triaged_f1, 4),
            "f1_gain": round(triaged_f1 - keep_f1, 4),
            "tuned_threshold": round(float(np.mean(thr_used)), 3) if thr_used else None,
            "review_tp": int(tp_keep), "review_fp": int(fp_keep),
            "n_variants": len(labeled), "n_chroms": len(chroms)}


def triage_vcf(vcf_in: str, vcf_out: str, model: ReviewTriageModel, sample: Optional[str] = None) -> dict:
    """Rewrite a VCF reclassifying REVIEW-tier records by the model. Returns a tally."""
    if not _HAVE_PYSAM:
        raise RuntimeError("pysam required")
    feat = extract_vcf_features(vcf_in, sample)
    scores = model.score(feat) if len(feat) else np.array([])
    score_by_key = {(r.chrom, r.pos, r.ref, r.alt): s for r, s in zip(feat.itertuples(), scores)}
    vf = pysam.VariantFile(vcf_in)
    hdr = vf.header.copy()
    for filt in ("PASS", "FAIL", "REVIEW"):
        if filt not in hdr.filters:
            hdr.filters.add(filt, None, None, "triaged")
    out = pysam.VariantFile(vcf_out, ("wz" if str(vcf_out).endswith(".gz") else "w"), header=hdr)
    tally = {"rescued": 0, "dropped": 0, "kept_review": 0, "untouched": 0}
    for rec in vf:
        cur = ";".join(rec.filter.keys()) if len(rec.filter.keys()) else "PASS"
        if "REVIEW" in cur.upper():
            alt = rec.alts[0] if rec.alts else ""
            s = score_by_key.get((rec.chrom, rec.pos, rec.ref, alt), 0.5)
            new = model.triage_filter(s, cur)
            rec.filter.clear(); rec.filter.add(new)
            tally["rescued" if new == "PASS" else "dropped" if new == "FAIL" else "kept_review"] += 1
        else:
            tally["untouched"] += 1
        out.write(rec)
    out.close()
    return tally


# ----------------------------------------------------------------- demo

def _demo_labeled(n=1200, seed=0) -> pd.DataFrame:
    """Synthetic VCF-feature set: real variants have high QD/MQ, balanced AD, low strand bias;
    artifacts have low QD, skewed AD, high FS. Tags ~half as REVIEW tier. (NOT real GIAB data.)"""
    rng = np.random.default_rng(seed)
    rows = []
    chroms = [f"chr{i}" for i in range(1, 7)]
    for i in range(n):
        real = rng.random() < 0.55
        qd = rng.normal(18 if real else 6, 4)
        fs = abs(rng.normal(2 if real else 25, 8))
        mq = rng.normal(58 if real else 40, 6)
        ab = np.clip(rng.normal(0.5 if real else 0.18, 0.12), 0, 1)
        sor = abs(rng.normal(1.0 if real else 3.5, 1.0))
        qual = rng.normal(200 if real else 45, 60)
        gq = rng.normal(60 if real else 25, 15)
        is_indel = int(rng.random() < 0.18)
        rows.append({"chrom": rng.choice(chroms), "pos": int(rng.integers(1, 1_000_000)),
                     "ref": "A", "alt": "G", "filter": "REVIEW" if rng.random() < 0.5 else "PASS",
                     "QUAL": qual, "QD": qd, "FS": fs, "SOR": sor, "MQ": mq,
                     "MQRankSum": rng.normal(0, 1), "ReadPosRankSum": rng.normal(0, 1),
                     "BaseQRankSum": rng.normal(0, 1), "DP": abs(rng.normal(35, 12)),
                     "GQ": gq, "allele_balance": ab, "is_indel": is_indel,
                     "indel_len": is_indel * int(rng.integers(1, 5)), "n_alt": 1,
                     "label": int(real)})
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description="REVIEW-tier variant triage model")
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--train-vcf"); ap.add_argument("--truth-vcf")
    ap.add_argument("--score-vcf"); ap.add_argument("--out-vcf")
    ap.add_argument("--model-out"); ap.add_argument("--out")
    a = ap.parse_args()

    if a.demo:
        labeled = _demo_labeled()
        print("DEMO mode - synthetic VCF features (NOT real GIAB data)")
        metrics = train_validate_by_chrom(labeled)
        print(json.dumps(metrics, indent=2, default=str))
        gain = (metrics["review_model_triaged_f1"] or 0) - (metrics["review_keep_all_f1"] or 0)
        print(f"\nREVIEW-tier F1: keep-all={metrics['review_keep_all_f1']} -> "
              f"model-triaged={metrics['review_model_triaged_f1']}  (delta {gain:+.4f})")
        if a.out:
            Path(a.out).write_text(json.dumps(metrics, indent=2, default=str))
        return 0

    if a.train_vcf and a.truth_vcf:
        feat = extract_vcf_features(a.train_vcf)
        labeled = label_against_truth(feat, a.truth_vcf)
        print(json.dumps(train_validate_by_chrom(labeled), indent=2, default=str))
        return 0

    print("provide --demo, or --train-vcf + --truth-vcf, or --score-vcf + --out-vcf")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
