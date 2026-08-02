#!/usr/bin/env python3
"""
cnv_ensemble_voting.py — the five-caller WES CNV consensus algorithm.

DISTILLED EXCERPT. This is a faithful, self-contained, runnable implementation of the voting
core from `pipeline/cnv_ensemble.py`. The production module additionally handles caller I/O and
format quirks, gene/target/exon annotation, mappability cross-referencing, population database
annotation (gnomAD-CNV / DGV / DECIPHER), control-sample handling, and CNV-VEP integration. The
algorithm below — merge, cluster, vote, assign confidence, report breakpoint uncertainty — is
the part worth reading.

WHY IT LOOKS LIKE THIS
----------------------
WES CNV calling is read-depth based: you only observe captured exons, so a deletion whose true
breakpoints fall in introns has invisible edges. Individual callers disagree constantly. Five
callers are used, deliberately mixing two resolutions so their failure modes don't correlate:

    bin-based    (250 bp bins)   ExomeDepth, CODEX2, panelcn.MOPS
    target-based (per-exon/gene) CNVkit, GATK-gCNV

Step 1 (within-caller merge) is the step people skip, and skipping it breaks everything: a
250 bp-bin caller emits ~20 adjacent calls for one biological CNV. Without merging first, that
caller casts 20 votes, drowns out the target-based callers, and drags the reported breakpoints
toward its own bin grid.

Run the worked example:

    python cnv_ensemble_voting.py --demo
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from typing import Iterable

# The five WES read-depth callers in the ensemble.
CALLERS = ("CNVkit", "ExomeDepth", "GATK-gCNV", "CODEX2", "panelcn.MOPS")

# Confidence tiers by caller support. Below MIN_CALLERS a cluster is retained as raw
# evidence but never emitted as a call.
CONSENSUS_CALLERS = 5
MIN_CALLERS = 3

# Adjacent same-caller calls closer than this are treated as one fragmented event.
MERGE_GAP_BP = 1000

# Reciprocal-overlap fraction required to cluster two calls that share no gene or target.
CLUSTER_RECIPROCAL_OVERLAP = 0.50


# --------------------------------------------------------------------- helpers

def _chrom_norm(chrom: str) -> str:
    return str(chrom).replace("chr", "").strip().upper()


def _direction(cnv_class: str) -> str:
    """Collapse caller-specific class names onto gain / loss."""
    text = str(cnv_class).lower()
    if any(k in text for k in ("del", "loss")):
        return "loss"
    if any(k in text for k in ("dup", "gain", "amp")):
        return "gain"
    return "neutral"


def _genes(call: dict) -> set[str]:
    raw = str(call.get("gene", "") or "")
    return {g.strip() for g in raw.replace(",", ";").split(";") if g.strip() and g.strip() != "."}


def _overlap_bp(a: dict, b: dict) -> int:
    return max(0, min(int(a["end"]), int(b["end"])) - max(int(a["start"]), int(b["start"])))


def _reciprocal_overlap(a: dict, b: dict) -> float:
    """Overlap as a fraction of the SHORTER interval, so a huge call can't swallow a small one."""
    shared = _overlap_bp(a, b)
    if shared <= 0:
        return 0.0
    len_a = int(a["end"]) - int(a["start"])
    len_b = int(b["end"]) - int(b["start"])
    shorter = min(len_a, len_b)
    return shared / shorter if shorter > 0 else 0.0


def _median_int(values: Iterable[int]) -> int:
    vals = sorted(int(v) for v in values)
    if not vals:
        return 0
    mid = len(vals) // 2
    if len(vals) % 2:
        return vals[mid]
    return int(round((vals[mid - 1] + vals[mid]) / 2))


def _majority(values: Iterable[str]) -> str:
    counts = Counter(v for v in values if v)
    return counts.most_common(1)[0][0] if counts else "."


def _mean(values: list[float]) -> float | None:
    vals = [v for v in values if v is not None]
    return sum(vals) / len(vals) if vals else None


# --------------------------------------------------- step 1: within-caller merge

def merge_within_caller(calls: list[dict]) -> list[dict]:
    """
    Collapse fragmented adjacent bins from the SAME caller into one event, BEFORE voting.

    Two calls merge when they share sample, caller, chromosome and direction, are not
    gene-disjoint, and either overlap or sit within MERGE_GAP_BP of each other.
    """
    merged: list[dict] = []
    block: list[dict] = []

    ordered = sorted(
        calls,
        key=lambda c: (c["sample"], c["caller"], _chrom_norm(c["chrom"]),
                       _direction(c["cnv_class"]), int(c["start"]), int(c["end"])),
    )

    for call in ordered:
        if not block:
            block = [call]
            continue
        if _mergeable(block[-1], call):
            block.append(call)
            continue
        merged.append(_merge_block(block))
        block = [call]

    if block:
        merged.append(_merge_block(block))
    return merged


def _mergeable(a: dict, b: dict) -> bool:
    if a["sample"] != b["sample"] or a["caller"] != b["caller"]:
        return False
    if _chrom_norm(a["chrom"]) != _chrom_norm(b["chrom"]):
        return False
    if _direction(a["cnv_class"]) != _direction(b["cnv_class"]):
        return False
    genes_a, genes_b = _genes(a), _genes(b)
    if genes_a and genes_b and not (genes_a & genes_b):
        return False
    return _overlap_bp(a, b) > 0 or (int(b["start"]) - int(a["end"])) <= MERGE_GAP_BP


def _merge_block(block: list[dict]) -> dict:
    """Weighted-merge a run of same-caller calls. Depth estimates weight by supporting bins."""
    if len(block) == 1:
        out = dict(block[0])
        out["merged_bins"] = 1
        return out

    weights = [max(1.0, float(c.get("n_points") or 1.0)) for c in block]
    total_w = sum(weights)

    def weighted(field: str) -> float | None:
        pairs = [(c.get(field), w) for c, w in zip(block, weights) if c.get(field) is not None]
        if not pairs or total_w <= 0:
            return None
        return sum(float(v) * w for v, w in pairs) / sum(w for _v, w in pairs)

    out = dict(block[0])
    out.update({
        "start": min(int(c["start"]) for c in block),
        "end": max(int(c["end"]) for c in block),
        "cnv_class": _majority([c["cnv_class"] for c in block]),
        "log2cr": weighted("log2cr"),
        "est_copy_number": weighted("est_copy_number"),
        "n_points": sum(float(c.get("n_points") or 0) for c in block),
        "gene": ";".join(sorted({g for c in block for g in _genes(c)})) or ".",
        "merged_bins": len(block),
    })
    return out


# ------------------------------------------------ step 2: cross-caller clustering

def cluster_across_callers(calls: list[dict]) -> list[list[dict]]:
    """
    Group calls that describe the same biological event.

    Bucketed by (sample, chromosome, direction) so unrelated events never compare. Within a
    bucket, a call joins a cluster if it shares a gene with it, or meets the reciprocal-overlap
    threshold.
    """
    clusters: list[list[dict]] = []
    buckets: dict[tuple[str, str, str], list[list[dict]]] = {}

    for call in sorted(calls, key=lambda c: (c["sample"], _chrom_norm(c["chrom"]),
                                             int(c["start"]), int(c["end"]))):
        key = (str(call["sample"]), _chrom_norm(call["chrom"]), _direction(call["cnv_class"]))
        bucket = buckets.setdefault(key, [])
        for cluster in bucket:
            if _same_cluster(cluster[0], call):
                cluster.append(call)
                break
        else:
            cluster = [call]
            clusters.append(cluster)
            bucket.append(cluster)

    return clusters


def _same_cluster(a: dict, b: dict) -> bool:
    if a["sample"] != b["sample"]:
        return False
    if _chrom_norm(a["chrom"]) != _chrom_norm(b["chrom"]):
        return False
    if _direction(a["cnv_class"]) != _direction(b["cnv_class"]):
        return False
    genes_a, genes_b = _genes(a), _genes(b)
    if genes_a and genes_b and (genes_a & genes_b):
        return True
    return _reciprocal_overlap(a, b) >= CLUSTER_RECIPROCAL_OVERLAP


# ------------------------------------------------------ step 3: one vote per caller

def caller_spans(cluster: list[dict]) -> dict[str, dict[str, int]]:
    """
    Collapse every call from one caller into ONE span — one caller, one vote.

    This is what stops a high-resolution caller from both out-voting the others and dragging
    the reported breakpoints onto its own bin grid.
    """
    spans: dict[str, dict[str, int]] = {}
    for call in cluster:
        caller = str(call.get("caller", "")).strip()
        if not caller:
            continue
        start, end = int(call["start"]), int(call["end"])
        if end < start:
            start, end = end, start
        if caller not in spans:
            spans[caller] = {"start": start, "end": end}
        else:
            spans[caller]["start"] = min(spans[caller]["start"], start)
            spans[caller]["end"] = max(spans[caller]["end"], end)
    return spans


# ------------------------------------------- step 4: consensus rows + breakpoint honesty

def build_consensus(calls: list[dict], min_callers: int = MIN_CALLERS) -> list[dict]:
    """
    Full pipeline: merge -> cluster -> vote -> assign confidence.

    Returns one row per biological event. Breakpoints are reported as a MEDIAN across per-caller
    spans, alongside the target-confirmed core (the region every caller agrees on) and the outer
    supported span (the union). WES cannot resolve exact breakpoints; the output says so instead
    of picking one number and implying precision it does not have.
    """
    merged = merge_within_caller([c for c in calls if _direction(c["cnv_class"]) != "neutral"])
    rows: list[dict] = []

    for cluster in cluster_across_callers(merged):
        spans = caller_spans(cluster)
        if not spans:
            continue

        supporting = sorted(spans)
        support = len(supporting)

        starts = [s["start"] for s in spans.values()]
        ends = [s["end"] for s in spans.values()]

        outer_start, outer_end = min(starts), max(ends)   # union: widest possible event
        core_start, core_end = max(starts), min(ends)     # intersection: every caller agrees
        start, end = _median_int(starts), _median_int(ends)

        # A degenerate median (callers barely overlapping) falls back to core, then to union.
        if end <= start:
            start, end = (core_start, core_end) if core_end > core_start else (outer_start, outer_end)

        if support >= CONSENSUS_CALLERS:
            confidence, ensemble_filter = "consensus", "PASS"
        elif support >= min_callers:
            confidence, ensemble_filter = "review", "REVIEW"
        else:
            confidence, ensemble_filter = "candidate", "LOW_SUPPORT"

        rows.append({
            "sample": cluster[0]["sample"],
            "chrom": cluster[0]["chrom"],
            "start": start,
            "end": end,
            "length_bp": abs(end - start),
            "cnv_class": _majority([c["cnv_class"] for c in cluster]),
            "gene": ";".join(sorted({g for c in cluster for g in _genes(c)})) or ".",
            "caller_support": support,
            "callers": ";".join(supporting),
            "confidence": confidence,
            "ensemble_filter": ensemble_filter,
            # --- breakpoint uncertainty, kept explicit ---
            "confirmed_start": core_start if core_end > core_start else "",
            "confirmed_end": core_end if core_end > core_start else "",
            "outer_start": outer_start,
            "outer_end": outer_end,
            "breakpoint_precision": _breakpoint_precision(core_start, core_end,
                                                          outer_start, outer_end),
            "est_copy_number": _mean([c.get("est_copy_number") for c in cluster]),
            "log2cr": _mean([c.get("log2cr") for c in cluster]),
            "source_interval_count": len(cluster),
        })

    return _apply_exclusivity(rows)


def _breakpoint_precision(core_start: int, core_end: int,
                          outer_start: int, outer_end: int) -> str:
    """How much do the callers actually agree on where this event begins and ends?"""
    if core_end <= core_start:
        return "low"
    outer_len = outer_end - outer_start
    core_len = core_end - core_start
    if outer_len <= 0:
        return "low"
    agreement = core_len / outer_len
    if agreement >= 0.90:
        return "high"
    if agreement >= 0.60:
        return "medium"
    return "low"


def _apply_exclusivity(rows: list[dict]) -> list[dict]:
    """
    A deletion and a duplication overlapping in the same sample cannot both be real.

    When they conflict, the better-supported call keeps PASS and the other is demoted. Rank is
    caller support first, then supporting bins, then effect size.
    """
    rows = [dict(r) for r in rows]
    for i, row in enumerate(rows):
        if row["ensemble_filter"] != "PASS":
            continue
        for j, other in enumerate(rows):
            if i == j or other["ensemble_filter"] != "PASS":
                continue
            if not _conflicts(row, other):
                continue
            if _rank(row) < _rank(other):
                row["ensemble_filter"] = "REVIEW;CONFLICT"
                row["confidence"] = "candidate"
                break
    return rows


def _conflicts(a: dict, b: dict) -> bool:
    if a["sample"] != b["sample"] or _chrom_norm(a["chrom"]) != _chrom_norm(b["chrom"]):
        return False
    if _direction(a["cnv_class"]) == _direction(b["cnv_class"]):
        return False
    if _overlap_bp(a, b) <= 0:
        return False
    genes_a, genes_b = _genes(a), _genes(b)
    if genes_a and genes_b and (genes_a & genes_b):
        return True
    return _reciprocal_overlap(a, b) >= 0.35


def _rank(row: dict) -> tuple:
    return (int(row.get("caller_support") or 0),
            int(row.get("source_interval_count") or 0),
            abs(float(row.get("log2cr") or 0.0)))


# ------------------------------------------------------------------------- demo

def _demo_calls() -> list[dict]:
    """
    Three scenarios in one synthetic sample:

      CYP2D6  a real deletion. All five callers see it, but the three bin-based callers emit
              fragmented 250 bp bins -> without within-caller merging they would cast 12 votes
              between them instead of 3.
      DPYD    a duplication seen by only two callers -> below threshold, not emitted as a call.
      TPMT    a deletion seen by four callers -> REVIEW, real enough to look at, not consensus.
    """
    calls: list[dict] = []

    # CYP2D6 deletion, chr22 - target-based callers give one call each
    for caller, (s, e) in {"CNVkit": (42_522_000, 42_528_000),
                           "GATK-gCNV": (42_522_500, 42_527_500)}.items():
        calls.append({"sample": "CARDIO_019", "caller": caller, "chrom": "chr22",
                      "start": s, "end": e, "cnv_class": "DEL", "gene": "CYP2D6",
                      "log2cr": -0.95, "est_copy_number": 1.0, "n_points": 8})

    # ...and bin-based callers fragment it into 250 bp bins
    for caller in ("ExomeDepth", "CODEX2", "panelcn.MOPS"):
        for bin_start in range(42_522_250, 42_527_750, 250):
            calls.append({"sample": "CARDIO_019", "caller": caller, "chrom": "chr22",
                          "start": bin_start, "end": bin_start + 250, "cnv_class": "DEL",
                          "gene": "CYP2D6", "log2cr": -0.90, "est_copy_number": 1.0,
                          "n_points": 1})

    # DPYD duplication - only two callers agree, should be dropped
    for caller, (s, e) in {"CNVkit": (97_450_000, 97_462_000),
                           "ExomeDepth": (97_451_000, 97_461_000)}.items():
        calls.append({"sample": "CARDIO_019", "caller": caller, "chrom": "chr1",
                      "start": s, "end": e, "cnv_class": "DUP", "gene": "DPYD",
                      "log2cr": 0.55, "est_copy_number": 3.0, "n_points": 5})

    # TPMT deletion - four callers, should land in REVIEW
    for caller, (s, e) in {"CNVkit": (18_130_000, 18_143_000),
                           "GATK-gCNV": (18_131_000, 18_142_000),
                           "CODEX2": (18_130_500, 18_144_000),
                           "panelcn.MOPS": (18_129_000, 18_142_500)}.items():
        calls.append({"sample": "CARDIO_019", "caller": caller, "chrom": "chr6",
                      "start": s, "end": e, "cnv_class": "DEL", "gene": "TPMT",
                      "log2cr": -0.70, "est_copy_number": 1.0, "n_points": 6})

    return calls


def main() -> int:
    ap = argparse.ArgumentParser(description="Five-caller WES CNV consensus")
    ap.add_argument("--demo", action="store_true", help="run the worked example")
    ap.add_argument("--min-callers", type=int, default=MIN_CALLERS)
    args = ap.parse_args()

    if not args.demo:
        ap.print_help()
        return 2

    raw = _demo_calls()
    merged = merge_within_caller(raw)

    print(f"raw caller calls             : {len(raw)}")
    print(f"after within-caller merge    : {len(merged)}")
    print("  (the 3 bin-based callers collapse from "
          f"{sum(1 for c in raw if c['caller'] in ('ExomeDepth', 'CODEX2', 'panelcn.MOPS'))} "
          "bins to 3 votes on CYP2D6)\n")

    rows = build_consensus(raw, min_callers=args.min_callers)
    for row in sorted(rows, key=lambda r: -r["caller_support"]):
        print(f"{row['gene']:<8} {row['chrom']}:{row['start']}-{row['end']}  "
              f"{row['cnv_class']:<4} support={row['caller_support']}/5  "
              f"{row['confidence']:<10} filter={row['ensemble_filter']:<16} "
              f"breakpoints={row['breakpoint_precision']}")

    dropped = [r for r in rows if r["caller_support"] < args.min_callers]
    print(f"\nemitted as calls : {len([r for r in rows if r['caller_support'] >= args.min_callers])}")
    print(f"below threshold  : {len(dropped)} (retained as raw evidence, not reported as CNVs)")
    print("\nfull rows:")
    print(json.dumps(rows, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
