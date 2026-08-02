# Code excerpts

Two components from the PGX pipeline, both self-contained and runnable. They're here because
they're the parts where the interesting decisions live — everything else is orchestration,
I/O and format wrangling.

Requires Python 3.10+, `numpy`, `pandas`, `scikit-learn`.

---

## `review_triage_model.py`

**Verbatim** from `pipeline/review_triage_dl.py` — the calibrated gradient-boosted tree that
re-scores the REVIEW quality tier. This is the component responsible for the +0.022 F1 gain on
GIAB HG002 ([results](../docs/RESULTS.md#triage-model)).

```bash
python review_triage_model.py --demo
```

The only change from the production file is replacing two non-ASCII characters in `print`
statements, which crash on Windows consoles under cp1252.

**What to look at:**

- `ReviewTriageModel.fit` — shallow boosted trees inside `CalibratedClassifierCV(method="isotonic")`.
  Calibration isn't decoration: the decision rule is stated in probabilities (rescue at ≥0.70,
  drop below 0.30), so the scores have to actually be probabilities.
- `train_validate_by_chrom` — validation splits by **whole chromosome**, not randomly, so nearby
  variants can't leak between train and test. The rescue threshold is re-tuned on the training
  fold each split and applied blind to the held-out fold, which forecloses threshold-peeking.
  Results are reported against two deliberately naive baselines (keep-all-REVIEW,
  drop-all-REVIEW).
- `extract_vcf_features` — allele balance derived from AD is the feature carrying the most signal:
  real heterozygotes sit near 0.5, artifacts are skewed.
- `triage_filter` — the reclassification rule, three lines, deliberately boring.

> ⚠️ `--demo` runs the full train/validate loop on **synthetic, easily separable** features to
> prove the plumbing end to end, so it reports AUROC ≈ 1.0. That number means nothing. The real
> measurements against GIAB HG002 are in [RESULTS.md](../docs/RESULTS.md).

---

## `cnv_ensemble_voting.py`

**Distilled** from `pipeline/cnv_ensemble.py` (~140 KB). A faithful, runnable implementation of
the voting core: merge → cluster → vote → assign confidence → report breakpoint uncertainty.

```bash
python cnv_ensemble_voting.py --demo
```

The production module also handles per-caller I/O and format quirks, gene/target/exon annotation,
mappability cross-referencing, population database annotation (gnomAD-CNV / DGV / DECIPHER),
control-sample handling and CNV-VEP integration. The algorithm is what's reproduced here.

**What to look at:**

- `merge_within_caller` — the step that's easy to skip and breaks everything if you do. A
  250 bp-bin caller emits ~20 adjacent calls for one biological CNV; without merging first it
  casts 20 votes, drowns out the target-based callers, and drags reported breakpoints onto its
  own bin grid. The demo makes this concrete: 74 raw calls → 11 merged events.
- `caller_spans` — enforces one caller, one vote.
- `build_consensus` — breakpoints reported as a **median across per-caller spans**, alongside the
  target-confirmed core (intersection) and outer supported span (union). WES can't resolve
  breakpoints that fall in unsequenced introns, so the output states its uncertainty instead of
  picking a number and implying precision it doesn't have.
- `_apply_exclusivity` — an overlapping deletion and duplication in one sample can't both be real;
  the better-supported call wins and the other is demoted rather than silently deleted.

The demo output shows all three regimes:

```
CYP2D6   chr22:42522250-42527750  DEL  support=5/5  consensus  filter=PASS
TPMT     chr6:18130250-18142750   DEL  support=4/5  review     filter=REVIEW
DPYD     chr1:97450500-97461500   DUP  support=2/5  candidate  filter=LOW_SUPPORT
```

---

## What isn't here

The full pipeline is ~40 modules: stage orchestration with checkpoint/resume, the four entry
points, VEP integration with a persistent hit cache, seven database clients with merge and
novelty tagging, ACMG criteria evaluation, star-allele and diplotype calling, the cohort analytics
suite, the ML responder workbench, the MegaVCF builder and its self-contained browser app, the
local web UI and server, preflight validation, and run manifesting.

That code is private. [ARCHITECTURE.md](../docs/ARCHITECTURE.md) describes how it fits together,
and the [live demos](../README.md#-try-it--two-live-interactive-demos) show what it produces.

Happy to walk through any of it in conversation.
