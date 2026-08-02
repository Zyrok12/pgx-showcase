# Results

Two independent validations of variant calling, plus an evaluation of the ML triage model.

- [Validation design](#validation-design)
- [Synthetic cohort](#synthetic-cohort)
- [GIAB HG002](#giab-hg002)
- [Triage model](#triage-model)
- [What the numbers mean](#what-the-numbers-mean)
- [Limitations](#limitations)

---

## Validation design

Building a new calling strategy is only interesting if you can show it doesn't miss variants and
doesn't invent them. Two complementary tests:

**1. Synthetic cohort with a controlled truth set.** A generator builds small PGx test cohorts —
controls plus cases carrying known SNVs and CNVs — writing FASTQs alongside truth tables
(`snv_truth.tsv`, `cnv_truth.bed`, `variant_truth.tsv`, `phenotype.tsv`). SNV reference alleles
are validated against the configured reference genome so the truth table can't silently disagree
with GRCh38. This gives exact expected-vs-detected counts for **both SNVs and CNVs**.

**2. GIAB HG002.** The NIST Genome in a Bottle consortium benchmark, run as whole-exome. Real
sequencing data, real error modes, an externally curated truth set. Calls are labelled by exact
`chrom:pos:ref:alt` match against the truth VCF, restricted to the high-confidence truth BED.

Scoring is F1 throughout:

$$F1 = \frac{2TP}{2TP + FP + FN}$$

CNVs were validated on the synthetic cohort only — no suitable truth BED for HG002 CNVs was
available.

---

## Synthetic cohort

Generated cohort: 259 SNVs and 24 CNVs across the truth set.

| | SNV / indel | CNV |
|---|---:|---:|
| Expected | 259 | 24 |
| True positives | 257 | 23 |
| **False positives** | **0** | **0** |
| False negatives | 2 | 1 |
| **F1** | **0.996** | **0.978** |

Both misses are explainable rather than random:

- The **2 missed SNVs** were inside a homozygous deletion. There is no sequence there to call
  from — these are *correct* non-calls, and counting them as misses is the conservative choice.
- The **1 missed CNV** was a heterozygous duplication. Het duplications shift read depth by only
  ~1.5× against noisy WES capture, and are the acknowledged hardest CNV class in exome data.

---

## GIAB HG002

Whole-exome, SNVs and indels.

| Caller | F1 |
|---|---:|
| GATK HaplotypeCaller | 0.902 |
| Strelka2 | 0.909 |
| DeepVariant | 0.911 |
| PGX ensemble (PASS tier only) | 0.907 |
| **PGX ensemble + ML triage** | **0.929** |

All standalone callers were run on the same data with the same reference, so the comparison is
like-for-like.

Two things are worth separating here, because conflating them would overstate the result:

- **Ensemble calling alone is not the win.** 0.907 sits right in the middle of the standalone
  callers. Combining callers improves *precision* — it is very good at not producing false
  positives — but on its own it doesn't beat DeepVariant on F1.
- **The triage model is the win.** +0.022 F1, taking the pipeline clear of every individual
  caller.

False positives remained at **0** throughout.

---

## Triage model

The [calibrated gradient-boosted tree](../code-excerpts/review_triage_model.py) re-scores the
REVIEW tier into P(real) and reclassifies: ≥0.70 → PASS, <0.30 → FAIL, otherwise stay REVIEW.

### Discrimination

| Metric | Value |
|---|---:|
| AUROC | 0.80 |
| AUPRC | 0.97 |

### Effect on calling

| Dataset | PASS tier only | + triage | Δ |
|---|---:|---:|---:|
| GIAB HG002 | 0.907 | **0.929** | **+0.022** |
| Synthetic cohort | 0.972 | **0.996** | **+0.024** |

### Reclassification volume (GIAB HG002)

| Action | Variants |
|---|---:|
| Rescued REVIEW → PASS | 1,373 |
| Dropped REVIEW → FAIL | 215 |
| Kept as REVIEW | 134 |

False positives stayed at 0 in both datasets after triage.

### Validation protocol

- **Chromosome-held-out splitting**, three folds — whole chromosomes are held out rather than
  random rows, so linkage and positional similarity can't leak between train and test. Random
  splitting here would inflate performance substantially and dishonestly.
- **The rescue threshold is tuned on the training fold** of each split and applied blind to the
  held-out fold. Tuning the threshold on the test set is the easiest way to fake a good number,
  and this design forecloses it.
- **Two deliberately naive baselines** — keep-all-REVIEW and drop-all-REVIEW — because a model is
  only useful relative to doing nothing.
- **Isotonic calibration** (`CalibratedClassifierCV`) is applied whenever there are enough
  positive and negative examples, because the decision rule is stated in probabilities. An
  uncalibrated score of 0.70 doesn't mean 70%, and thresholding on it would be arbitrary.

---

## What the numbers mean

**Hard filters throw away real biology.** GATK's own thresholds cap at 0.902 — a meaningful share
of what they discard is genuine, just borderline. That's most damaging exactly where this
pipeline is aimed: pharmacogenes sit beside pseudogenes and near-identical paralogs, so their
true variants routinely carry mediocre mapping quality. Filtering hard on MQ preferentially
deletes the variants the study exists to find.

**A tiered filter plus a learned second opinion beats a better threshold.** The interesting
comparison isn't ensemble-vs-single-caller, it's *tier-and-rescue* vs *cut-once*. The REVIEW tier
holds 1,722 GIAB calls that a hard filter would have decided about blindly; the model gets 1,373
of them right in the rescue direction, and the F1 gain follows.

**AUROC 0.80 with AUPRC 0.97 is the useful shape.** Ranking quality is only fair, but precision
among promoted variants is excellent — when the model says PASS it is almost never promoting an
artifact. Given that a false positive here can mean a wrong biological conclusion in a clinical
research context, that asymmetry is the right one to have.

**Zero false positives, twice.** This was the explicit design goal rather than a happy accident,
and it drove the consensus requirements throughout — two-caller agreement for SNVs, ≥3/5 for
CNVs, and a conservative rescue threshold.

**The obvious extrapolation.** A shallow tree ensemble on ~14 tabular features produced +0.022
F1. Richer features and more labelled genomes should go considerably further. WES variant calling
has lagged WGS precisely because there's less signal per locus, so this is where learned
approaches have the most headroom.

---

## Limitations

Stated plainly, because the numbers above are validation evidence and not clinical evidence:

- **No clinical validation.** This is research-use-only software. Establishing diagnostic
  accuracy requires a formal clinical validation study that has not been performed.
- **CNVs were validated on synthetic data only.** No GIAB CNV truth BED was available for the
  exome. The synthetic F1 of 0.978 reflects a generator whose noise model is necessarily simpler
  than real capture data.
- **The synthetic cohort is small.** 259 SNVs and 24 CNVs. It demonstrates correctness, not
  robustness at scale.
- **GIAB is one genome.** HG002 is one individual of one ancestry. Callers — and learned models
  especially — behave differently across populations, and a model trained on one genome's error
  modes may not transfer.
- **The triage model inherits its training data's biases.** Labels come from GIAB, so it has
  learned what artifacts look like *in that data*.
- **WES CNV breakpoints are estimates.** Read-depth calling cannot resolve breakpoints that fall
  in unsequenced introns. Single-exon and single-target events can be real but need orthogonal
  confirmation.
- **Small cohorts weaken CNV reference structure.** Callers building a panel of normals from the
  cohort itself degrade when the cohort is small.
- **Association results are research signals.** Ancestry structure, capture bias, missingness and
  multiple testing all have to be weighed before interpretation.
- **Star-allele calling is limited** by WES capture, paralogy, phasing and structural variation.

Any clinical-grade conclusion requires orthogonal validation and expert sign-out.
