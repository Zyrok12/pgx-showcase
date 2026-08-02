# Methods

Scientific background and the analysis methods the pipeline implements. For the engineering
design see [ARCHITECTURE.md](ARCHITECTURE.md); for validation numbers see [RESULTS.md](RESULTS.md).

- [Background: why pharmacogenomics needs its own pipeline](#background-why-pharmacogenomics-needs-its-own-pipeline)
- [Star-allele nomenclature](#star-allele-nomenclature)
- [ACMG/AMP classification](#acmgamp-classification)
- [Calling metrics](#calling-metrics)
- [The triage model](#the-triage-model)
- [Cohort analytics](#cohort-analytics)
- [Visualisation](#visualisation)
- [References](#references)

---

## Background: why pharmacogenomics needs its own pipeline

Pharmacogenomics aims to match drug and dose to a patient's genotype. Specific variants make
individuals ultrarapid or poor metabolisers of particular drugs, and a large share of therapeutic
failures and adverse drug reactions trace back to giving a standard dose to a non-standard
metaboliser. Roughly 95% of people carry at least one high-risk pharmacogenomic variant. The
clinical footprint is broad: oncology, psychiatry, cardiology, pain management.

Two properties make pharmacogenes awkward for a generic germline pipeline:

1. **They are heavily multi-allelic.** Sites routinely carry three or more alleles, so
   multi-allelic records must be split into independent normalised records before any catalogue
   lookup will match.
2. **They sit in hard genomic neighbourhoods.** Pharmacogenes are flanked by pseudogenes and
   near-identical paralogs (CYP2D6/CYP2D7 being the canonical case). Reads map ambiguously,
   mapping quality drops, and generic quality filters preferentially delete true variants in
   exactly the genes under study.

Both drive design decisions elsewhere in the pipeline: `bcftools norm -m -both` before enrichment,
a dedicated caller pair (DeepVariant + Octopus) for pharmacogene regions, and a tiered filter with
a learned rescue step instead of a single hard cut.

## Star-allele nomenclature

Pharmacogenomics describes genotypes as **haplotypes** — combinations of sequence variants — using
star-allele notation: `CYP2D6*1/*1` (homozygous), `CYP2D6*1/*4` (heterozygous). `*1` is the
reference allele assigned when no other variant is detected.

Two consequences that matter for implementation:

- Only variants altering the protein sequence (amino-acid change or length change) define a star
  allele. A star allele can contain many variants but is named for its functional consequence.
- **Different variants can define the same star allele** — two variants producing an identical
  consequence, e.g. the same premature stop, map to the same allele.

So star-allele calling is not a per-variant lookup; it's haplotype resolution followed by a
catalogue match. Each diplotype maps to an activity score and a metaboliser phenotype, which is
what CPIC guidelines are keyed on. The pipeline calls diplotypes with PyPGx (Cyrius for CYP2D6
where available) and reports diplotype, phenotype, activity score, calling source and status.

## ACMG/AMP classification

Variant impact follows the ACMG/AMP framework (Richards et al. 2015): **Pathogenic, Likely
Pathogenic, VUS, Likely Benign, Benign**, assigned by combining evidence codes drawn from
population, computational, functional and segregation data.

| Direction | Strength tiers |
|---|---|
| Pathogenic | Very Strong (PVS1) · Strong (PS1–4) · Moderate (PM1–6) · Supporting (PP1–5) |
| Benign | Stand-Alone (BA1) · Strong (BS1–4) · Supporting (BP1–7) |

Classification is retrieved from GeneBe where evidence isn't directly derivable, and — importantly
for a research tool — the MegaVCF Explorer lets a reviewer **add or remove criteria and see the
classification recompute live**. Curation is a judgement call informed by context the pipeline
doesn't have; the interface treats it that way rather than presenting a fixed verdict.

## Calling metrics

The tiering thresholds ([table](ARCHITECTURE.md#quality-tiering)) act on standard GATK annotations.
What each one is actually detecting:

| Metric | What it measures | Failure mode it catches |
|---|---|---|
| **QD** — Quality by Depth | variant confidence normalised by depth | high QUAL driven by extreme depth rather than per-base accuracy |
| **FS** — Fisher Strand | Phred-scaled p-value for strand bias | alt allele seen mostly on one strand → likely artifact |
| **SOR** — Strand Odds Ratio | strand bias as an odds ratio | preferred over FS at high depth, where Fisher's test over-triggers |
| **MQ** — Mapping Quality | RMS mapping quality of supporting reads | reads that couldn't be placed confidently |
| **MQRankSum** | alt vs ref mapping-quality comparison | negative → alt-supporting reads map worse than ref-supporting ones |
| **ReadPosRankSum** | position of the allele within reads | alt allele clustered at read ends, where error rates spike |

Normalising QUAL by depth (QD) rather than thresholding QUAL directly is what keeps deep,
low-quality pileups from passing.

## The triage model

**Problem.** PASS/REVIEW/FAIL tiering is a better shape than a single cut, but REVIEW is a mixed
bag: real variants that are borderline because of mapping complexity, low-MQ regions or uneven
WES capture, sitting alongside genuine artifacts. No fixed threshold separates them cleanly,
because the same QD value means different things at different MQ, depth and allele balance.

**Approach.** Learn the decision. Score each REVIEW call as P(real) and reclassify on that.

**Features** (14, all from standard VCF annotations — no new computation required):

```
QUAL, QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum, BaseQRankSum, DP, GQ,
allele_balance (alt fraction from AD), is_indel, indel_len, n_alt
```

Allele balance is doing a lot of work: real heterozygotes sit near 0.5 and homozygotes near 1.0,
while artifacts are skewed. Indel status and length are included because indels carry different
quality distributions and would otherwise pull the thresholds around.

**Labels are free.** Every call in a GIAB run is TP/FP-labelable by exact `chrom:pos:ref:alt`
match against the NIST truth set — label = 1 if present in the truth VCF and inside the
high-confidence BED, 0 otherwise. A large supervised training set with no manual curation.

**Model.**

```python
base = HistGradientBoostingClassifier(
    max_depth=3, learning_rate=0.05, max_iter=300, l2_regularization=1.0
)
model = CalibratedClassifierCV(base, method="isotonic", cv=3)
```

Shallow boosted trees fit the problem: tabular features, non-linear interactions between quality
metrics, and a sequential fit where each tree corrects its predecessors' residuals. `max_depth=3`
plus L2 regularisation keeps individual trees from memorising the training genome.

**Calibration is not optional here.** The decision rule is stated in probabilities — rescue at
P ≥ 0.70, drop below 0.30 — so the scores have to *be* probabilities. Isotonic regression fits a
monotonic mapping from raw score to empirical frequency on held-out validation folds, applied
whenever there are enough examples of each class.

**Decision rule.**

| P(real) | Action |
|---|---|
| ≥ 0.70 | rescue REVIEW → **PASS** |
| 0.30 – 0.70 | keep as **REVIEW** |
| < 0.30 | drop REVIEW → **FAIL** |

Calls not in the REVIEW tier are never touched.

**Validation** is by chromosome-held-out splitting with the threshold re-tuned per training fold —
see [RESULTS.md](RESULTS.md#validation-protocol) for why, and
[`review_triage_model.py`](../code-excerpts/review_triage_model.py) for the implementation.

## Cohort analytics

Four families, all cohort-level.

### Population structure

- **PCA** on the genotype dosage matrix (`sklearn` `StandardScaler` → decomposition), to detect
  ancestry stratification that would otherwise produce spurious associations.
- **Identity-by-state (IBS)** — mean proportion of variants at identical dosage between each pair
  of samples, rendered as a heatmap. Unexpectedly high IBS flags duplicates, relatives or
  mislabelled samples.
- **Runs of homozygosity (ROH)** — extended homozygous stretches, indicating consanguinity,
  uniparental disomy or hemizygous deletion. Reported as per-sample ROH burden.

### Quality control

- Per-sample heterozygosity rate, call rate, mean depth, duplication rate, **Ti/Tv ratio**
  (transitions vs transversions; healthy WES sits near 2.0 — a deviation signals a systematic
  calling problem before any biology is interpreted).
- Coverage: mean and median depth, breadth at 1×/10×/20×, mean-to-median ratio. Uneven coverage is
  a common source of false association.
- **Hardy–Weinberg equilibrium** per variant, observed vs expected genotype counts.
- **Allele frequency spectrum** as a MAF histogram, showing the rare-variant load.
- **Variant density** — counts per chromosome, and the top 20 genes by variant count.

### Association

**PLINK2 logistic regression with Firth penalisation** is the primary engine, testing responders
vs non-responders (or any case/control contrast: toxicity, adverse event, dose requirement).
Firth penalisation matters because pharmacogenomic cohorts are small and carrier counts are often
tiny — standard logistic regression separates completely and produces infinite odds ratios
exactly where the interesting signal is.

With complete sex metadata, autosomes 1–22 plus chrX/PAR are tested with GRCh38 PAR handling;
chrY and MT are excluded by design. Outputs: p-values (histogram), Bonferroni threshold (Manhattan
plot), odds ratios with 95% CI (forest plot), Benjamini–Hochberg q-values, gene-signal plot, and a
case/control carrier summary. Fisher exact carrier and allele tests are written as a companion
table but never drive the primary plots.

CNVs enter the genotype matrix as variant-like carrier events, so they're tested alongside small
variants:

| CNV state | Dosage |
|---|---:|
| no event | 0 |
| heterozygous deletion/duplication | 1 |
| homozygous deletion / high amplification | 2 |

### LD and PGx readiness

- **Linkage disequilibrium** as pairwise r². A cohort can carry >50,000 variants, so the heatmap
  is restricted to the top-ranked significant variants; full pairwise r² is written to CSV so any
  variant's LD partners can be looked up.
- **PGx readiness heatmap** — per-gene coverage adequacy across the pharmacogene panel, answering
  "can I actually trust the calls in this gene for this cohort?" before interpretation starts.

### ML responder workbench

An exploratory layer separate from the primary association analysis: regularised logistic
baselines, random forest, a shallow neural classifier, calibration curves, VAE-style patient
embeddings (with a deterministic linear fallback when PyTorch is absent), and **permutation sanity
checks** to confirm that reported performance collapses when labels are shuffled.

## Visualisation

- **Plotly** for interactive figures — clicking a point on the Manhattan plot opens a detail card
  with that variant's PGx, ACMG, caller and medication evidence.
- **Matplotlib** for static, publication-ready figures.
- **SciPy** for the underlying statistics.

Both are produced from the same run, so the figure in a paper and the figure a researcher clicked
through are the same analysis.

## References

Core literature behind the methods above.

1. McLaren W, *et al.* **The Ensembl Variant Effect Predictor.** *Genome Biology* 17:122 (2016).
   [doi:10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4)
2. Richards S, *et al.* **Standards and guidelines for the interpretation of sequence variants:
   a joint consensus recommendation of ACMG and AMP.** *Genetics in Medicine* 17(5):405–424 (2015).
   [doi:10.1038/gim.2015.30](https://doi.org/10.1038/gim.2015.30)
3. Tayeh MK, *et al.* **Clinical pharmacogenomic testing and reporting: a technical standard of
   the ACMG.** *Genetics in Medicine* (2022).
   [doi:10.1016/j.gim.2021.12.009](https://doi.org/10.1016/j.gim.2021.12.009)
4. Whirl-Carrillo M, *et al.* **An evidence-based framework for evaluating pharmacogenomics
   knowledge for personalized medicine.** *Clinical Pharmacology & Therapeutics* 110(3):563–572 (2021).
5. Donnelly RS, *et al.* **Decoding pharmacogenomic test interpretation and application to patient
   care.** *JACCP* 7(6):581–588 (2024). [doi:10.1002/jac5.1958](https://doi.org/10.1002/jac5.1958)
6. Van der Auwera GA, *et al.* **From FastQ data to high confidence variant calls: the GATK best
   practices pipeline.** *Current Protocols in Bioinformatics* 43:11.10.1–33 (2013).
   [doi:10.1002/0471250953.bi1110s43](https://doi.org/10.1002/0471250953.bi1110s43)
7. Kim S, *et al.* **Strelka2: fast and accurate calling of germline and somatic variants.**
   *Nature Methods* 15:591–594 (2018). [doi:10.1038/s41592-018-0051-x](https://doi.org/10.1038/s41592-018-0051-x)
8. Poplin R, *et al.* **A universal SNP and small-indel variant caller using deep neural
   networks.** *Nature Biotechnology* 36:983–987 (2018). [doi:10.1038/nbt.4235](https://doi.org/10.1038/nbt.4235)
9. Zook J, *et al.* **Extensive sequencing of seven human genomes to characterize benchmark
   reference materials.** *Scientific Data* 3:160025 (2016). [doi:10.1038/sdata.2016.25](https://doi.org/10.1038/sdata.2016.25)
10. Danecek P, *et al.* **Twelve years of SAMtools and BCFtools.** *GigaScience* 10(2):giab008 (2021).
    [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)
11. Abdelwahab O, Torkamaneh D. **Artificial intelligence in variant calling: a review.**
    *Frontiers in Bioinformatics* 5:1574359 (2025). [doi:10.3389/fbinf.2025.1574359](https://doi.org/10.3389/fbinf.2025.1574359)
12. Chen X, Wang M, Zhang H. **The use of classification trees for bioinformatics.**
    *WIREs Data Mining and Knowledge Discovery* 1(1):55–63 (2011). [doi:10.1002/widm.14](https://doi.org/10.1002/widm.14)
13. Vaulet T, *et al.* **Gradient boosted trees with individual explanations: an alternative to
    logistic regression for viability prediction in the first trimester of pregnancy.**
    *Computer Methods and Programs in Biomedicine* 213:106520 (2022). [doi:10.1016/j.cmpb.2021.106520](https://doi.org/10.1016/j.cmpb.2021.106520)
14. Riccio C, *et al.* **Variant effect predictors: a systematic review and practical guide.**
    *Human Genetics* 143(5):625–634 (2024). [doi:10.1007/s00439-024-02670-5](https://doi.org/10.1007/s00439-024-02670-5)
15. Behjati S, Tarpey PS. **What is next generation sequencing?** *Archives of Disease in
    Childhood — Education and Practice* 98(6):236–238 (2013). [doi:10.1136/archdischild-2013-304340](https://doi.org/10.1136/archdischild-2013-304340)
