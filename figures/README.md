# Figures

Screenshots of the PGX pipeline UI and outputs. **Every figure shows synthetic data** — the
`demo_pgx_showcase` cohort (128 synthetic samples labelled `CARDIO_*`, `CONTROL_*`, `ONCO_*`).
No patient data appears in any figure.

Most of these are interactive in the real thing — try the
[live demos](../README.md#-try-it--two-live-interactive-demos).

---

### The application

| | |
|---|---|
| **[01 — Home](01-main-ui-home.png)** | The local web UI. Runs are organised into projects; everything executes on the machine hosting the server, so genomic data never leaves the building. |
| **[02 — Variant browser preview](02-variant-browser-preview.png)** | Per-run result view. Tabs for overview, variant browser, pharmacogenomics, medications, genome map, analytics, QC and provenance. |
| **[03 — Variant detail](03-variant-detail-preview.png)** | One variant, all evidence in one place: annotation, ACMG/AMP criteria, caller confirmation, clinical evidence, PGx drug guidance, population frequency, cohort statistics, source links. |
| **[04 — Genome map](04-genome-map.png)** | Cohort variant density per chromosome, binned and coloured by dominant VEP impact. |

### Analytics

| | |
|---|---|
| **[05](05-analytics-figure-gallery-a.png)** / **[06](06-analytics-figure-gallery-b.png)** — Figure gallery | The publication-ready matplotlib output: Manhattan, QQ, forest, p-value histogram, gene signal, ClinVar and impact breakdowns, sample QC, case/control stats, metabolizer activity heatmap, CNV caller support heatmap, PGx readiness heatmap, ML feature importance, ROC and precision-recall, VAE patient embedding. |
| **[07 — Interactive Manhattan](07-interactive-manhattan.png)** | Click any point and a detail card opens with that variant's rsID, sample, impact, ACMG class, p-value, odds ratio, affected medications, star allele, catalog phenotype, evidence sources and caller support. |
| **[08 — Medication evidence explorer](08-medication-evidence-explorer.png)** | Inverts the view: start from a drug, see which variants in the cohort could disrupt it and in which samples. This is the view a pharmacogenomics researcher actually starts from. |

### MegaVCF Explorer

The cohort-level output — one self-contained HTML file holding every variant from every sample.

| | |
|---|---|
| **[09 — Variant browser](09-megavcf-explorer-variant-browser.png)** | Faceted filtering by gene, rsID, position, type, PGx drug, gnomAD AF, impact, ACMG criteria, database and status. Reviewer comments and filtered CSV export are built in. |
| **[10 — ACMG evidence editor](10-variant-detail-acmg-editor.png)** | Add or remove ACMG criteria and the classification **recomputes live**. Curation is a judgement call informed by context the pipeline doesn't have, so the interface treats it that way. |
| **[11 — PGx evidence](11-variant-detail-pgx-evidence.png)** | Criteria met vs opposed, then PharmGKB / CPIC (with level) / DPWG drug guidance, PharmVar star allele and function, ClinVar and LOVD classification. |
| **[12 — Population and cohort](12-variant-detail-population-cohort.png)** | Allele frequency by population with popmax, consequence and impact, OMIM disease context, cohort p-value / odds ratio / mean depth, actionability, and every contributing evidence source. |
| **[13 — Diplotypes and phenotypes](13-diplotype-phenotype.png)** | Per-sample star-allele diplotypes with metaboliser phenotype, activity score, calling source, call status and CPIC drug guidance. |
