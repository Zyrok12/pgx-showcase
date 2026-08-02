# Publishing this repository

Steps to get it live, plus how to refresh the demos later.

---

## 1. Scope of what is published

Intellectual property in the pipeline belongs to **Bruno Young de Castro**, so publishing this
showcase is his call.

Two things keep the published surface narrow regardless: no pipeline source is included beyond
the two excerpts in `code-excerpts/`, and no patient data appears anywhere in the repository —
every figure, demo page and table is generated from synthetic cohorts.

---

## 2. Publishing checklist

Already done in this repository:

- GitHub username, LinkedIn URL and contact email filled in throughout.
- `.gitattributes` forces LF on `*.sh`, so `refresh_demo.sh` survives a Windows checkout.
- `docs/.nojekyll` stops Pages running Jekyll over the demo HTML.
- Demo pages copied from a run marked `"synthetic_demo": true`, with local paths scrubbed.

Before any push, sanity-check what you are about to make public:

```bash
git status --short
du -sh .
```

Expect roughly 3–4 MB and no `.vcf`, `.bam` or `.fastq` files anywhere. `.gitignore` blocks them,
but look anyway — this is the step where a mistake becomes permanent and public.

---

## 3. Create the repository and push

Create an **empty** repository on github.com first — no README, no `.gitignore`, no license, or
the first push is rejected as non-fast-forward.

```bash
cd pgx-showcase
git init
git add -A
git commit -m "PGX Pipeline showcase: docs, results, figures, live demos, code excerpts"
git branch -M main
git remote add origin https://github.com/Zyrok12/pgx-showcase.git
git push -u origin main
```

> On Windows, do **not** wrap `git remote get-url origin` in a script with
> `$ErrorActionPreference = 'Stop'`. Git writes "No such remote" to stderr, which PowerShell 5.1
> promotes to a terminating error. Use `git remote` and test whether the output contains `origin`.

---

## 4. Enable GitHub Pages

In the repository: **Settings → Pages → Build and deployment**

- **Source:** Deploy from a branch
- **Branch:** `main`, folder **`/docs`**

Save. The site appears at `https://Zyrok12.github.io/pgx-showcase/` within a minute or two.

Verify these three load:

| URL | |
|---|---|
| `/pgx-showcase/` | landing page |
| `/pgx-showcase/demo/mega_vcf_explorer.html` | MegaVCF Explorer |
| `/pgx-showcase/demo/analytics_dashboard_interactive.html` | analytics dashboard |

Add the Pages URL to the repository's **About** panel so it shows at the top right — that's the
first thing a visitor clicks.

---

## 5. Suggested repository settings

- **Description:** *End-to-end pharmacogenomics pipeline for WES cohorts — ensemble variant
  calling, ML quality triage, 5-caller CNV consensus, database enrichment and interactive cohort
  analytics.*
- **Topics:** `bioinformatics` `pharmacogenomics` `genomics` `variant-calling` `cnv`
  `machine-learning` `wes` `python` `gatk` `deepvariant`
- Pin it on your GitHub profile.

---

## 6. Refreshing the demos after a new run

`scripts/refresh_demo.sh` re-copies the demo pages and figures from a completed pipeline run.

```bash
./scripts/refresh_demo.sh ~/PGX/runs/demo_pgx_showcase
```

It refuses to run against a run whose `run_manifest.json` is not marked `"synthetic_demo": true`,
so a real patient run can't be published by accident. Override only if you are certain:

```bash
./scripts/refresh_demo.sh --allow-non-synthetic /path/to/run
```

After refreshing, re-check for identifiers before committing:

```bash
grep -roE '\b[A-Z][0-9]{1,3}_S[0-9]+\b' docs/demo/ | sort -u
```

That pattern matches the non-synthetic sample-naming convention used in the private workspace
(a letter, digits, `_S`, digits). It must return nothing. Note that the pattern is described
rather than the IDs themselves — don't paste real sample IDs into this file to document them.

---

## One wording nit worth fixing

The interactive analytics dashboard carries this subtitle, baked into the generated HTML:

> *Synthetic data only. Designed for demos with higher-ups. It shows product capabilities, not
> clinical validation.*

"Demos with higher-ups" was written for an internal audience and reads oddly on a public
portfolio site. Change it in the generator (`scripts/generate_showcase_run.py` in the private
repo) and re-run the demo, then `refresh_demo.sh` — rather than hand-editing the copy here, which
the next refresh would overwrite. Something like *"Synthetic data only — demonstrates product
capabilities, not clinical validation"* keeps the disclaimer and drops the framing.

---

## What was already removed

The Home-screen screenshot (`figures/01-main-ui-home.png`) originally included a *Recent activity*
table listing non-synthetic sample IDs. It was cropped above that table before being added here.
If you ever regenerate that screenshot, crop it the same way — or take it on a machine whose run
history contains only demo runs.
