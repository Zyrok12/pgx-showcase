#!/usr/bin/env bash
#
# refresh_demo.sh - rebuild docs/demo/ from a completed PGX pipeline run.
#
#   ./scripts/refresh_demo.sh ~/PGX/runs/demo_pgx_showcase
#
# Refuses to run against anything not marked as a synthetic demo, so a real
# patient run cannot be published by accident. See docs/PUBLISHING.md.
#
set -euo pipefail

ALLOW_NON_SYNTHETIC=0
if [[ "${1:-}" == "--allow-non-synthetic" ]]; then
  ALLOW_NON_SYNTHETIC=1
  shift
fi

RUN="${1:-}"
if [[ -z "$RUN" ]]; then
  echo "usage: $0 [--allow-non-synthetic] <path-to-run-directory>" >&2
  exit 2
fi

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS="$RUN/results"
DEST="$REPO/docs/demo"

[[ -d "$RESULTS" ]] || { echo "error: no results/ under $RUN" >&2; exit 1; }

# ---------------------------------------------------------------- safety gate
MANIFEST="$RESULTS/run_manifest.json"
if [[ $ALLOW_NON_SYNTHETIC -eq 0 ]]; then
  if [[ ! -f "$MANIFEST" ]]; then
    echo "error: $MANIFEST not found - cannot confirm this run is synthetic." >&2
    echo "       Re-run with --allow-non-synthetic only if you are certain." >&2
    exit 1
  fi
  if ! grep -q '"synthetic_demo"[[:space:]]*:[[:space:]]*true' "$MANIFEST"; then
    echo "error: run is not marked \"synthetic_demo\": true in run_manifest.json." >&2
    echo "       This repository is public. Refusing to copy possible patient data." >&2
    exit 1
  fi
fi

echo "source : $RESULTS"
echo "target : $DEST"

rm -rf "$DEST"
mkdir -p "$DEST"/{analytics,advanced_analytics/ml_workbench,advanced_analytics/cnv,advanced_analytics/research_dl,tables}

# --- self-contained interactive pages -------------------------------------
cp "$RESULTS/mega_vcf/mega_vcf_explorer.html"                "$DEST/"
cp "$RESULTS/analytics/analytics_dashboard_interactive.html" "$DEST/"

# --- static figure dashboard -----------------------------------------------
# It lives in analytics/ and uses paths relative to that directory
# (bare *.png plus ../advanced_analytics/...), so it must stay there to resolve.
cp "$RESULTS/analytics/analytics_dashboard.html"        "$DEST/analytics/"
cp "$RESULTS"/analytics/*.png                           "$DEST/analytics/"
cp "$RESULTS"/advanced_analytics/*.png                  "$DEST/advanced_analytics/"           2>/dev/null || true
cp "$RESULTS"/advanced_analytics/ml_workbench/*.png     "$DEST/advanced_analytics/ml_workbench/"  2>/dev/null || true
cp "$RESULTS"/advanced_analytics/cnv/*.png              "$DEST/advanced_analytics/cnv/"           2>/dev/null || true
cp "$RESULTS"/advanced_analytics/research_dl/*.png      "$DEST/advanced_analytics/research_dl/"   2>/dev/null || true

# The dashboard references this one by bare filename even though the pipeline writes it to
# advanced_analytics/, so it needs a copy alongside the HTML or the image renders broken.
cp "$RESULTS/advanced_analytics/metabolizer_activity_heatmap.png" "$DEST/analytics/" 2>/dev/null || true

# --- small summary tables ---------------------------------------------------
for f in \
  "advanced_analytics/ml_workbench/ml_classifier_metrics.tsv" \
  "cnv/ensemble/consensus_cnv_segments.tsv" \
  "analytics/variant_significance_summary.tsv" ; do
  [[ -f "$RESULTS/$f" ]] && cp "$RESULTS/$f" "$DEST/tables/"
done
# The manifest records the absolute output path of the run. Scrub it - a public repo has no
# business disclosing the directory layout of the machine that produced it.
if [[ -f "$MANIFEST" ]]; then
  sed -E 's#(/home/[^"]*|/mnt/[a-z]/[^"]*)#<redacted>#g' "$MANIFEST" > "$DEST/tables/run_manifest.json"
fi

# ---------------------------------------------------------------- post-check
echo
echo "--- scanning for real sample identifiers ---"
# Matches the non-synthetic sample naming convention: letter, digits, "_S", digits.
if HITS=$(grep -rhoE '\b[A-Z][0-9]{1,3}_S[0-9]+\b' "$DEST" | sort -u) && [[ -n "$HITS" ]]; then
  echo "WARNING: possible real sample IDs found in the copied demo:" >&2
  echo "$HITS" >&2
  echo "Review before committing." >&2
  exit 1
fi
echo "clean - no real-format sample IDs found."

echo
echo "--- scanning for leaked local paths ---"
if PATHS=$(grep -rhoE '(/home/[a-z]+|/mnt/[a-z]/Users/[A-Za-z]+)' "$DEST" | sort -u) && [[ -n "$PATHS" ]]; then
  echo "WARNING: local filesystem paths found in the copied demo:" >&2
  echo "$PATHS" >&2
  exit 1
fi
echo "clean - no local paths found."

echo
echo "--- refreshed ---"
find "$DEST" -type f | sed "s|$DEST/||" | sort
echo
echo "$(find "$DEST" -type f | wc -l) files, $(du -sh "$DEST" | cut -f1)"
