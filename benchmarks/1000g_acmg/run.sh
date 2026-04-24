#!/usr/bin/env bash
# 1000G + ACMG-SF-v3.2 benchmark orchestrator.
#
# Full pipeline:
#   1. fetch.py  → subsetted acmg.phased.vcf.gz (one-shot, ~hours)
#   2. snakemake → per-(sample, gene) hapli analyze + aggregate (scale-out)
#   3. analyze.py → paper table + PNGs comparing to gnomAD v4 constraint
#
# Usage:
#   ./benchmarks/1000g_acmg/run.sh smoke     # 5 samples × 3 genes
#   ./benchmarks/1000g_acmg/run.sh full      # 2504 × 79
#
# The smoke mode uses a tiny synthetic fixture from `smoke_fixture.py` and
# validates end-to-end in under a minute; use it before committing changes.

set -eo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
MODE="${1:-full}"

case "$MODE" in
    smoke)
        echo "==> smoke test (synthetic fixture, 5 samples × 3 genes)"
        uv run python3 "$REPO/benchmarks/1000g_acmg/smoke_test.py"
        ;;
    full)
        if [ -z "${SF_CONFIG:-}" ]; then
            SF_CONFIG="$REPO/benchmarks/1000g_acmg/config.yaml"
        fi
        if [ ! -f "$SF_CONFIG" ]; then
            echo "ERROR: Snakemake config not found at $SF_CONFIG" >&2
            echo "Copy config.template.yaml and fill in paths, or set SF_CONFIG." >&2
            exit 1
        fi
        echo "==> full run with config $SF_CONFIG"
        CORES="${CORES:-16}"
        cd "$REPO"
        uv run snakemake \
            -s "$REPO/workflows/Snakefile" \
            --configfile "$SF_CONFIG" \
            --cores "$CORES" \
            --rerun-incomplete \
            --keep-going
        echo "==> computing per-gene compound-het-LoF allele frequencies"
        OUT_DIR="$(uv run python3 -c "import yaml; print(yaml.safe_load(open('$SF_CONFIG'))['output_dir'])")"
        uv run python3 "$REPO/benchmarks/1000g_acmg/analyze.py" \
            --per-sample "$OUT_DIR/aggregate/per_sample.tsv" \
            --per-gene   "$OUT_DIR/aggregate/per_gene.tsv" \
            --gene-list  "$REPO/benchmarks/1000g_acmg/acmg_sf_v3_2.tsv" \
            --out        "$REPO/benchmarks/1000g_acmg/results/"
        ;;
    *)
        echo "Usage: $0 {smoke|full}" >&2
        exit 1
        ;;
esac
