#!/usr/bin/env bash
# HPRC R2 × Mode B benchmark orchestrator.
#
# Usage:
#   ./benchmarks/hprc/run.sh smoke     # synthetic 3-sample Mode B fixture
#   ./benchmarks/hprc/run.sh full      # whatever is in config.yaml's `haps:` block
#
# The smoke mode exercises the Mode B path (Liftoff on assembled FASTAs,
# no VCF) without downloading any real HPRC data; use before committing.

set -eo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
MODE="${1:-full}"

case "$MODE" in
    smoke)
        echo "==> smoke test (synthetic 3-sample Mode B fixture)"
        uv run python3 "$REPO/benchmarks/hprc/smoke_test.py"
        ;;
    full)
        if [ -z "${SF_CONFIG:-}" ]; then
            SF_CONFIG="$REPO/benchmarks/hprc/config.yaml"
        fi
        if [ ! -f "$SF_CONFIG" ]; then
            echo "ERROR: Snakemake config not found at $SF_CONFIG" >&2
            echo "Copy config.template.yaml, fetch assemblies (fetch.py), and fill in paths." >&2
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
        echo "==> computing per-gene presence + LoF distributions"
        OUT_DIR="$(uv run python3 -c "import yaml; print(yaml.safe_load(open('$SF_CONFIG'))['output_dir'])")"
        uv run python3 "$REPO/benchmarks/hprc/analyze.py" \
            --per-sample "$OUT_DIR/aggregate/per_sample.tsv" \
            --per-gene   "$OUT_DIR/aggregate/per_gene.tsv" \
            --sample-manifest "$REPO/benchmarks/hprc/sample_list.tsv" \
            --out        "$REPO/benchmarks/hprc/results/"
        ;;
    *)
        echo "Usage: $0 {smoke|full}" >&2
        exit 1
        ;;
esac
