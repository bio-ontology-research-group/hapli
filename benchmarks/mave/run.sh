#!/usr/bin/env bash
# End-to-end MAVE benchmark: fetch → score → analyse on a MaveDB scoreset.
#
# Usage:
#   ./run.sh <urn> [<lof-threshold>] [<lof-is=low|high>]
#
# Example:
#   ./run.sh urn:mavedb:00000081-a-1 2 high          # BRCA1 RING (Starita 2015)
#
# Outputs:
#   data/mave/<urn>.scores.csv
#   data/mave/<urn>.reference.fa
#   data/mave/<urn>.scored.csv
#   benchmarks/mave/results/<urn>.auroc.json
#   benchmarks/mave/results/<urn>.auroc.png
#
# The heavy step is `score.py` — ~30 s on CPU for a 1.8k-aa reference with the
# 8M model. Repeat runs of the same URN hit the ESM per-position-log-probs
# cache (~$HAPLI_CACHE_DIR) and finish in seconds.
set -euo pipefail

URN="${1:?first arg: MaveDB scoreset URN (e.g. urn:mavedb:00000081-a-1)}"
THRESHOLD="${2:-2}"
LOF_IS="${3:-high}"

# Always run relative to the project root so relative paths (benchmarks/, data/)
# resolve sanely regardless of caller cwd.
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

SAFE=$(echo "$URN" | tr ':/' '__')
DATA=data/mave
RESULTS=benchmarks/mave/results
mkdir -p "$DATA" "$RESULTS"

echo "==> fetch"
uv run python3 benchmarks/mave/fetch.py "$URN" -o "$DATA"

echo "==> score"
uv run python3 benchmarks/mave/score.py \
    --scores "$DATA/${SAFE}.scores.csv" \
    --reference "$DATA/${SAFE}.reference.fa" \
    --out "$DATA/${SAFE}.scored.csv"

echo "==> analyse (LoF threshold=$THRESHOLD, lof_is=$LOF_IS)"
uv run python3 benchmarks/mave/analyze.py \
    --scored "$DATA/${SAFE}.scored.csv" \
    --lof-threshold "$THRESHOLD" \
    --lof-is "$LOF_IS" \
    --out-dir "$RESULTS"

echo
echo "Summary: $RESULTS/${SAFE}.scored.auroc.json"
cat "$RESULTS/${SAFE}.scored.auroc.json"
