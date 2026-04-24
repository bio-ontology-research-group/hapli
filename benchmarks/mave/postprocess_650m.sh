#!/usr/bin/env bash
# Post-process the 650M-param GB1 cluster run.
# Runs after `gb1_650m.sbatch` lands on unimatrix01:~/hapli-bench/results/.
#
# Usage: ./postprocess_650m.sh [JOB_ID]
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"

DATA="$REPO/data/mave"
RESULTS="$REPO/benchmarks/mave/results"
mkdir -p "$DATA" "$RESULTS"

echo "==> rsync 650M results from cluster"
rsync -e "ssh -o StrictHostKeyChecking=no" -av \
    "unimatrix01:~/hapli-bench/results/gb1.esm_650m_full.csv" \
    "$DATA/"

scored="$DATA/gb1.esm_650m_full.csv"
[ -s "$scored" ] || { echo "no result file: $scored" >&2; exit 1; }
echo "  fetched $(wc -l < "$scored") rows ($(stat -c %s "$scored") bytes)"

echo
echo "==> mave AUROC analysis (S_additive vs S_joint vs residual)"
uv run python3 benchmarks/mave/analyze_doubles.py \
    --scored "$scored" \
    --category missense_multi \
    --lof-threshold 0.1 \
    --lof-is low \
    --out "$RESULTS/gb1_650m_doubles.auroc.json"

echo
echo "==> experimental-epistasis correlation (residual vs E_exp)"
uv run python3 benchmarks/mave/experimental_epistasis.py \
    --scores "$DATA/urn_mavedb_00000105-a-1.scores.csv" \
    --scored "$scored" \
    --out "$RESULTS/gb1_650m_experimental_epistasis.json"

echo
echo "==> position-stratified epistasis"
uv run python3 benchmarks/mave/stratified_epistasis.py \
    --scored "$scored" \
    --bins 1,5,10,30,200 \
    --out "$RESULTS/gb1_650m_stratified.json"

echo
echo "==> paper figure (4-panel)"
uv run python3 benchmarks/mave/figure_paper.py \
    --scored "$scored" \
    --scores "$DATA/urn_mavedb_00000105-a-1.scores.csv" \
    --out "$RESULTS/gb1_650m_paper_figure.png" \
    --summary "$RESULTS/gb1_650m_paper_figure.json"

echo
echo "==> done. Headline AUROC + ρ:"
uv run python3 -c "
import json
a = json.load(open('$RESULTS/gb1_650m_doubles.auroc.json'))
e = json.load(open('$RESULTS/gb1_650m_experimental_epistasis.json'))
print(f'  rows kept (missense_multi):  {a[\"n_rows\"]}')
print(f'  S_additive AUROC:            {a[\"auroc\"][\"S_additive\"]:.4f}')
print(f'  S_joint    AUROC:            {a[\"auroc\"][\"S_joint\"]:.4f}')
print(f'  residual   AUROC:            {a[\"auroc\"][\"residual\"]:.4f}')
print(f'  ΔAUROC (joint − additive):   {a[\"auroc\"][\"S_joint\"] - a[\"auroc\"][\"S_additive\"]:+.4f}')
print()
ce = e['headline_correlations']
print(f'  Spearman(residual, E_exp):       {ce[\"residual_vs_E_exp\"][\"spearman\"]:+.4f}')
print(f'  Spearman(S_additive, E_exp):     {ce[\"s_additive_vs_E_exp\"][\"spearman\"]:+.4f}')
print(f'  ratio |residual| / |S_additive|: {abs(ce[\"residual_vs_E_exp\"][\"spearman\"]) / max(abs(ce[\"s_additive_vs_E_exp\"][\"spearman\"]), 0.001):.2f}x')
"
