# Phase 5-A — MaveDB benchmark

The primary non-circular anchor for the paper. Tests how well hapli's ESM2-based
scoring recovers *experimental* functional scores from deep-mutational-scanning
assays, using single- and double-mutant libraries from MaveDB.

## Layout

```
fetch.py      # download a MaveDB scoreset + its UniProt reference protein
score.py      # compute per-variant ESM2 log-odds on the reference
analyze.py    # AUROC / AUPR against a MAVE-derived LoF label
run.sh        # one-liner orchestrator
results/      # auroc.json + auroc.png per scoreset
```

## Quickstart (BRCA1 RING, Starita 2015)

```bash
./run.sh urn:mavedb:00000081-a-1 2 high
```

Expected: AUROC ≈ 0.69 (at `lof-threshold=2, lof-is=high`) with the default
`esm2_t6_8M_UR50D` checkpoint on CPU. Larger checkpoints (650M, 3B) raise this
to the 0.75–0.80 range typically reported in the ESM-1v / ProteinGym literature.

## What `run.sh` does

1. **fetch.py** calls `api.mavedb.org/v1/score-sets/<urn>/scores` and writes the
   raw scores CSV, plus the metadata JSON and the UniProt reference FASTA
   (discovered from the scoreset's target-gene externalIdentifiers).
2. **score.py** does *one* ESM2 forward pass on the reference protein (cached
   to `~/.cache/hapli/esm` by sha256) and reads off, for each variant row,
   `log p(alt | ref_context) − log p(ref | ref_context)`. This is the
   single-variant baseline — equivalent to ESM-1v's wildtype-marginals score.
3. **analyze.py** builds a binary LoF label from the MAVE functional score
   (the polarity — whether LOW or HIGH MAVE means LoF — depends on the assay
   and is a CLI switch), computes AUROC / AUPR, writes a PR+ROC figure.

## Scoreset-specific conventions

MaveDB scores do not have a universal polarity. Check the scoreset's
`methodText` before picking `--lof-is`:

| Scoreset family | `--lof-is` | `--lof-threshold` | Rationale |
|---|---|---|---|
| BRCA1 RING/BRCT depletion (Starita 2015) | `high` | `2` | score = # replicates depleted; high = LoF |
| BRCA1 SGE (Findlay 2018) | `low`  | `-1` | function_score; negative = LoF |
| TP53 saturation (Giacomelli 2018) | `low` | `0` | z-score relative to wildtype |
| PTEN phosphatase (Mighell 2018) | `low` | `-0.5` | relative stability score |

## Phase 5-A2 — Double-mutant benchmark (GB1 pairwise, Olson 2014)

**Status: running.** `score_haplotypes.py` and `analyze_doubles.py` cover any
HGVS-p shape (single / double / multi missense, nonsense, frameshift, range
delins, SV-like N-terminal deletion) in one pass and produce a wide results
CSV (category, identity, hap_length, length_changed, premature_stop_at,
S_additive, S_joint, residual, flagged, hapli_score_heuristic).

Run on **GB1 IgG-FC binding (urn:mavedb:00000105-a-1, Olson 2014, pairwise
saturation scan at positions 10/38, 11/52, 28/31, ...)**, 2000 double mutants:

| Predictor              | AUROC   | Spearman |p-value|
|---|---|---|---|
| S_additive (baseline)  | 0.633   | −0.22    | 6e-24 |
| **S_joint (haplotype)**| **0.757** | **−0.45** | 4e-99 |
| residual               | 0.755   | −0.45    | 1e-101 |

- n = 2000 haplotypes, 953 LoF at `mave_score < 0.1`
- mean |residual| = 11.8 log-units — the haplotype-aware signal is large
- **+0.12 AUROC** improvement over variant-by-variant, **~2× Spearman**

The residual tracks S_joint almost perfectly — practically all the gain
from going joint over additive is in the residual, which is the claim the
paper is making: variant-by-variant baselines miss a substantial fraction
of the signal that's in the joint haplotype context.

Figure: `benchmarks/mave/results/gb1_doubles.auroc.png`

Run yourself:
```bash
uv run python3 benchmarks/mave/score_haplotypes.py \
    --scores data/mave/urn_mavedb_00000105-a-1.scores.csv \
    --reference data/mave/urn_mavedb_00000105-a-1.reference.fa \
    --out data/mave/gb1.esm_2000.csv \
    --max-rows 2000 --with-esm

uv run python3 benchmarks/mave/analyze_doubles.py \
    --scored data/mave/gb1.esm_2000.csv \
    --category missense_multi --lof-threshold 0.1 --lof-is low \
    --out benchmarks/mave/results/gb1_doubles.auroc.json
```

Compute budget: ~100 ms per unique haplotype with the 8M-param tiny ESM2 on
CPU. Full 315 k GB1 double-mutant catalogue is ~9 hours single-core; GPU
or 650M-param model gets you to meaningful AUROC improvements (typically
+0.05 to +0.15 over the 8M result) at ~20 min / 1000 rows.

## A note on PSD95 PDZ3 (urn:mavedb:00000053-a-1)

We attempted to add PSD95 PDZ3 (McLaughlin 2012, 648 022 variants, broader
position-pair distance distribution than GB1) as a second pairwise scan
to strengthen the position-stratified analysis. The MaveDB v1 `/scores`
endpoint returns HTTP 504 (Gateway Time-out) for this scoreset — the
single-response design appears not to scale to that record count. The
canonical workaround is to fetch the supplementary scores file from the
McLaughlin 2012 paper's data deposit (PDB / Zenodo) and convert it to
the MaveDB CSV layout. This is documented as a follow-on; the GB1
benchmark above is sufficient for the paper's headline epistasis claim.

## Phase 5-A3 — targets for the next sprint

Still to do for the paper's full epistasis story:

- **Experimental epistasis target** = `s(v1,v2) − [s(v1) + s(v2) − s(wt)]`.
  Requires both single- and double-mutant scores in the same experimental
  fitness scale. GB1 Olson 2014 has this (the single-mutant scan is a
  sibling scoreset not yet downloaded here). Target correlation: hapli's
  residual vs. experimental epistasis.
- Extension to the other classic pairwise scoresets: Araya 2012 hYAP65
  WW domain, PSD95 PDZ3 (McLaughlin 2012), Diss & Lehner 2018 Fos-Jun,
  Faure 2022/2024 mixed-effects saturation doubles.
- Stratified AUROC: positions close in sequence (< 10 aa) vs. far
  (> 30 aa) — epistasis signal should concentrate in the close-by pairs.

## Interpreting the result

- `roc_auc` > 0.7 on single mutants: hapli's ESM2 integration is working.
- `roc_auc` ≈ 0.5: flipped polarity — check `--lof-is`.
- `roc_auc` < 0.4: polarity almost certainly flipped; invert `--lof-is`.
- `roc_auc` between 0.55 and 0.65: expected for small ESM2 models on hard targets.

## Known limitations (this pass)

- Nonsense variants (`p.*`, `p.Ter`) are dropped because ESM2's amino-acid
  vocabulary has no stop-codon token. LOFTEE/AlphaMissense-flagged LoF calls
  cover those cases instead; the diploid aggregator's `premature_stop_at` is
  the hapli-native signal.
- Synonymous variants (`p.=`) are dropped — they're the ESM2 identity case
  and uninformative.
- HGVS protein notation with multi-site changes (`p.[A1V;C3T]`) is not yet
  parsed — this is where double-mutant scoresets will live, handled in the
  Phase 5-A2 script.
