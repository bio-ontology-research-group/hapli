# Paper-cases benchmarks

Artifact-producing scripts that run the full hapli pipeline on every
`data/paper_cases/` bundle and emit per-case summary TSVs. These are the
source of paper Table 1 and paper Table 2 / Figure 2 (controlled-residual
validation).

## Scripts

### `run_summary.py` → `summary.tsv` (paper Table 1)

Runs `hapli analyze` on each of the 12 cases and emits one row per case:
per-haplotype Liftoff presence, diploid score, compound-het-LoF flag,
hap1/hap2 identity and premature-stop, frameshift-restored flag,
per-variant consequence categories, and (if `--with-esm` was used) the
joint-residual + flagged column.

```bash
uv run python3 benchmarks/paper_cases/run_summary.py \
    --cases-dir data/paper_cases \
    --out benchmarks/paper_cases/summary.tsv
```

Key rows in the current `summary.tsv`:

| case                         | hap1_presence | hap1_score | compound_het_lof | frameshift_restored |
|---                           |---            |---         |---               |---                  |
| 04 compound_het_lof          | uncertain     | 0.00       | **True**         | —                   |
| 05 frameshift_rescue         | intact        | 0.85       | —                | **True**            |
| 08 symbolic_del_whole_gene   | **deleted**   | 0.00       | —                | —                   |
| 10 inversion_body            | low_identity  | 0.89       | —                | —                   |
| 11 tandem_dup_small_gene     | duplicated    | 1.00       | —                | —                   |

### `run_epistasis.py` → `epistasis.tsv` (paper Figure 2)

Runs `hapli analyze --with-esm` on the subset of cases that have ≥ 2 coding
substitutions on one haplotype; records the S_additive / S_joint / residual
triple.

```bash
uv run python3 benchmarks/paper_cases/run_epistasis.py \
    --cases-dir data/paper_cases \
    --out benchmarks/paper_cases/epistasis.tsv \
    --esm-model esm2_t6_8M_UR50D
```

Current numbers (8M model, CPU, <10 s per case):

| case                         | n_vars | S_additive | S_joint | residual   | flagged |
|---                           |---     |---         |---      |---         |---      |
| 05 frameshift_rescue (hero)  | 6      | -58.48     | -1.43   | **+57.05** | ✓       |
| 06 compound_missense_pocket  | 2      | -21.91     | -9.70   | **+12.20** | ✓       |
| 07 splice_donor_cryptic      | 0      | 0.00       | 0.00    | 0.00       | —       |
| 10 inversion_body            | 0      | 0.00       | 0.00    | 0.00       | —       |

Cases 07 (splice-only) and 10 (SV) correctly emit n=0 — the residual is
defined for ≥ 2 coding substitutions on the same transcript. The residual
interpretation: positive = joint protein is *less* deleterious than the
sum of per-variant log-odds (potential rescue); negative = joint protein
is *more* deleterious (potential synergistic disruption).
