# Hapli — Haplotype-aware Gene-Function Assessment

`hapli` reframes variant interpretation as a gene-level question — *is this
gene present and functional on each of an individual's two haplotypes?* —
instead of variant-by-variant HGVS consequence calling. Motivating cases
the tool handles correctly that variant-by-variant tools do not:

- **Compound-het LoF**: two LoF variants on different haplotypes →
  `evidence.diploid.compound_het_lof = True` (the actionable AR signal
  bcftools/csq can't aggregate).
- **Frameshift + rescue**: +1 / -1 indel pair on the same haplotype →
  protein is structurally restored after a small substitution window;
  `evidence.protein[*].frameshift_region.restored = True`.
- **Compound missense in a binding pocket**: two individually-benign
  missense variants on the same haplotype → ESM2 epistasis residual
  `evidence.epistasis[*].flagged = True`.
- **SV gene deletion**: symbolic `<DEL>` removing a whole gene →
  Liftoff reports `presence = deleted`, diploid score = 0.

## Install

```bash
uv sync
# optional ML extra (torch + fair-esm) for the epistasis residual:
uv sync --extra ml
# external binaries (must be on PATH): bcftools (≥1.15), samtools, minimap2
# Liftoff installs in its own env to avoid pysam pin conflict:
uv tool install git+https://github.com/agshumate/Liftoff
```

## Two operating modes

### Mode A — phased VCF input

```bash
uv run hapli analyze \
    --gene BRCA1 --sample HG002 \
    --vcf phased.vcf.gz \
    --reference grch38.fa \
    --gff gencode.v45.gff3 \
    --output-dir results/ \
    [--with-esm]                                # ESM2 epistasis residual
    [--esm-model esm2_t33_650M_UR50D]           # which ESM2 checkpoint (default 8M)
    [--alphamissense-table am_scores.tsv.gz]    # per-variant pathogenicity baseline
    [--clinvar-vcf clinvar.vcf.gz]              # CLNSIG annotation per consequence
    [--gnomad-constraint constraint.tsv]        # pLI / MisZ priors
    [--clingen-dosage clingen.tsv]              # haploinsufficiency priors
```

Pipeline: `bcftools consensus` → Liftoff → `bcftools csq -p a` → protein
extraction + diff → optional ESM2 epistasis → diploid report.

### Mode B — pre-assembled haplotypes (HPRC, hifiasm, Verkko)

```bash
uv run hapli assess \
    --gene BRCA1 --sample HG002 \
    --hap1 HG002.hap1.fa --hap2 HG002.hap2.fa \
    --reference grch38.fa --gff gencode.v45.gff3 \
    --output-dir results/ \
    [--with-esm] [--gnomad-constraint ...] [--clingen-dosage ...]
```

No VCF; variant-equivalent information recovered from the protein diff.

### Other subcommands

```bash
uv run hapli explore  results/HG002_BRCA1_analysis.json   # Textual TUI
uv run hapli interpret results/HG002_BRCA1_analysis.json   # LLM narrative (OPENROUTER_API_KEY)
uv run hapli version
```

## Output schema (v2)

The output JSON keeps the v1 top-level keys (`gene`, `sample`, `region`,
`transcripts`) for back-compat with the TUI and LLM consumer, and adds a
new `evidence` block:

```
evidence.presence       — Liftoff per-haplotype call: intact/partial/low_identity/
                          duplicated/deleted/uncertain
evidence.consequence    — bcftools csq per-haplotype per-transcript HGVS records
evidence.protein        — reference-vs-haplotype protein diff: identity, premature_stop,
                          frameshift_region, substitutions
evidence.epistasis      — additive vs. joint ESM2 log-likelihood residual (--with-esm)
evidence.diploid        — two-number summary (hap1_score, hap2_score) + compound_het_lof
                          flag + constraint priors. Deliberately NOT a single diploid call.
evidence.missense_agg   — AlphaMissense per-transcript aggregation (--alphamissense-table)
evidence.splice / lof   — placeholders for SpliceAI / LOFTEE integration (deferred to v2)
evidence.consequence[*].clnsig — ClinVar clnsig annotation (--clinvar-vcf)
```

## Test data generation

```bash
# 12 named motivating cases used by the paper + benchmarks:
uv run python3 scripts/generate_paper_cases.py --all --out data/paper_cases
# All 12 wired: 01 synonymous pair, 02 single benign missense, 03 single stop-gain,
#   04 compound-het LoF, 05 frameshift rescue, 06 compound missense pocket,
#   07 splice donor + cryptic, 08 symbolic <DEL> whole-gene, 09 upstream DEL shift,
#   10 <INV> gene body, 11 <DUP:TANDEM>, 12 in-frame 3 bp del.
```

## Benchmarks

```bash
# 5-A — MAVE single-mutant baseline (BRCA1 RING, Starita 2015):
benchmarks/mave/run.sh urn:mavedb:00000081-a-1 2 high
#   → AUROC 0.69 with the tiny 8M ESM2 (CPU)

# 5-A2 — MAVE double-mutant epistasis (GB1 pairwise, Olson 2014):
uv run python3 benchmarks/mave/score_haplotypes.py \
    --scores GB1.scores.csv --reference GB1.reference.fa \
    --out gb1.scored.csv --max-rows 2000 --with-esm
uv run python3 benchmarks/mave/analyze_doubles.py \
    --scored gb1.scored.csv --lof-threshold 0.1 --lof-is low \
    --out results.json
#   → AUROC: S_joint 0.76 vs S_additive 0.63 (8M model on 2K rows)
#     Spearman vs experimental epistasis E_exp: residual |0.42| vs additive |0.21|

# 5-B — head-to-head with bcftools/csq on curated paper cases:
uv run python3 benchmarks/gnomad_flips/run_head_to_head.py --with-esm
#   → table of disagreement axes (compound_het_lof, epistasis_flagged,
#     presence=deleted) hapli emits but csq cannot

# 5-D — synthetic SV-haplotype shape coverage:
uv run python3 benchmarks/sv/run.py --out benchmarks/sv/results
#   → 12 cases (single missense, double missense, nonsense, frameshift,
#     range del, N-terminal del, compound LoF+missense), 45/45 assertions
```

## Repository layout

```
hapli/
  cli/                 — Typer CLI: analyze (Mode A) + assess (Mode B) + explore + interpret
  workflow/pipeline.py — orchestrator
  external/            — wrappers for bcftools (consensus, csq), liftoff, alphamissense, constraint TSVs
  core/                — schema (pydantic v2 + fallback shim), HGVS parser, protein extraction + diff
  interpretation/      — ESM2 scoring, epistasis residual, diploid aggregator, LLM
benchmarks/
  mave/                — MaveDB integration: fetch, score, analyze (singles + doubles + experimental epistasis)
  gnomad_flips/        — head-to-head with bcftools/csq on curated paper cases
  sv/                  — synthetic SV-haplotype shape coverage
  hprc/                — Phase 5-C runbook (200-sample scale demo on HPRC R2)
scripts/
  generate_paper_cases.py — synthetic case fixtures used by tests + benchmarks
data/
  test/                — small synthetic fixture used by the showboat walkthrough
  paper_cases/         — generated motivating cases (gitignored)
  mave/                — downloaded MaveDB scoresets (gitignored)
docs/
  walkthrough.md       — showboat-driven executable end-to-end walkthrough (regression-tested)
tests/                 — pytest suite, currently 163 tests
```

## Status

| Phase | Description | Status |
|---|---|---|
| 0 | SV-corruption fix + tests + showboat walkthrough | ✅ |
| 1 | External tool integration (bcftools, Liftoff, csq, AlphaMissense lookup) | ✅ |
| 2 | Haplotype-resolved protein extraction + diff | ✅ |
| 3 | ESM2 additive-vs-joint epistasis residual | ✅ |
| 4 | Diploid two-number report + gnomAD/ClinGen priors | ✅ |
| 5-A | MAVE single-mutant baseline | ✅ (AUROC 0.69) |
| 5-A2 | MAVE double-mutant epistasis | ✅ (AUROC 0.76, ρ vs E_exp = 0.42) |
| 5-A2-pub | Full GB1 (315k) with ESM2 650M on GPU | running on unimatrix01 node005 (resumable; ~120K rows done on first cluster submission, resume-in-progress) |
| 5-B | hapli vs bcftools/csq head-to-head | ✅ on curated cases (1000G runbook documented) |
| 5-C | HPRC 200-sample scale demo | runbook ready, needs cluster compute |
| 5-D | SV-haplotype shape coverage | ✅ (45/45 assertions) |
