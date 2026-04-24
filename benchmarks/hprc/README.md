# Phase 5-C — HPRC Release 2 scale demo

**Status:** runbook only. Each of the 200+ HPRC individuals produces two
3 Gb diploid-resolved assemblies; running the full pipeline per sample is
hours of wall time on a 16-core + 1 GPU workstation. Scaffolding the protocol
so the campaign can be fired off unattended once compute is provisioned.

## Claim under test

On real, non-simulated, fully phased diploid assemblies of healthy-ish
individuals, hapli's gene-level function calls are internally consistent
with published population LoF allele frequencies from gnomAD.

This is a *scale demo*, not a ground-truth comparison. The goal is to show
that:

1. The whole pipeline runs end-to-end on real human assemblies (not just
   the synthetic fixtures).
2. Gene-level LoF counts across ~200 individuals are in the right ballpark
   compared to gnomAD's per-gene LoF allele frequency times 2 × 200 = 400
   chromosomes, for the subset of genes where HPRC samples overlap with
   gnomAD's ancestry distribution.

## Required inputs

| Artifact | Source | Size |
|---|---|---|
| HPRC R2 phased assemblies (≥200 samples × 2 haps) | `s3://human-pangenomics/working/HPRC_R2/` | ~1.2 TB |
| GRCh38 reference FASTA | Ensembl FTP | 3 GB |
| GENCODE v45 GFF3 | `ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/` | 100 MB |
| gnomAD v4 LoF allele frequencies | `gnomad.broadinstitute.org/downloads` | 1 GB |

Mode B (no VCF — use assembly directly) applies here. `hapli assess` is
the Mode B CLI: it takes `--hap1 <fasta>` + `--hap2 <fasta>` + a reference
FASTA + GFF3 and runs Liftoff → protein extraction → protein diff → diploid
aggregation, with no `bcftools consensus` / `bcftools csq` pass (consequences
are inferred from the protein diff instead). End-to-end coverage:
`tests/test_paper_cases.py::test_mode_b_assess_produces_schema_v2_on_pre_assembled_haplotypes`.
For HPRC scale-out, the per-sample command is therefore:

```bash
uv run main.py assess --gene <GENE> --sample <HGxxx> \
    --hap1 HGxxx.hap1.fa --hap2 HGxxx.hap2.fa \
    --reference GRCh38.fa --gff gencode.v45.gff3 \
    --output-dir results/HGxxx/
```

and the aggregator `hapli aggregate` collapses the per-sample JSONs into
per-gene allele-frequency + compound-het-LoF TSVs.

## Scripts

| file | purpose |
|---|---|
| [`sample_list.tsv`](sample_list.tsv) | Curated 20-sample HPRC R2 manifest (hap1/hap2 URLs + population labels). Replace with the full-cohort manifest via `fetch.py --refresh-manifest` for the 200+ sample run |
| [`fetch.py`](fetch.py) | Downloads + unzips + faidx's hap1/hap2 FASTAs for selected samples; emits a YAML `haps:` block for Snakemake config |
| [`config.template.yaml`](config.template.yaml) | Mode B Snakemake config template — copy, fill `haps:` from `fetch.py --print-config`, run |
| [`run.sh`](run.sh) | Orchestrator: `smoke` or `full` |
| [`analyze.py`](analyze.py) | Post-aggregation: per-gene presence distribution + per-(gene, super-population) LoF allele frequency |
| [`smoke_test.py`](smoke_test.py) | End-to-end validation on a 3-sample synthetic Mode B fixture (HG_A / HG_B / HG_C across AFR / EAS / EUR) |

## Running

### Smoke test first

```bash
./benchmarks/hprc/run.sh smoke
```

Builds 3 synthetic "assembled" haplotype pairs (HG_A carries a whole-gene
deletion of a toy GACMG1 on hap1; HG_B reference-identical; HG_C a single
SNV), runs Mode B Snakemake → `hapli aggregate` → analyze.py, and asserts
that the Axis-1 SV-presence signal reaches Table A. Runs in <30 s; run
before committing.

### Full run (real HPRC data)

```bash
# 1. Refresh the manifest from the canonical HPRC index (optional):
uv run python3 benchmarks/hprc/fetch.py --refresh-manifest \
    --out-manifest benchmarks/hprc/sample_list.tsv

# 2. Fetch assemblies. Start with 20 samples (~24 GB) for a proof-of-concept;
#    scale to the full manifest for the paper:
uv run python3 benchmarks/hprc/fetch.py \
    --manifest benchmarks/hprc/sample_list.tsv \
    --out-dir data/hprc/assemblies/ \
    --print-config > benchmarks/hprc/haps_block.yaml

# 3. Download GRCh38 reference + GENCODE v45 GFF3 (see 1000g_acmg/README).

# 4. Copy config template, paste haps block, fill paths:
cp benchmarks/hprc/config.template.yaml benchmarks/hprc/config.yaml
# edit reference, gff, (optional) gnomad_constraint; replace `haps:` block
# with the output from step 2.

# 5. Run (Liftoff is the dominant cost; target 16+ cores):
CORES=32 ./benchmarks/hprc/run.sh full
```

Expected runtime: ~30 min per haplotype per gene-panel run (Liftoff
dominates). For 20 samples × 2 haplotypes × 79 genes ≈ 24 h on 32 cores
if genes are fanned out per invocation. Snakemake parallelises
`(sample, gene)` wildcards naturally.

## Expected outputs

Under `results/hprc/`:

```
results/hprc/
├── <SAMPLE>_<GENE>_analysis.json
├── <SAMPLE>_<GENE>_hap{1,2}.pep.fa
├── <SAMPLE>_<GENE>_hap{1,2}.lifted.gff3
├── ...
└── aggregate/
    ├── per_sample.tsv
    └── per_gene.tsv
```

`analyze.py` reads those and writes:

```
benchmarks/hprc/results/
├── per_gene_presence.tsv               (per-gene × status counts — paper Table 3)
├── per_gene_per_superpop_lof_af.tsv    (per-(gene, super-pop) LoF AF)
└── tally.json                          (paper abstract numbers)
```

## Success criteria

- Runs to completion on a 16-core + 1 GPU workstation in <24 hours for
  the full 200-sample cohort.
- Per-gene LoF frequency in HPRC correlates (Pearson > 0.4) with gnomAD
  v4 LoF allele frequency for a curated gene set (~2,000 constrained
  genes, pLI>0.9).
- Produces a per-sample JSON (one file per (sample, gene)) that clinical
  consumers can ingest for downstream interpretation.

## Notes

- HPRC is heavily enriched for African ancestry (HG002, HG003, HG004 are
  Ashkenazi trio; many others are YRI / Yoruban). gnomAD's allele
  frequencies are dominated by NFE (non-Finnish European). Expect some
  population-stratification noise in the correlation.
- Mode B support — i.e. `hapli assess --hap1 hap1.fa --hap2 hap2.fa
  --gff ref.gff3 --reference ref.fa` without a VCF — is a prerequisite.
  Currently (Phase 1–4) hapli's `analyze` subcommand requires a VCF;
  add `assess` before this benchmark can run for real.
