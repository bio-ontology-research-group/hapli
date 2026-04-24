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

## Protocol (pseudocode)

```
# One-time setup
aws s3 sync s3://human-pangenomics/working/HPRC_R2/assemblies/ assemblies/
download_grch38 + gencode_v45

# Per sample × haplotype — parallelisable via Snakemake
for sample in hprc_samples:
    for hap_idx in 1, 2:
        liftoff_out=lifted/{sample}.hap{hap_idx}.gff3
        liftoff -g anno/gencode.v45.gff3 -o $liftoff_out \
            --copies --polish assemblies/{sample}.hap{hap_idx}.fa ref/grch38.fa

        # Per-gene protein extraction + diff vs. reference protein
        python3 -m hapli.core.protein ... 

# Aggregate
python3 aggregate_per_gene.py \
    lifted/ protein_diffs/ \
    --out hprc_gene_lof_counts.tsv

# Compare to gnomAD
python3 compare_gnomad.py hprc_gene_lof_counts.tsv gnomad_v4_constraint.tsv \
    --out hprc_vs_gnomad.tsv
```

## Scripts to write (none yet exist)

- `fetch_hprc.sh`         — AWS S3 sync with resume
- `run_sample.smk`        — Snakemake rule for one sample
- `aggregate_per_gene.py` — population-level LoF count per gene
- `compare_gnomad.py`     — correlation HPRC LoF freq vs. gnomAD LoF freq

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
