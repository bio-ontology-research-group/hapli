# Phase 5-B — gnomAD compound-variant flip benchmark

**Status:** **runnable head-to-head on the curated paper-cases**, plus the
runbook for the full 1000G-scale Danecek 2017 reproduction.

## Quickstart — runnable head-to-head (now)

```bash
uv run python3 scripts/generate_paper_cases.py --all --out data/paper_cases --force
uv run python3 benchmarks/gnomad_flips/run_head_to_head.py --with-esm
```

For each curated paper case (compound-het LoF, frameshift+rescue, compound
missense pocket, SV deletion, synonymous control), runs `hapli analyze`
with the full pipeline (bcftools consensus → Liftoff → bcftools csq →
protein extraction → ESM2 epistasis → diploid aggregator) and projects
both the bcftools/csq view (per-haplotype consequence labels) and the
hapli view (per-haplotype score, `compound_het_lof` flag, ESM-residual
flag) to a comparable surface.

Result on the curated set:

| case                          | csq.hap1                | csq.hap2     | hapli h1/h2 | chl | epi | hapli-only signal               |
|---|---|---|---|---|---|---|
| 01_synonymous_pair            | synonymous              | synonymous   | 1.00 / 1.00 | .   | .   | —                                |
| 04_compound_het_lof           | stop_gained             | frameshift   | 0.00 / 0.49 | ✓   | .   | compound_het_lof_only_in_hapli   |
| 05_frameshift_rescue          | inframe_altering,missense | —          | 0.85 / 1.00 | .   | ✓   | epistasis_only_in_hapli          |
| 06_compound_missense_pocket   | missense                | —            | 0.96 / 1.00 | .   | ✓   | epistasis_only_in_hapli          |
| 08_symbolic_del_whole_gene    | —                       | —            | 0.00 / 1.00 | .   | .   | presence_hap1=deleted_only_in_hapli |

The four "hapli-only signal" rows are the gnomAD-flip-style disagreements
that exist by design — bcftools/csq has no diploid aggregation, no
epistasis residual, and no presence/SV call:

  - **04**: csq emits two LoF tags on different haps; hapli additionally
    emits `compound_het_lof=True` at the diploid level (the actionable
    autosomal-recessive signal).
  - **05**: csq's haplotype-aware machinery correctly tags the +1/-1 pair
    as `inframe_altering` (this is csq's strength); hapli's ESM-residual
    additionally flags the rescue at the protein-functional level
    (residual = +57 log-units for the 8M model).
  - **06**: csq emits two missense tags; hapli's epistasis residual flags
    the joint effect that variant-by-variant cannot reach.
  - **08**: csq emits nothing (the gene was removed by the SV);
    Liftoff reports `presence=deleted` and the diploid score=0.

Output: `benchmarks/gnomad_flips/results/head_to_head.json`.

## Claim under test

Haplotype-level gene-function assessment produces different calls from
variant-by-variant consequence tools on a curated set of compound-variant
loci. Specifically:

1. Reproduce BCFtools/csq's published 501-of-5019 flip figure from
   Danecek 2017 on 1000 Genomes Phase 3 phased VCF.
2. On that flip set, run hapli (Modes A and A+ESM) and produce a per-gene
   agreement matrix:
     `{hapli, Haplosaurus, BCFtools/csq, LOFTEE-per-variant}`
3. Enumerate hapli-specific flips (cases where hapli disagrees with ALL
   three competitors) and manually classify them against ClinVar / MAVE
   ground truth.

## Required inputs

| Artifact | Source | Size |
|---|---|---|
| 1000G Phase 3 phased VCF | `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/` | ~80 GB |
| GRCh37 reference FASTA | same FTP | 3 GB |
| GENCODE v19 / Ensembl GFF3 | Ensembl FTP | 50 MB |
| BCFtools/csq published flip list | supplementary to Danecek 2017 | < 1 MB |
| LOFTEE plugin + human_ancestor.fa.gz | gnomAD broadinstitute releases | 5 GB |
| VEP cache for GRCh37 | Ensembl VEP install | 15 GB |

## Protocol (pseudocode)

```
# Download + index inputs once
download_1000g_phase3 → 1000g/chr*.vcf.gz
download_grch37_ref   → ref/grch37.fa
download_gencode_v19  → anno/gencode.v19.gff3

# For each compound-variant site in the published flip list:
for site in flip_list:
    # Subset VCF to ±10 kb around the compound site
    bcftools view -r {chrom}:{start-10000}-{end+10000} 1000g/chr{chrom}.vcf.gz > case.vcf.gz

    # Run all four tools on the same (sample, gene)
    for sample in samples_carrying(compound_site):
        hapli analyze --gene {gene} --sample {sample} --vcf case.vcf.gz \
            --reference ref/grch37.fa --gff anno/gencode.v19.gff3 \
            --output-dir out/hapli/{sample}_{gene} --with-esm

        bcftools csq -p a -f ref/grch37.fa -g anno/gencode.v19.gff3 \
            -s {sample} case.vcf.gz > out/csq/{sample}_{gene}.vcf

        vep --plugin Haplosaurus -i case.vcf.gz -o out/haplo/{sample}_{gene}.txt

        vep --plugin LOFTEE -i case.vcf.gz -o out/loftee/{sample}_{gene}.txt

# Build per-(sample, gene) agreement matrix
python3 agreement.py out/ > flips.tsv

# Classify hapli-specific flips
python3 classify.py flips.tsv --against clinvar/ --against mavedb-cache/
```

## Scripts to write (none yet exist — placeholders)

- `fetch_inputs.sh`     — idempotent data downloader
- `subset_flips.py`     — extract per-site VCF slices
- `run_competitors.sh`  — run Haplosaurus + csq + LOFTEE
- `agreement.py`        — merge per-tool outputs into a matrix
- `classify.py`         — hapli-specific flips vs. ground truth

## Success criteria

- Reproduce Danecek 2017's 501 flips within ±10%.
- Identify at least one hapli-only flip that a human reviewer can defend
  as genuinely different from the per-variant / Haplosaurus consequence
  concatenation.
- Provide a short table (flip count × gene × disagreement type) for the paper.

## Out of scope for v1

- Rare-variant compound-het cases in populations under-represented in 1000G.
  The Saudi Human Genome Program cohort handles this as a discussion vignette.
- Structural variants in the flip set (1000G Phase 3 SVs are sparse). HPRC
  (Phase 5-C) is the better SV benchmark.
