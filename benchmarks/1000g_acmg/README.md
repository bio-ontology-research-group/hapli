# 1000 Genomes Phase 3 × ACMG-SF-v3.2 — population-scale Axis 2 benchmark

**Claim under test.** Across a public phased population cohort, hapli emits
genotype-level signals (Axis 2: `compound_het_lof`; Axis 1: SV presence calls)
that bcftools/csq emits only as per-variant records without aggregation. On
the ACMG Secondary Findings v3.2 gene panel (Miller 2023, *Genet Med*, 79
curated clinically actionable genes), we quantify how many of each signal
hapli produces across 2504 individuals, and whether the per-gene LoF allele
frequency is calibrated against gnomAD v4.

**Axis** mapping:

| hapli signal            | Aggregation axis      | csq equivalent             |
|---                      |---                    |---                         |
| `compound_het_lof=True` | Axis 2 (diploid)      | none (per-variant only)    |
| `presence=deleted`      | Axis 1 (SV presence)  | silent on symbolic records |
| `presence=duplicated`   | Axis 1                | silent                     |
| `presence=low_identity` | Axis 1                | silent                     |

**Cohort and data.**

- 1000 Genomes Phase 3 NYGC high-coverage phased WGS (Byrska-Bishop 2022,
  2504 unrelated + 698 related samples, GRCh38-aligned, Shapeit2-phased,
  SNV + INDEL + SV merged).
- URL: <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/>
- Access: public, no dbGaP required.
- Reference: GRCh38 primary assembly + GENCODE v45 GFF3.

**Gene panel.** ACMG SF v3.2 (79 genes) in
[`acmg_sf_v3_2.tsv`](acmg_sf_v3_2.tsv); columns include
`disease_category` and `inheritance` (AD / AR / XL) for downstream
per-gene disagreement reporting.

**Scripts.**

| file | purpose |
|---|---|
| [`acmg_sf_v3_2.tsv`](acmg_sf_v3_2.tsv) | ACMG SF v3.2 gene list (79 rows) |
| [`fetch.py`](fetch.py) | Region-subset NYGC phased VCFs to the ACMG gene intervals (one shot, ~hours) |
| [`config.template.yaml`](config.template.yaml) | Snakemake config template — copy, fill paths, run |
| [`run.sh`](run.sh) | Orchestrator: `smoke` or `full` mode |
| [`analyze.py`](analyze.py) | Post-aggregation per-gene summary + comparison to gnomAD |
| [`smoke_test.py`](smoke_test.py) | End-to-end validation on a 5-sample × 3-gene synthetic fixture |

## Running

### Smoke test (first)

```bash
./benchmarks/1000g_acmg/run.sh smoke
```

Builds a synthetic 5-sample × 3-gene fixture (GACMG1 = AR compound-het
carrier on S1, GACMG3 = AD SV-deleted on S3/hap1, GACMG2 = intact),
runs the full Snakemake pipeline + `hapli aggregate`, and asserts that
both the Axis-2 (compound-het-LoF) and Axis-1 (presence=deleted) signals
land in the aggregator output. Runs in <1 min; use before committing.

### Full run

```bash
# 1. Download references + annotation once:
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gff3.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip gencode.v45.annotation.gff3.gz hg38.fa.gz
samtools faidx hg38.fa

# 2. Subset NYGC VCFs to ACMG regions (~1–4 h, ~200 MB final VCF):
uv run python3 benchmarks/1000g_acmg/fetch.py \
    --genes benchmarks/1000g_acmg/acmg_sf_v3_2.tsv \
    --gff gencode.v45.annotation.gff3 \
    --out-dir data/1000g_acmg/ \
    --threads 4

# 3. Copy config template and fill in absolute paths:
cp benchmarks/1000g_acmg/config.template.yaml benchmarks/1000g_acmg/config.yaml
# edit reference, gff, vcf, (optional) gnomad_constraint paths

# 4. Run (recommend 16+ cores on a cluster):
CORES=32 ./benchmarks/1000g_acmg/run.sh full
```

Expected runtime: ~2 h / sample × 2504 samples × 79 genes. Snakemake
fans out over `(sample, gene)` wildcards; target a cluster with 32+
cores. `with_esm: false` in the default config — the bulk-DMS null
result in paper §3.3 means we don't want to burn GPU time on a
signal that doesn't show on this cohort either.

## Expected outputs

Under `results/1000g_acmg/`:

```
results/1000g_acmg/
├── <SAMPLE>_<GENE>_analysis.json        (2504 × 79 × 1)
├── <SAMPLE>_<GENE>_hap{1,2}.pep.fa
├── <SAMPLE>_<GENE>_haplotypes.fa
├── ...
└── aggregate/
    ├── per_sample.tsv                   (≈ 198K rows: sample × gene)
    └── per_gene.tsv                     (79 rows: one per ACMG gene)
```

`analyze.py` reads the two aggregator TSVs and writes:

```
benchmarks/1000g_acmg/results/
├── per_gene_summary.tsv                 (79 rows; paper Table 2)
├── tally.json                           (paper abstract figures)
└── lof_af_vs_oe_lof.png                 (paper Figure 3)
```

## Known caveats

- **Statistical phasing.** 1000G NYGC uses Shapeit2 + duoHMM pedigree
  phasing; it is not trio-exact. A small fraction of "compound-het-LoF"
  flags may be phasing errors (both LoF on the same haplotype) rather
  than true compound-het. Trio-confirmed cases are retained; the
  singleton-pedigree samples carry a caveat that is transparent in the
  §3.3-equivalent paper text.

- **SV call completeness.** The NYGC merged SV+INDEL VCF is not
  exhaustive for every SV class. Symbolic records hapli does not
  support (`<BND>`, `<CNV>`, `<INS>` without ALT sequence) are
  rejected at preflight with clear Mode-B recommendations; they appear
  as `not_run` in the aggregator for affected genes rather than
  silently failing.

- **Smoke-test artefacts.** The synthetic 3-gene fixture uses 300-bp
  toy genes — Liftoff's k-mer-based mapping is brittle on this scale
  (a single SNV can make the tiny gene unmappable). Real ACMG genes
  (10–100 kb, multi-exon) are well above Liftoff's noise floor; the
  smoke test's purpose is to exercise the pipeline-wiring integrity,
  not the per-gene call accuracy on 300-bp sequences.
