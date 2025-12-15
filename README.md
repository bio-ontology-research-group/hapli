# Hapli: Genotype-Centric Variant Analysis

Hapli is a toolkit for rethinking genomic variation. Instead of analyzing isolated variants, Hapli constructs complete **haplotypes** for an individual and tests whether genomic features (Genes, Transcripts, Exons) are functionally preserved within those haplotypes.

## Philosophy

1.  **Haplotypes over Variants**: Two variants might individually be benign but together destructive (or vice versa). Hapli generates the full sequence context.
2.  **Hierarchical Alignment**: Features are aligned hierarchically (Gene -> Transcript -> Exon). If a gene aligns well, we check its children. This improves accuracy and performance.
3.  **Functional Preservation**: The goal is to determine if the *functional unit* (e.g., the CDS) is intact in the new genotype.

## Installation

This project is managed with `uv`.

```bash
# Install dependencies
uv sync
```

Alternatively, standard pip installation:
```bash
pip install .
```

**External Dependencies:**
*   `minimap2` (must be in your PATH) for alignment.
*   `bcftools`/`tabix` (recommended for VCF handling).

## Usage

Hapli operates in two main steps:

### 1. Generate Haplotypes
Create full haplotype sequences (FASTA) from a reference genome and a phased VCF.

```bash
uv run main.py generate-haplotypes \
    --vcf path/to/phased.vcf.gz \
    --reference path/to/ref.fa \
    --region chr1:1000-2000 \
    --sample sample_001 \
    --output haplotypes.fa
```

### 2. Hierarchical Alignment
Align a specific gene and its sub-features to the generated haplotypes.

```bash
uv run main.py align-gene \
    --haplotypes haplotypes.fa \
    --gff features.gff3 \
    --reference ref.fa \
    --gene GeneName \
    --output results.json
```

## Automated Workflow

A `Snakemake` workflow is provided in `workflows/` to automate analysis across multiple samples and genes.

```bash
uv run snakemake -s workflows/Snakefile --configfile workflows/config.yaml --cores 4
```

## Output Format

The results are saved as JSON, preserving the feature hierarchy:

```json
{
  "sample_001_hap1": {
    "feature_id": "ENSG000...",
    "feature_type": "gene",
    "identity": 0.99,
    "mapq": 60,
    "children": [
      {
        "feature_id": "ENST000...",
        "feature_type": "transcript",
        "children": [
           {
             "feature_id": "exon1",
             "identity": 1.0
           }
        ]
      }
    ]
  }
}
```

## Development

Run tests:
```bash
uv run pytest
```