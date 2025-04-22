# GRCh38 Genome Processing Use Case

This document provides instructions for downloading the GRCh38 reference genome, generating variation graphs (GFA files) from the reference, and processing VCF files from the 1000 Genomes Project.

## Prerequisites

The following tools are required for these scripts:

- Python 3.7 or later
- BioPython (`pip install biopython`)
- tqdm for progress bars (`pip install tqdm`)
- pyfaidx for FASTA handling (`pip install pyfaidx`)
- vg toolkit for graph operations (https://github.com/vgteam/vg#installation)

For the VCF to GFA conversion:
- pysam for VCF processing (`pip install pysam`)
- gfapy for GFA manipulation (`pip install gfapy`)

Install the Python packages:

```bash
pip install biopython tqdm pyfaidx pysam gfapy
```

Install the vg toolkit following the instructions at: https://github.com/vgteam/vg#installation

## Downloading GRCh38 Reference Genome

The `download_grch38.py` script downloads the GRCh38 reference genome FASTA file from Ensembl or NCBI.

```bash
python scripts/download_grch38.py --output-dir data/reference --extract --source ensembl
```

Options:
- `--output-dir`: Directory to save the downloaded file (default: `data/reference`)
- `--extract`: Extract the downloaded gzip file
- `--source`: Source to download from (`ensembl` or `ncbi`, default: `ensembl`)

## Converting FASTA to GFA

The `fasta_to_gfa.py` script converts a FASTA file to GFA format.

```bash
python scripts/fasta_to_gfa.py --input data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output data/reference/GRCh38.gfa
```

Options:
- `--input, -i`: Input FASTA file path (required)
- `--output, -o`: Output GFA file path (default: input_path.gfa)
- `--use-external`: Use external tool (vg) for conversion if available

## Downloading VCF Files from 1000 Genomes Project

The `download_1000g_vcf.py` script downloads phased and unphased VCF files from the 1000 Genomes Project.

```bash
python scripts/download_1000g_vcf.py --output-dir data/vcf --phased 3 --unphased 3 --extract
```

Options:
- `--output-dir`: Directory to save the downloaded files (default: `data/vcf`)
- `--phased`: Number of phased VCF files to download (default: 3, max: 3)
- `--unphased`: Number of unphased VCF files to download (default: 3, max: 3)
- `--extract`: Extract the downloaded gzip files

## Converting VCF Files to GFA

The `vcf_to_gfa_converter.py` script converts VCF files to GFA format, either jointly or separately.

### Joint Conversion (all VCFs into a single GFA):

```bash
python scripts/vcf_to_gfa_converter.py \
  --vcf data/vcf/phased_1_ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf \
  --vcf data/vcf/phased_2_ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf \
  --reference data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output data/graphs/joint_graph.gfa \
  --mode joint
```

### Separate Conversion (one GFA per VCF):

```bash
python scripts/vcf_to_gfa_converter.py \
  --vcf data/vcf/phased_1_ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf \
  --vcf data/vcf/phased_2_ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf \
  --reference data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output data/graphs/ \
  --mode separate
```

Options:
- `--vcf, -v`: VCF file paths (can specify multiple times)
- `--reference, -r`: Reference genome FASTA file
- `--output, -o`: Output GFA file or directory
- `--mode`: Convert VCFs `joint` or `separate` (default: `joint`)
- `--region`: Restrict to region (e.g., "chr1:1000-2000")
- `--memory`: Memory limit in GB for vg (default: 16)

## Complete Workflow Example

Here's a complete workflow example:

```bash
# Create directories
mkdir -p data/reference data/vcf data/graphs

# Download GRCh38 reference (chromosome 22 only for faster processing)
python scripts/download_grch38.py --output-dir data/reference --extract

# Convert reference to GFA
python scripts/fasta_to_gfa.py --input data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output data/graphs/GRCh38_reference.gfa

# Download VCF files (3 phased, 3 unphased)
python scripts/download_1000g_vcf.py --output-dir data/vcf --phased 3 --unphased 3 --extract

# Convert VCFs to GFA separately
python scripts/vcf_to_gfa_converter.py \
  --vcf data/vcf/phased_1_ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf \
  --vcf data/vcf/unphased_1_HG00096.chrom22.ILLUMINA.bwa.GBR.low_coverage.20130415.bam.genotypes.vcf \
  --reference data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output data/graphs/ \
  --region "22" \
  --mode separate

# Convert VCFs to GFA jointly
python scripts/vcf_to_gfa_converter.py \
  --vcf data/vcf/phased_1_ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf \
  --vcf data/vcf/phased_2_ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf \
  --reference data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output data/graphs/joint_chr21_chr22.gfa \
  --mode joint
```

## Notes on Large Genome Processing

- The GRCh38 reference genome is large (~3GB), so downloading and processing may take time.
- Converting the entire genome to GFA can be memory-intensive.
- For testing, consider restricting to specific chromosomes or regions using the `--region` option.
- The vg toolkit requires significant memory for larger genomes. Adjust the `--memory` parameter as needed.
