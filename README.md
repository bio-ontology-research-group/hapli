# Haplotype Annotation Tool

## Purpose

This tool annotates features (like genes, exons) onto haplotypes represented within a Graphical Fragment Assembly (GFA) file. It uses a reference genome annotation (GFF3) and the GFA structure to project annotations onto specific paths within the graph.

## Installation

### Dependencies

The tool requires the following Python packages:
- PyGFA (for GFA parsing)
- Biopython (for GFF3 and FASTA parsing)
- NetworkX (for feature relationship graphs)
- PyYAML (for configuration files)
- argparse (for command-line arguments)

Install dependencies using:

```bash
pip install -r requirements.txt
```

## Usage

```bash
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta
```

Or using a configuration file:

```bash
python -m src.main --config-file src/config.yaml
```

## Path Selection

The tool allows you to select specific paths (haplotypes) for annotation:

### By Sample Name

```bash
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --sample-names sample1,sample2
```

### By Haplotype

```bash
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --haplotype-ids 1
```

### By Specific Path IDs

```bash
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --path-ids path1,path2
```

### Combined Selection

You can combine selection criteria:

```bash
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --sample-names sample1,sample2 --haplotype-ids 1
```

This will select haplotype 1 from both sample1 and sample2.

## Analysis Features

The tool includes advanced analysis capabilities to evaluate how features are represented across different haplotypes:

### Impact Classification

Features are classified according to their representation in different haplotypes:

- **Present**: Feature is present with high sequence identity (default: >90%)
- **Absent**: Feature could not be aligned to the path
- **Modified**: Feature has sequence changes but is largely intact
- **Truncated**: Feature is present but with significant length reduction
- **Expanded**: Feature is present but with significant length increase
- **Fragmented**: Feature is split across multiple locations

### Variant Detection

The tool identifies sequence changes within aligned features:

- **SNPs**: Single nucleotide polymorphisms
- **Insertions**: Additional sequence not present in the reference
- **Deletions**: Sequence missing compared to the reference
- **Complex Changes**: Substitutions involving multiple bases

### Feature Reconciliation

When child features (like exons) don't align properly within their parent features (like genes):

- Automatic adjustment of child features to fit within parent boundaries
- Identification of potentially incorrect feature boundaries
- Confidence scores for reconciliation decisions

### Analysis Output

Analysis results can be exported in multiple formats:

```bash
# Export as TSV (default)
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file results.tsv

# Export as JSON
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file results.json --output-format json
```

## Configuration

The tool uses a hierarchical configuration system:

1.  **Command-line arguments:** Highest priority. Any setting provided via a command-line flag overrides other sources.
2.  **Configuration file:** Settings can be specified in a YAML file using the `-c` or `--config-file` argument.
3.  **Defaults:** Built-in default values are used if a setting is not provided via CLI or config file.

**Required Parameters:**

These parameters *must* be provided either via the configuration file or command-line arguments:

*   `gfa_file`: Path to the input GFA file.
    *   CLI: `--gfa-file /path/to/your.gfa`
    *   YAML: `gfa_file: /path/to/your.gfa`
*   `gff3_file`: Path to the input GFF3 annotation file.
    *   CLI: `--gff3-file /path/to/your.gff3`
    *   YAML: `gff3_file: /path/to/your.gff3`
*   `reference_fasta`: Path to the reference FASTA file (used by GFF3).
    *   CLI: `--reference-fasta /path/to/your.fasta`
    *   YAML: `reference_fasta: /path/to/your.fasta`

**Optional Parameters:**

*   `output_file`: Path where the output annotations should be saved.
    *   CLI: `--output-file /path/to/output.tsv`
    *   YAML: `output_file: /path/to/output.tsv`
*   `log_level`: Sets the logging verbosity. Options: DEBUG, INFO, WARNING, ERROR, CRITICAL. (Default: INFO)
    *   CLI: `--log-level DEBUG`
    *   YAML: `log_level: DEBUG`
*   `sample_names`: Comma-separated list of sample names to include.
    *   CLI: `--sample-names sample1,sample2`
    *   YAML: `sample_names: [sample1, sample2]`
*   `haplotype_ids`: Comma-separated list of haplotype IDs to include.
    *   CLI: `--haplotype-ids 1,2`
    *   YAML: `haplotype_ids: [1, 2]`
*   `path_ids`: Comma-separated list of specific path IDs to include.
    *   CLI: `--path-ids path1,path2`
    *   YAML: `path_ids: [path1, path2]`

**Example YAML Configuration File (`config.yaml`):**

```yaml
gfa_file: data/example.gfa
gff3_file: data/example.gff3
reference_fasta: data/reference.fasta
output_file: output/annotations.tsv
log_level: INFO
sample_names: [sample1, sample2]
haplotype_ids: [1]
```

