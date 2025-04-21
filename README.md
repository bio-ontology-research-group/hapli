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

**Example YAML Configuration File (`config.yaml`):**

```yaml
gfa_file: data/example.gfa
gff3_file: data/example.gff3
reference_fasta: data/reference.fasta
output_file: output/annotations.tsv
log_level: INFO
```

