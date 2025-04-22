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
- rdflib (for RDF report generation)
- pyshex (for RDF validation)

Install dependencies using:

```bash
pip install -r requirements.txt
```

## Command-Line Usage

```bash
python -m src.main [options]
```

### Required Parameters

These parameters must be provided either via command-line or in a config file:

* `--gfa-file PATH`: Path to the input GFA file
* `--gff3-file PATH`: Path to the input GFF3 annotation file
* `--reference-fasta PATH`: Path to the reference FASTA file
* `--config-file PATH`: Load settings from a YAML configuration file

### Output Options

* `--output-file PATH`: Path where output should be saved (default: stdout)
* `--output-format FORMAT`: Format for output reports (tsv, json, rdf) (default: tsv)
* `--rdf-format FORMAT`: RDF serialization if output_format is rdf (turtle, xml, json-ld, ntriples)

### Path Selection

* `--sample-names LIST`: Comma-separated list of sample names to include
* `--haplotype-ids LIST`: Comma-separated list of haplotype IDs to include
* `--path-ids LIST`: Comma-separated list of specific path IDs to include

### Parallelization Options

* `--num-workers N`: Number of parallel workers to use (default: number of CPU cores)
* `--batch-size N`: Number of features to process in each batch (default: 100)
* `--pool-type TYPE`: Worker pool type ('process' or 'thread') (default: process)

### Debug and Logging

* `--log-level LEVEL`: Logging verbosity (DEBUG, INFO, WARNING, ERROR, CRITICAL) (default: INFO)
* `--save-intermediate`: Save intermediate data for later visualization
* `--intermediate-dir PATH`: Directory to save intermediate data (default: intermediate_data)

### Examples

```bash
# Basic usage with required files
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta

# Use a configuration file
python -m src.main --config-file config.yaml

# Enable parallel processing
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --num-workers 4

# Save intermediate data for later visualization
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --save-intermediate --intermediate-dir intermediate_data
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

## Reporting

The tool provides rich reporting capabilities to visualize and analyze the annotated features:

### Single Haplotype Reports

Generate detailed reports for a single haplotype:

```bash
# Generate a TSV report (default)
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file results.tsv

# Generate a JSON report
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file results.json --output-format json

# Generate an RDF report in Turtle format
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file results.ttl --output-format rdf --rdf-format turtle
```

### Comparative Reports

Generate reports comparing features across multiple haplotypes:

```bash
# Generate a comparative TSV report for multiple haplotypes
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --path-ids path1,path2 --output-file comparative.tsv --comparative

# Generate a comparative JSON report
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --path-ids path1,path2 --output-file comparative.json --output-format json --comparative

# Generate a comparative RDF report in JSON-LD format
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --path-ids path1,path2 --output-file comparative.jsonld --output-format rdf --rdf-format json-ld --comparative
```

### Filtering Reports

Filter reports to focus on specific impact types or genomic regions:

```bash
# Filter by impact type (PRESENT, MODIFIED, ABSENT, etc.)
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file present_features.tsv --impact-types PRESENT

# Filter by genomic region
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file region_features.tsv --region-start 5000 --region-end 10000

# Combine filters
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file filtered.tsv --impact-types PRESENT,MODIFIED --region-start 5000 --region-end 10000
```

### RDF Reports

The tool supports multiple RDF serialization formats for semantic web applications:

```bash
# Generate RDF in Turtle format
python -m src.main --output-file report.ttl --output-format rdf --rdf-format turtle

# Generate RDF in XML format
python -m src.main --output-file report.rdf --output-format rdf --rdf-format xml

# Generate RDF in JSON-LD format
python -m src.main --output-file report.jsonld --output-format rdf --rdf-format json-ld

# Generate RDF in N-Triples format
python -m src.main --output-file report.nt --output-format rdf --rdf-format ntriples
```

### RDF Schema Validation

Validate RDF reports against the provided ShEx schema:

```bash
# Validate a report against the ShEx schema
python -m src.validate_rdf --rdf-file report.ttl --shex-file src/reporting/schemas/shex_schema.shex
```

### Advanced Comparative Analysis

Identify consensus and discriminating features across multiple haplotypes:

```bash
# Identify features present in at least 90% of haplotypes
python -m src.find_consensus --path-ids path1,path2,path3,path4 --threshold 0.9 --output-file consensus_features.json

# Identify features that discriminate between haplotypes
python -m src.find_discriminating --path-ids path1,path2,path3,path4 --output-file discriminating_features.json
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
*   `output_format`: Format for output reports (tsv, json, rdf).
    *   CLI: `--output-format json`
    *   YAML: `output_format: json`
*   `rdf_format`: RDF serialization format if output_format is rdf (turtle, xml, json-ld, ntriples).
    *   CLI: `--rdf-format turtle`
    *   YAML: `rdf_format: turtle`
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
*   `comparative`: Generate a comparative report (for multiple paths).
    *   CLI: `--comparative`
    *   YAML: `comparative: true`
*   `impact_types`: Filter by impact types.
    *   CLI: `--impact-types PRESENT,MODIFIED`
    *   YAML: `impact_types: [PRESENT, MODIFIED]`
*   `region_start`: Start position for region filtering.
    *   CLI: `--region-start 5000`
    *   YAML: `region_start: 5000`
*   `region_end`: End position for region filtering.
    *   CLI: `--region-end 10000`
    *   YAML: `region_end: 10000`

**Example YAML Configuration File (`config.yaml`):**

```yaml
gfa_file: data/example.gfa
gff3_file: data/example.gff3
reference_fasta: data/reference.fasta
output_file: output/annotations.tsv
output_format: json
rdf_format: turtle
log_level: INFO
sample_names: [sample1, sample2]
haplotype_ids: [1]
comparative: true
impact_types: [PRESENT, MODIFIED]
region_start: 5000
region_end: 10000
```

