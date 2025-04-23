# Haplotype Annotation Tool

## Purpose

This tool annotates features (like genes, exons) onto haplotypes represented within a Graphical Fragment Assembly (GFA) file. It uses a reference genome annotation (GFF3) and the GFA structure to project annotations onto specific paths within the graph, leveraging sequence alignment to accurately place features.

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
- **mappy** (Python wrapper for minimap2 sequence aligner)

Install dependencies using:

```bash
pip install -r requirements.txt
```
**Note:** `mappy` requires `minimap2` to be installed and accessible in your system's PATH, or it needs to compile it during installation, which might require a C compiler. Refer to the [mappy documentation](https://github.com/lh3/mappy) for details.

## Command-Line Usage

```bash
python -m src.main [options]
```

### Required Parameters

These parameters must be provided either via command-line or in a config file:

* `--gfa-file PATH`: Path to the input GFA file
* `--gff3-file PATH`: Path to the input GFF3 annotation file
* `--reference-fasta PATH`: Path to the reference FASTA file

### Configuration File

* `--config-file PATH`: Load settings from a YAML configuration file. CLI arguments override file settings.

### Output Options

* `--output-file PATH`: Path where output should be saved (default: stdout).
* `--output-format FORMAT`: Format for output reports (`tsv`, `json`, `rdf`). (default: `tsv`)
* `--rdf-format FORMAT`: RDF serialization if `output_format` is `rdf` (`turtle`, `xml`, `json-ld`, `ntriples`). (default: `turtle`)

### Path Selection

* `--sample-names LIST`: Comma-separated list of sample names to include.
* `--haplotype-ids LIST`: Comma-separated list of haplotype IDs to include.
* `--path-ids LIST`: Comma-separated list of specific path IDs to include.

### Alignment Options

*   `--minimap-preset STR`: Preset for minimap2 alignment (e.g., `splice`, `map-pb`, `map-ont`, `asm5`). Determines alignment parameters. (default: `splice`)

### Parallelization Options (Currently Affect Placeholder Components)

* `--num-workers N`: Number of parallel workers to use. (default: number of CPU cores)
* `--batch-size N`: Number of features/tasks to process in each batch. (default: 100)
* `--pool-type TYPE`: Worker pool type (`process` or `thread`). (default: `process`)

### Analysis & Reporting Options

*   `--comparative`: Generate a comparative report across multiple selected paths. (default: False)
*   `--impact-types LIST`: Filter report by comma-separated impact types (e.g., `PRESENT,ABSENT`).
*   `--region-start INT`: Start position for region filtering (based on path coordinates).
*   `--region-end INT`: End position for region filtering (based on path coordinates).

### Debug and Logging

* `--log-level LEVEL`: Logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). (default: `INFO`)
* `--save-intermediate`: Save intermediate data (config, selected paths, raw alignments TSV). (default: False)
* `--intermediate-dir PATH`: Directory to save intermediate data. (default: `intermediate_data`)

### Examples

```bash
# Basic usage with required files
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta

# Use a configuration file
python -m src.main --config-file config.yaml

# Specify a different minimap2 preset for PacBio HiFi reads and save JSON output
python -m src.main --gfa-file data/hifi.gfa --gff3-file data/anno.gff3 --reference-fasta data/ref.fasta \
    --minimap-preset map-hifi --output-format json --output-file output/hifi_annotations.json

# Save intermediate data, including raw alignments
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta \
    --save-intermediate --intermediate-dir intermediate_output
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

This will select haplotype 1 from both sample1 and sample2. If no selection criteria are provided, the tool attempts to process all paths found in the GFA file.

## Alignment Strategy

The tool employs a hierarchical, two-phase alignment strategy using **minimap2** (via the `mappy` Python wrapper) to accurately place reference features onto the selected GFA paths:

1.  **Path Sequence Extraction**: The full sequence for each selected GFA path is constructed by concatenating its segment sequences. Paths for which sequences cannot be constructed are skipped.
2.  **Phase 1: Parent Feature Alignment**:
    *   Top-level features (typically 'gene' features, or those without parents in the GFF hierarchy) are identified.
    *   The reference sequence of each parent feature is extracted from the reference FASTA.
    *   Each parent feature sequence is aligned against the full sequence of each selected path using `minimap2` with the configured preset (`--minimap-preset`).
    *   Successful alignments (based on minimap2's scoring and filtering) are recorded, including their coordinates on the path, CIGAR string, mapping quality, etc. A single feature may align multiple times to the same path.
3.  **Phase 2: Child Feature Alignment**:
    *   For each successfully aligned *instance* of a parent feature from Phase 1:
        *   The corresponding region of the path sequence is extracted.
        *   Child features (e.g., 'exon', 'mRNA') associated with that parent in the GFF hierarchy are identified.
        *   The reference sequence of each child feature is extracted.
        *   Each child feature sequence is aligned against the *extracted parent region sequence* using `minimap2` (potentially with different parameters optimized for shorter sequences, although currently uses the same preset).
        *   Successful child alignments are recorded, and their coordinates are mapped back to the full path's coordinate system. They are linked to their parent feature and the specific parent alignment instance they fall within via qualifiers in the resulting feature object.

This hierarchical approach ensures that child features are primarily searched for within the genomic context established by their parent feature's alignment, improving accuracy and efficiency.

Alignment results, including coordinates, CIGAR strings, scores, and parent-child relationships, are stored internally as modified `Bio.SeqFeature` objects and used for downstream analysis (Impact Classification, Variant Detection, Reconciliation). Raw alignment details can be saved to a TSV file using the `--save-intermediate` flag.

## Analysis Features (Placeholders Active)

The tool includes advanced analysis capabilities to evaluate how features are represented across different haplotypes. **Note: Currently, Impact Classification, Variant Detection, and Reconciliation use placeholder implementations.**

### Impact Classification (Placeholder)

Features are classified according to their representation in different haplotypes based on alignment results:

- **Present**: Feature is aligned with seemingly good metrics (placeholder logic uses mock identity/coverage).
- **Absent**: Feature could not be aligned to the path.
- **Modified**: Feature is aligned but metrics suggest differences (placeholder logic).
- **Truncated**: (Not implemented in placeholder)
- **Expanded**: (Not implemented in placeholder)
- **Fragmented**: (Not implemented in placeholder)

### Variant Detection (Placeholder)

The tool identifies sequence changes within aligned features:

- **SNPs**: (Placeholder may mock based on CIGAR)
- **Insertions**: (Placeholder may mock based on CIGAR)
- **Deletions**: (Placeholder may mock based on CIGAR)
- **Complex Changes**: (Placeholder mocks as 'COMPLEX' if CIGAR has I/D/X)

### Feature Reconciliation (Placeholder)

Handles cases where child features don't align properly within their parent features:

- Placeholder currently performs no reconciliation.

## Reporting (Using Placeholders)

The tool provides rich reporting capabilities to visualize and analyze the annotated features. **Note: Reporting currently uses placeholder analysis results.**

### Single Haplotype Reports

Generate detailed reports for a single haplotype (if multiple paths are analyzed without `--comparative`, only the first path's summary is reported):

```bash
# Generate a TSV report (default)
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file results.tsv

# Generate a JSON report
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file results.json --output-format json

# Generate an RDF report in Turtle format
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file results.ttl --output-format rdf --rdf-format turtle
```

### Comparative Reports

Generate reports comparing features across multiple haplotypes (requires `--comparative` flag):

```bash
# Generate a comparative TSV report for multiple haplotypes
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --path-ids path1,path2 --output-file comparative.tsv --comparative

# Generate a comparative JSON report
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --path-ids path1,path2 --output-file comparative.json --output-format json --comparative

# Generate a comparative RDF report in JSON-LD format
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --path-ids path1,path2 --output-file comparative.jsonld --output-format rdf --rdf-format json-ld --comparative
```

### Filtering Reports

Filter reports to focus on specific impact types or genomic regions (based on placeholder analysis):

```bash
# Filter by impact type (PRESENT, MODIFIED, ABSENT, etc.)
python -m src.main --gfa-file data/example.gfa --gff3-file data/example.gff3 --reference-fasta data/reference.fasta --output-file present_features.tsv --impact-types PRESENT

# Filter by genomic region (based on aligned coordinates on path)
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
# (Assuming a validation script or tool exists)
# Example using pyshex:
# pip install pyshex
# python -m pyshex.validate -s src/reporting/schemas/shex_schema.shex -d report.ttl -n <start_node> -st <start_shape>
```

### Advanced Comparative Analysis (Placeholders)

Identify consensus and discriminating features across multiple haplotypes:

```bash
# Placeholder commands - functionality not yet implemented
# python -m src.find_consensus --path-ids path1,path2,path3,path4 --threshold 0.9 --output-file consensus_features.json
# python -m src.find_discriminating --path-ids path1,path2,path3,path4 --output-file discriminating_features.json
```

## Configuration

The tool uses a hierarchical configuration system:

1.  **Command-line arguments:** Highest priority.
2.  **Configuration file:** Settings specified in a YAML file (`-c` or `--config-file`).
3.  **Defaults:** Built-in default values (see `src/config.py` or use `--help`).

**Required Parameters:**

*   `gfa_file`: Path to GFA file. (CLI: `--gfa-file`, YAML: `gfa_file`)
*   `gff3_file`: Path to GFF3 file. (CLI: `--gff3-file`, YAML: `gff3_file`)
*   `reference_fasta`: Path to reference FASTA. (CLI: `--reference-fasta`, YAML: `reference_fasta`)

**Optional Parameters (Selection):**

*   `output_file`: Output path. (CLI: `--output-file`, YAML: `output_file`)
*   `output_format`: `tsv`, `json`, `rdf`. (CLI: `--output-format`, YAML: `output_format`)
*   `rdf_format`: `turtle`, `xml`, `json-ld`, `ntriples`. (CLI: `--rdf-format`, YAML: `rdf_format`)
*   `log_level`: `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`. (CLI: `--log-level`, YAML: `log_level`)
*   `sample_names`: Comma-separated samples. (CLI: `--sample-names`, YAML: `sample_names: [s1, s2]`)
*   `haplotype_ids`: Comma-separated haplotypes. (CLI: `--haplotype-ids`, YAML: `haplotype_ids: [1, 2]`)
*   `path_ids`: Comma-separated paths. (CLI: `--path-ids`, YAML: `path_ids: [p1, p2]`)
*   `comparative`: Generate comparative report. (CLI: `--comparative`, YAML: `comparative: true`)
*   `impact_types`: Filter by impact types. (CLI: `--impact-types`, YAML: `impact_types: [PRESENT]`)
*   `region_start`: Region filter start. (CLI: `--region-start`, YAML: `region_start`)
*   `region_end`: Region filter end. (CLI: `--region-end`, YAML: `region_end`)
*   `minimap_preset`: Minimap2 preset. (CLI: `--minimap-preset`, YAML: `minimap_preset`)
*   `save_intermediate`: Save intermediate files. (CLI: `--save-intermediate`, YAML: `save_intermediate: true`)
*   `intermediate_dir`: Directory for intermediate files. (CLI: `--intermediate-dir`, YAML: `intermediate_dir`)
*   `num_workers`: Number of parallel workers. (CLI: `--num-workers`, YAML: `num_workers`)
*   `batch_size`: Batch size for parallel tasks. (CLI: `--batch-size`, YAML: `batch_size`)
*   `pool_type`: Parallel pool type (`process` or `thread`). (CLI: `--pool-type`, YAML: `pool_type`)


**Example YAML Configuration File (`config.yaml`):**

```yaml
# Input files (Required)
gfa_file: data/example.gfa
gff3_file: data/example.gff3
reference_fasta: data/reference.fasta

# Output settings
output_file: output/annotations.json
output_format: json
# rdf_format: turtle # Only needed if output_format is rdf

# Logging level
log_level: INFO

# Path selection (choose one or combine)
# path_ids: [sampleA_h1, sampleB_h1]
sample_names: [sampleA, sampleB]
haplotype_ids: [1]

# Alignment settings
minimap_preset: splice

# Reporting settings
comparative: true # Generate comparison between selected paths

# Debugging
save_intermediate: true
intermediate_dir: intermediate_output

# Analysis Filters (applied during reporting)
# impact_types: [PRESENT, MODIFIED]
# region_start: 5000
# region_end: 10000

# Parallelization (currently affects placeholders)
# num_workers: 8
# pool_type: process
```
```
DEVELOPER.md
