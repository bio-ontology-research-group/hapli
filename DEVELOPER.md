# Developer Guide

## Implementation Details

### Configuration Management (`src/config.py`)

The application uses a dedicated `Config` class (`src/config.py`) to manage configuration settings.

**Key Features:**

*   **Hierarchical Loading:** Configuration is loaded in the following order of precedence (highest first):
    1.  Command-line arguments (parsed using `argparse`).
    2.  YAML configuration file (specified via `-c` or `--config-file`, loaded using `PyYAML`).
    3.  Default values defined within the `Config` class.
*   **Libraries Used:**
    *   `argparse`: For parsing command-line arguments.
    *   `PyYAML` (dependency): For parsing YAML configuration files.
    *   `os`: For path manipulations and existence checks.
*   **Validation:** After loading, the `Config` class automatically validates:
    *   The presence of all required parameters (`gfa_file`, `gff3_file`, `reference_fasta`).
    *   The existence of the files specified for these required parameters.
    *   An appropriate `ConfigurationError` is raised if validation fails.
*   **Access:** Configuration values can be accessed using the `get(key, default)` method of a `Config` instance. Resource file paths are specifically tracked and can be retrieved using `get_resource_files()`.
*   **Extensibility:** New configuration parameters can be added by:
    1.  Adding a default value in `DEFAULT_CONFIG` (if applicable).
    2.  Adding a corresponding `add_argument` call in `_setup_argparser`.
    3.  Adding the parameter key to `REQUIRED_PARAMS` if it's mandatory.

### Path Analysis Module (`src/path_analysis.py`)

The Path Analysis module provides functionality to:

1. Extract paths from a parsed GFA file
2. Identify potential haplotype relationships between paths
3. Group paths by sample based on naming conventions or metadata
4. Select specific paths for annotation

**Key Components:**

- **PathAnalyzer Class:** The main class that handles all path analysis.
  - **load_gfa()**: Loads paths from a GFA object
  - **group_paths_by_sample()**: Groups paths by inferred sample name
  - **identify_haplotypes()**: Identifies potential haplotype relationships
  - **select_paths()**: Selects paths based on sample, haplotype, or direct path IDs
  - **get_path_segments()**: Retrieves segments for a specific path

**Haplotype Naming Patterns:**

The module recognizes several common naming patterns for haplotypes in GFA paths:
- `sample_hap1`, `sample_hap2`
- `sample_h1`, `sample_h2`
- `sample.1`, `sample.2`, `sample_1`, `sample_2`
- `hap1_sample`, `hap2_sample`
- `h1_sample`, `h2_sample`

**Sample and Haplotype Inference:**

- Sample names are extracted from path IDs based on the above patterns
- Paths without recognizable patterns are considered as their own sample group
- The module also looks for metadata tags (SM for sample, HP for haplotype) if available

### Analysis Module (`src/analysis/`)

The Analysis module provides comprehensive functionality to evaluate and characterize feature alignments:

#### Impact Classifier (`src/analysis/impact_classifier.py`)

* Categorizes aligned features into impact types (present, absent, modified, etc.)
* Uses sequence identity and coverage thresholds to determine status
* Considers feature length changes to identify truncations and expansions

#### Variant Detector (`src/analysis/variant_detector.py`)

* Identifies sequence changes between reference and aligned features
* Detects SNPs, insertions, and deletions
* Uses CIGAR strings when available for precise variant calling
* Provides filtering and summarization methods

#### Feature Reconciler (`src/analysis/reconciliation.py`)

* Handles cases where child features don't align within parent boundaries
* Implements strategies like adjusting child features, suggesting parent adjustments, or orphaning children
* Maintains confidence scores for each reconciliation decision
* Preserves feature hierarchies when possible

#### Summary Generator (`src/analysis/summary_generator.py`)

* Aggregates results from impact classification, variant detection, and reconciliation
* Provides comprehensive feature-level and path-level summaries
* Supports export to TSV and JSON formats
* Computes statistics on feature impacts and variant types

#### Impact Classification System

The Impact Classification System uses multiple criteria to categorize features:

1. **Sequence Identity**: Percent of matching bases between reference and aligned sequences
2. **Coverage**: Percent of the reference feature covered by the alignment
3. **Length Ratio**: Ratio of aligned feature length to reference feature length

These metrics are combined to determine the most appropriate classification:

| Classification | Identity | Coverage | Length Ratio | Description |
|---------------|----------|----------|--------------|-------------|
| Present       | High     | High     | ~1.0         | Feature present with no significant changes |
| Absent        | N/A      | N/A      | N/A          | Feature couldn't be aligned at all |
| Modified      | Medium   | High     | ~1.0         | Feature present with sequence modifications |
| Truncated     | Any      | Any      | <0.9         | Feature present but significantly shorter |
| Expanded      | Any      | Any      | >1.1         | Feature present but significantly longer |
| Fragmented    | Any      | Any      | Any          | Feature split across multiple locations |

The default thresholds are:
- High identity: ≥90%
- High coverage: ≥80%
- Length ratio tolerance: ±10%

These thresholds can be customized through the API.

### Reporting Module (`src/reporting/`)

The Reporting module provides comprehensive functionality to generate reports from analysis results:

#### Report Generator (`src/reporting/report_generator.py`)

* Generates formatted reports from analysis summaries
* Supports multiple output formats (TSV, JSON, RDF)
* Provides filtering by impact type and genomic region
* Handles single-haplotype reports

#### Comparative Report Generator (`src/reporting/comparative_report.py`)

* Generates reports comparing features across multiple haplotypes
* Identifies consensus features present in most haplotypes
* Identifies discriminating features that differentiate haplotypes
* Provides comprehensive comparison metrics

#### Formatters (`src/reporting/formatters/`)

The module includes multiple formatters for different output formats:

* **TSV Formatter**: Generates tab-separated values reports for easy import into spreadsheets
* **JSON Formatter**: Generates structured JSON reports for programmatic processing
* **RDF Formatter**: Generates Resource Description Framework data with support for:
  * Turtle (.ttl)
  * RDF/XML (.rdf)
  * JSON-LD (.jsonld)
  * N-Triples (.nt)

#### RDF Schema (`src/reporting/schemas/shex_schema.shex`)

* ShEx (Shape Expressions) schema that defines the structure of RDF data
* Used for validation of generated RDF reports
* Follows semantic web best practices
* Incorporates standard ontologies like FALDO for genomic positions

#### Report Filtering

The reporting module supports filtering to focus on relevant data:

* **Impact Type Filtering**: Focus on specific impact types (PRESENT, MODIFIED, etc.)
* **Region Filtering**: Extract features within a specific genomic region
* **Combined Filtering**: Apply multiple filters simultaneously

#### RDF Data Model

The RDF data model represents genomic features and their relationships:

* **Feature Representation**: Uses the FALDO ontology for genomic positions
* **Variant Representation**: Structured representation of sequence variants
* **Path Relationships**: Connects features to their containing paths
* **Impact Classification**: Categorizes feature impacts with metrics
* **Comparative Relationships**: Links feature occurrences across different paths

### Parsers Module (`src/parsers/`)

The application uses a dedicated parsers module to handle various input file formats:

#### GFA Parser (`src/parsers/gfa_parser.py`)

* Uses the PyGFA library to parse Graphical Fragment Assembly (GFA) files
* Handles both GFA1 and GFA2 formats
* Provides methods to access segments, paths, and sequences

#### GFF3 Parser (`src/parsers/gff_parser.py`)

* Uses Biopython's BCBio.GFF module to parse GFF3 annotation files
* Indexes features by ID and type for quick lookup
* Provides methods to retrieve features by ID or type

#### FASTA Parser (`src/parsers/fasta_parser.py`)

* Uses Biopython's Bio.SeqIO module to parse FASTA sequence files
* Provides methods to retrieve sequences by ID
* Includes utility methods for sequence manipulation

#### Feature Relationship Graph (`src/parsers/feature_graph.py`)

* Creates a directed graph of parent-child relationships between features
* Uses NetworkX for the graph structure
* Provides methods to access ancestors, descendants, and orphaned features
* Handles the hierarchical nature of genomic features (e.g., gene → mRNA → exon)

## Code Organization

The codebase is organized into the following structure:

```
src/
├── config.py          # Configuration management
├── path_analysis.py   # Path and haplotype analysis
├── parsers/           # Input file parsers
│   ├── __init__.py
│   ├── gfa_parser.py  # GFA file parser
│   ├── gff_parser.py  # GFF3 file parser
│   ├── fasta_parser.py  # FASTA file parser
│   └── feature_graph.py  # Feature relationship graph
├── analysis/          # Feature analysis components
│   ├── __init__.py
│   ├── impact_classifier.py  # Feature impact classification
│   ├── variant_detector.py   # Sequence variant detection
│   ├── reconciliation.py     # Feature hierarchy reconciliation
│   └── summary_generator.py  # Analysis summary generation
├── reporting/         # Reporting components
│   ├── __init__.py
│   ├── report_generator.py     # Single-path report generation
│   ├── comparative_report.py   # Multi-path comparative reports
│   ├── formatters/             # Output format implementations
│   │   ├── __init__.py
│   │   ├── tsv_formatter.py    # TSV output format
│   │   ├── json_formatter.py   # JSON output format
│   │   └── rdf_formatter.py    # RDF output in multiple serializations
│   └── schemas/                # Validation schemas
│       └── shex_schema.shex    # ShEx schema for RDF validation
└── ...
```

## Testing

Each module has corresponding unit tests in the `tests/` directory:

* `tests/test_config.py` - Tests for configuration management
* `tests/test_path_analysis.py` - Tests for path/haplotype analysis
* `tests/test_parsers/` - Tests for file parsers
* `tests/test_analysis/` - Tests for analysis components
* `tests/test_reporting.py` - Tests for reporting components
* `tests/expected_outputs/` - Expected output files for validation

Run the tests using:

```bash
python -m unittest discover tests
```

To run a specific test file:

```bash
python -m unittest tests.test_reporting
```

## RDF Report Validation

The RDF reports generated by the tool can be validated against the provided ShEx schema:

```bash
# Install pyshex if not already installed
pip install pyshex

# Run the validation
python -m pyshex.validate -s src/reporting/schemas/shex_schema.shex -n http://example.org/haplo/AnnotationReport report.ttl
```

## Contributing

(Guidelines for contributing to the project will go here.)

