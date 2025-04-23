# Developer Guide

## Implementation Details

### Parallelization Module (`src/parallel/`)

The Parallelization module provides a flexible framework for concurrent execution:

#### Task Manager (`src/parallel/task_manager.py`)

* Implements both process-based and thread-based worker pools
* Includes smart task chunking for efficient workload distribution
* Provides progress tracking with ETA estimation
* Ensures proper resource management and cleanup

**Key Components:**
- `TaskChunker`: Divides tasks into balanced chunks
- `ProgressTracker`: Monitors and reports execution progress
- `ProcessWorkerPool`: Implements multiprocessing-based parallelism
- `ThreadWorkerPool`: Implements threading-based parallelism
- `execute_parallel()`: Convenience function for simple parallel execution

**Usage Example:**
```python
from src.parallel.task_manager import execute_parallel

# Process-based parallelism (best for CPU-bound tasks)
results = execute_parallel(
    func=my_function,
    tasks=data_items,
    num_workers=4,  # Number of processes
    pool_type='process'
)

# Thread-based parallelism (best for I/O-bound tasks)
results = execute_parallel(
    func=my_function,
    tasks=data_items,
    num_workers=8,  # Number of threads
    pool_type='thread'
)
```

#### Parallel Alignment (`src/parallel/parallel_alignment.py`)

* Parallelizes minimap2 executions across multiple processes
* Implements smart batching to reduce overhead
* Handles result aggregation from parallel workers
* Provides progress tracking during alignment

**Key Components:**
- `BatchAlignmentTask`: Encapsulates data for independent alignment tasks
- `ParallelAligner`: High-level interface for parallel feature alignment

**Usage Example:**
```python
from src.parallel.parallel_alignment import ParallelAligner

aligner = ParallelAligner(
    num_workers=4,       # Number of parallel workers
    batch_size=50,       # Features per batch
    pool_type='process', # Use processes for computation-heavy alignment
    minimap_preset='splice'  # Passed to minimap2
)

aligned_features = aligner.align_features(features, reference_seq)
```

#### Hierarchical Executor (`src/parallel/hierarchical_executor.py`)

* Executes tasks respecting hierarchical dependencies
* Uses a directed acyclic graph (DAG) to manage execution order
* Optimizes by running independent branches in parallel
* Handles synchronization between dependent tasks

**Key Components:**
- `Task`: Represents a task with dependencies in the execution graph
- `HierarchicalExecutor`: Manages task dependencies and parallel execution
- `execute_hierarchical_tasks()`: Convenience function for task execution

**Usage Example:**
```python
from src.parallel.hierarchical_executor import HierarchicalExecutor

# Create executor
executor = HierarchicalExecutor(num_workers=4)

# Add tasks with dependencies
executor.add_task('load_data', load_function)
executor.add_task('parse', parse_function, dependencies=['load_data'])
executor.add_task('filter', filter_function, dependencies=['parse'])
executor.add_task('analyze', analyze_function, dependencies=['filter'])

# Execute respecting dependencies
results = executor.execute()
```

### Main Application Module (`src/main.py`)

The main application module provides the command-line interface and orchestrates the entire annotation pipeline:

#### HaplotypeAnnotationTool Class

The core class that implements the annotation pipeline:

* **Initialization and Setup:**
  * `__init__()`: Initializes the tool components
  * `configure_logging()`: Sets up logging with the specified verbosity
  * `load_config()`: Loads configuration from command-line args and/or config file
  * `initialize_components()`: Creates parsers, analyzers, and processors

* **Data Processing:**
  * `load_input_files()`: Loads GFA, GFF3, and FASTA files
  * `select_paths()`: Selects paths to annotate based on config criteria
  * `save_intermediate_data()`: Saves intermediate results for visualization

* **Pipeline Execution:**
  * `run()`: Runs the full annotation pipeline
  * Handles errors and returns appropriate exit codes

#### Main Entry Point

* `main()`: Command-line entry point function that creates a tool instance and runs it
* Returns exit code (0 for success, non-zero for error)

#### Dependency Management

* The module follows dependency injection patterns for better testability
* Components are initialized with appropriate parallelization settings
* Configuration is centralized and validated before pipeline execution

#### Error Handling

* User-friendly error messages with appropriate logging levels
* Graceful handling of configuration errors vs. runtime errors
* Detailed error information in debug mode

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
*   **Main Module Integration:** The configuration is loaded by the HaplotypeAnnotationTool and used to control:
    1.  Logging verbosity
    2.  Input file paths
    3.  Path/haplotype selection
    4.  Parallelization settings
    5.  Output formatting
    6.  Intermediate data saving

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

### Alignment Module (`src/alignment/`)

The Alignment module provides functionality for sequence alignment using minimap2:

#### Minimap2 Wrapper (`src/alignment/minimap_wrapper.py`)

* Provides a Python interface to the minimap2 aligner via the mappy library
* Configures alignment parameters for different sequence types
* Handles various input formats (strings, SeqRecords, files)
* Processes alignment results into a standardized format

**Key Components:**
- `MinimapAligner`: Main class for minimap2 alignment
- `AlignmentResult`: Data class for storing alignment results

**Alignment Parameters:**
- **Presets:** Predefined parameter sets for different sequence types
  - `map-pb`/`map-ont`: PacBio/Oxford Nanopore genomic reads
  - `map-hifi`: PacBio HiFi genomic reads
  - `sr`: Short genomic paired-end reads
  - `splice`: Long-read spliced alignment
  - `asm5`/`asm10`/`asm20`: Assembly to reference alignment
  - `cdna`: cDNA alignment
- **Custom Parameters:**
  - `k`: k-mer size (11-15 for short reads, 15-19 for genomic)
  - `w`: minimizer window size
  - `min_intron_len`/`max_intron_len`: intron length bounds for splice mapping
  - `scoring`: tuple of (match, mismatch, gap_open, gap_extend)

**Usage Example:**
```python
from src.alignment.minimap_wrapper import MinimapAligner

# Create an aligner with appropriate preset
aligner = MinimapAligner(preset="map-ont")  # For Oxford Nanopore reads

# Load reference sequence
aligner.load_reference(reference_seq)
# Or load from file
aligner.load_reference_file("reference.fasta")

# Align a query sequence
alignments = aligner.align_sequence(query_seq)

# Process alignment results
for aln in alignments:
    print(f"Alignment score: {aln.score}")
    print(f"Identity: {aln.identity:.2f}")
    print(f"CIGAR: {aln.cigar}")
    print(f"Query region: {aln.query_start}-{aln.query_end}")
    print(f"Target region: {aln.target_start}-{aln.target_end}")
```

#### Alignment Result Processing (`src/alignment/alignment_result.py`)

* Provides standardized data structures for alignment results
* Converts minimap2 alignments to a comprehensive format
* Calculates detailed alignment statistics (identity, coverage, gaps, etc.)
* Supports serialization/deserialization for result storage
* Includes text-based alignment visualization

**Key Components:**
- `AlignmentResult`: Comprehensive representation of an alignment
- `AlignmentStatistics`: Statistical metrics for alignments
- `CigarOperation`: Representation of CIGAR string operations
- `AlignmentType`: Enumeration of alignment quality categories

**Data Structure:**
The `AlignmentResult` class contains:
- Basic alignment information (query/target names, sequences)
- Alignment coordinates (start/end positions)
- Alignment details (score, CIGAR string, strand)
- Computed statistics (identity, coverage, matches, mismatches, gaps)
- Visualization data (aligned sequences with gap characters)

**Usage Example:**
```python
from src.alignment.alignment_result import AlignmentResult

# Create from minimap2 alignment
result = AlignmentResult.from_minimap2(alignment, query_seq, target_seq)

# Access alignment information
print(result.get_summary())
print(f"Identity: {result.statistics.identity:.2%}")
print(f"Coverage: {result.statistics.coverage:.2%}")

# Serialize to JSON
json_data = result.to_json()
with open("alignment_result.json", "w") as f:
    f.write(json_data)

# Deserialize from JSON
with open("alignment_result.json") as f:
    data = json.load(f)
    restored_result = AlignmentResult.from_dict(data)
```

#### Terminal Display (`src/alignment/terminal_display.py`)

* Renders alignments as text-based representations
* Shows sequence matches/mismatches with ASCII characters
* Uses color coding (via ANSI escape codes) for terminal display
* Includes options for compact vs. detailed views

**Key Components:**
- `AlignmentDisplay`: Main class for rendering alignments
- `display_alignment_result()`: Convenience function for single alignment display
- `print_alignment_result()`: Function to directly print an alignment to console

**Display Features:**
- Color-coded sequence representation (green for matches, red for gaps, etc.)
- Match/mismatch indicators (|, ., space)
- Position numbering for reference
- Compact summary tables for multiple alignments
- Detailed statistics display

**Terminal Visualization:**
The alignment visualization shows:
```
Alignment: query1 to target1
Score: 60.0
Identity: 95.00%
Coverage: 100.00%
Matches: 95, Mismatches: 5, Gaps: 0
Query: 0-100 (Forward)
Target: 0-100

      0 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
        ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      0 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT

     60 ACGTACGTACGTACGTACGTACGT
        ||||||||.||||.||||||.|||
     60 ACGTACGTTCGTTAGTACGTTCGT
```

**Usage Example:**
```python
from src.alignment.terminal_display import AlignmentDisplay, print_alignment_result

# Direct printing of an alignment
print_alignment_result(alignment_result, use_color=True, detailed=True)

# More control with the AlignmentDisplay class
display = AlignmentDisplay(use_color=True, line_width=100)

# Display a single alignment
output = display.display_alignment(alignment_result, detailed=True)
print(output)

# Display multiple alignments in compact form
output = display.display_compact_summary(alignment_results)
print(output)
```

**Test Data:**
The test suite uses artificially generated test data based on the PhiX174 viral genome:
- Reference sequence (~5.4kb)
- Extracted regions (100-500bp segments)
- Variant sequences with:
  - SNPs (single nucleotide polymorphisms)
  - Insertions (1-10bp)
  - Deletions (1-10bp)
  - Complex combinations of the above

The test data is automatically generated when running the tests if it doesn't exist.

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
├── parallel/          # Parallel processing components
│   ├── __init__.py
│   ├── task_manager.py         # Worker pool implementation
│   ├── parallel_alignment.py   # Parallel feature alignment
│   └── hierarchical_executor.py  # Dependency-based execution
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

