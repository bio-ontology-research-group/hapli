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
├── parsers/           # Input file parsers
│   ├── __init__.py
│   ├── gfa_parser.py  # GFA file parser
│   ├── gff_parser.py  # GFF3 file parser
│   ├── fasta_parser.py  # FASTA file parser
│   └── feature_graph.py  # Feature relationship graph
├── ...
```

## Testing

Each module has corresponding unit tests in the `tests/` directory:

* `tests/test_config.py` - Tests for configuration management
* `tests/test_gfa_parser.py` - Tests for GFA parsing
* `tests/test_gff_parser.py` - Tests for GFF3 parsing
* `tests/test_fasta_parser.py` - Tests for FASTA parsing
* `tests/test_feature_graph.py` - Tests for feature relationship graphs

Run the tests using:

```bash
python -m unittest discover tests
```

## Contributing

(Guidelines for contributing to the project will go here.)

