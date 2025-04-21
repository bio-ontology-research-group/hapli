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

(Detailed explanation of the algorithms, data structures, and design choices used in the *main annotation tool* will go here.)


## Contributing

(Guidelines for contributing to the project will go here.)

