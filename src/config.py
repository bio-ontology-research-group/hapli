import argparse
import yaml
import os
import sys # Import the sys module
from typing import Dict, Any, Optional, List

# Default configuration values
DEFAULT_CONFIG: Dict[str, Any] = {
    "log_level": "INFO",
    "output_file": None,
}

# Required configuration parameters
REQUIRED_PARAMS: List[str] = ["gfa_file", "gff3_file", "reference_fasta"]

class ConfigurationError(Exception):
    """Custom exception for configuration errors."""
    pass

class Config:
    """
    Manages configuration settings for the application.

    Loads settings from YAML files and allows overrides via command-line arguments.
    Validates the presence of required parameters.
    """
    def __init__(self):
        self._settings: Dict[str, Any] = DEFAULT_CONFIG.copy()
        self._resource_files: Dict[str, Optional[str]] = {
            "gfa_file": None,
            "gff3_file": None,
            "reference_fasta": None
        }

    def load(self, args: Optional[List[str]] = None):
        """
        Loads configuration from files and CLI arguments.

        Args:
            args: Optional list of command-line arguments (primarily for testing).
                  If None, sys.argv[1:] will be used.
        """
        parser = self._setup_argparser()
        # Use sys.argv[1:] if args is None, otherwise use provided args
        cli_args = args if args is not None else sys.argv[1:]
        parsed_args = parser.parse_args(cli_args)


        # 1. Load from config file if specified
        if parsed_args.config_file:
            if not os.path.exists(parsed_args.config_file):
                raise ConfigurationError(f"Config file not found: {parsed_args.config_file}")
            try:
                with open(parsed_args.config_file, 'r') as f:
                    file_config = yaml.safe_load(f)
                    if file_config:  # Check if file is not empty
                         self._settings.update(file_config)
            except yaml.YAMLError as e:
                raise ConfigurationError(f"Error parsing config file {parsed_args.config_file}: {e}")
            except Exception as e:
                 raise ConfigurationError(f"Error reading config file {parsed_args.config_file}: {e}")


        # 2. Override with CLI arguments (only those explicitly provided)
        cli_overrides = {
            # Use cli_args for iterating over what was actually passed
            key: value for key, value in vars(parsed_args).items()
            # Check if the argument was actually provided in cli_args or has a default value that isn't None
            if value is not None and key != 'config_file'
        }
        self._settings.update(cli_overrides)


        # 3. Validate required parameters
        self._validate_required()

        # 4. Register resource files
        self._register_resources()

    def _setup_argparser(self) -> argparse.ArgumentParser:
        """Sets up the argument parser."""
        parser = argparse.ArgumentParser(description="Haplotype Annotation Tool Configuration")

        # Configuration file option
        parser.add_argument(
            "-c", "--config-file",
            type=str,
            help="Path to the YAML configuration file."
        )

        # Resource files (can be set in config file or CLI)
        parser.add_argument(
            "--gfa-file",
            type=str,
            help="Path to the GFA file."
        )
        parser.add_argument(
            "--gff3-file",
            type=str,
            help="Path to the GFF3 annotation file."
        )
        parser.add_argument(
            "--reference-fasta",
            type=str,
            help="Path to the reference FASTA file."
        )

        # Other options (can be set in config file or CLI)
        parser.add_argument(
            "--output-file",
            type=str,
            help="Path to the output annotation file."
        )
        parser.add_argument(
            "--log-level",
            type=str,
            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            # Default is now handled by DEFAULT_CONFIG, avoid setting default here
            # to distinguish between explicitly set and default values.
            help="Logging level."
        )

        return parser

    def _validate_required(self):
        """Checks if all required parameters are set."""
        missing = [param for param in REQUIRED_PARAMS if self._settings.get(param) is None]
        if missing:
            raise ConfigurationError(f"Missing required configuration parameters: {', '.join(missing)}")

        # Additional validation: Check if files exist
        for param in REQUIRED_PARAMS:
             filepath = self._settings.get(param)
             # Make sure filepath is not None before checking existence
             if filepath and not os.path.exists(filepath):
                 # Warning might be better than error, depending on tool stage
                 # For now, raise error if required file doesn't exist after config load
                 raise ConfigurationError(f"Required file not found: {param} = {filepath}")


    def _register_resources(self):
        """Updates the resource file paths based on the loaded configuration."""
        for key in self._resource_files:
            if key in self._settings:
                self._resource_files[key] = self._settings[key]

    def get(self, key: str, default: Any = None) -> Any:
        """Retrieves a configuration value."""
        # Use the internal settings dictionary, falling back to the provided default
        return self._settings.get(key, default)

    def get_all(self) -> Dict[str, Any]:
        """Retrieves all configuration settings."""
        return self._settings.copy()

    def get_resource_files(self) -> Dict[str, Optional[str]]:
        """Retrieves the registered resource file paths."""
        return self._resource_files.copy()

# Singleton instance (optional, but common pattern)
# config = Config()

# Example Usage (if run directly, though typically imported)
if __name__ == "__main__":
    try:
        config_manager = Config()
        # Pass sys.argv explicitly to load method for clarity
        config_manager.load(sys.argv[1:]) # Loads from arguments passed to the script
        print("Configuration loaded successfully:")
        print(config_manager.get_all())
        print("\nResource files:")
        print(config_manager.get_resource_files())

        # Access specific settings
        gfa = config_manager.get("gfa_file")
        print(f"\nGFA File: {gfa}")

    except ConfigurationError as e:
        # Use sys.stderr correctly now that sys is imported
        print(f"Configuration Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)
