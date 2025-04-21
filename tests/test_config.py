import unittest
import os
import yaml
import tempfile
from unittest.mock import patch
import sys

# Dynamically adjust path to import from src
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(project_root, 'src'))

from config import Config, ConfigurationError # noqa: E402

class TestConfig(unittest.TestCase):

    def setUp(self):
        """Set up temporary files and directories for tests."""
        self.test_dir = tempfile.TemporaryDirectory()
        self.config_dir = self.test_dir.name

        # Create dummy resource files for existence checks
        self.dummy_gfa = self._create_dummy_file("test.gfa")
        self.dummy_gff = self._create_dummy_file("test.gff3")
        self.dummy_fasta = self._create_dummy_file("test.fasta")

        # Create a valid config file
        self.valid_config_path = os.path.join(self.config_dir, "valid_config.yaml")
        self.valid_config_data = {
            "gfa_file": self.dummy_gfa,
            "gff3_file": self.dummy_gff,
            "reference_fasta": self.dummy_fasta,
            "output_file": "output.txt",
            "log_level": "DEBUG",
        }
        with open(self.valid_config_path, 'w') as f:
            yaml.dump(self.valid_config_data, f)

        # Create an invalid config file (missing required param)
        self.invalid_config_path = os.path.join(self.config_dir, "invalid_config.yaml")
        self.invalid_config_data = {
            "gfa_file": self.dummy_gfa,
            # "gff3_file": "missing", # Required param missing
            "reference_fasta": self.dummy_fasta,
        }
        with open(self.invalid_config_path, 'w') as f:
            yaml.dump(self.invalid_config_data, f)

        # Create config file with non-existent file paths
        self.nonexistent_files_config_path = os.path.join(self.config_dir, "nonexistent_files.yaml")
        self.nonexistent_files_data = {
             "gfa_file": os.path.join(self.config_dir, "non_existent.gfa"),
             "gff3_file": self.dummy_gff,
             "reference_fasta": self.dummy_fasta,
        }
        with open(self.nonexistent_files_config_path, 'w') as f:
            yaml.dump(self.nonexistent_files_data, f)


    def tearDown(self):
        """Clean up temporary directory."""
        self.test_dir.cleanup()
        # Remove src from sys.path if it was added
        if os.path.join(project_root, 'src') in sys.path:
            sys.path.remove(os.path.join(project_root, 'src'))


    def _create_dummy_file(self, filename: str) -> str:
        """Helper to create an empty dummy file."""
        path = os.path.join(self.config_dir, filename)
        with open(path, 'w') as f:
            f.write("") # Create empty file
        return path

    def test_load_valid_config_file(self):
        """Tests loading settings solely from a valid YAML file."""
        config = Config()
        args = ["--config-file", self.valid_config_path]
        config.load(args)

        self.assertEqual(config.get("gfa_file"), self.dummy_gfa)
        self.assertEqual(config.get("gff3_file"), self.dummy_gff)
        self.assertEqual(config.get("reference_fasta"), self.dummy_fasta)
        self.assertEqual(config.get("output_file"), "output.txt")
        self.assertEqual(config.get("log_level"), "DEBUG") # Overrides default

        resources = config.get_resource_files()
        self.assertEqual(resources["gfa_file"], self.dummy_gfa)
        self.assertEqual(resources["gff3_file"], self.dummy_gff)
        self.assertEqual(resources["reference_fasta"], self.dummy_fasta)


    def test_missing_required_parameter_in_file(self):
        """Tests that loading raises an error if a required parameter is missing."""
        config = Config()
        args = ["--config-file", self.invalid_config_path]
        with self.assertRaisesRegex(ConfigurationError, "Missing required configuration parameters: gff3_file"):
            config.load(args)

    def test_missing_required_parameter_no_file_no_cli(self):
        """Tests error when required params are not provided at all."""
        config = Config()
        args = [] # No config file, no CLI args for required params
        with self.assertRaisesRegex(ConfigurationError, "Missing required configuration parameters: gfa_file, gff3_file, reference_fasta"):
            config.load(args)


    def test_cli_overrides_file(self):
        """Tests that CLI arguments override config file settings."""
        config = Config()
        cli_output = os.path.join(self.config_dir, "cli_output.txt")
        cli_gfa = self._create_dummy_file("cli.gfa") # Need another dummy file that exists

        args = [
            "--config-file", self.valid_config_path,
            "--log-level", "WARNING", # Override file's DEBUG
            "--output-file", cli_output, # Override file's output.txt
            "--gfa-file", cli_gfa # Override file's test.gfa
        ]
        config.load(args)

        self.assertEqual(config.get("log_level"), "WARNING")
        self.assertEqual(config.get("output_file"), cli_output)
        self.assertEqual(config.get("gfa_file"), cli_gfa)
        # Check that unspecified params still come from the file
        self.assertEqual(config.get("gff3_file"), self.dummy_gff)
        self.assertEqual(config.get("reference_fasta"), self.dummy_fasta)

        resources = config.get_resource_files()
        self.assertEqual(resources["gfa_file"], cli_gfa) # Check resource registration updated


    def test_cli_only_config(self):
        """Tests providing all required config via CLI."""
        config = Config()
        cli_output = "cli_only_out.txt"
        args = [
            "--gfa-file", self.dummy_gfa,
            "--gff3-file", self.dummy_gff,
            "--reference-fasta", self.dummy_fasta,
             "--output-file", cli_output,
             "--log-level", "INFO" # Uses defaults if not set, check required only
        ]
        config.load(args)

        self.assertEqual(config.get("gfa_file"), self.dummy_gfa)
        self.assertEqual(config.get("gff3_file"), self.dummy_gff)
        self.assertEqual(config.get("reference_fasta"), self.dummy_fasta)
        self.assertEqual(config.get("output_file"), cli_output)
        # Check default value if not provided by CLI either
        self.assertEqual(config.get("log_level"), "INFO") # Should be default


    def test_non_existent_config_file(self):
         """Tests error when specified config file does not exist."""
         config = Config()
         args = ["--config-file", os.path.join(self.config_dir, "non_existent_config.yaml")]
         with self.assertRaisesRegex(ConfigurationError, "Config file not found"):
             config.load(args)


    def test_non_existent_resource_file(self):
        """Tests error when a required resource file specified in config does not exist."""
        config = Config()
        args = ["--config-file", self.nonexistent_files_config_path]
        # The error message checks the specific parameter name and the path
        expected_error_regex = r"Required file not found: gfa_file = .*non_existent.gfa"
        with self.assertRaisesRegex(ConfigurationError, expected_error_regex):
            config.load(args)


if __name__ == "__main__":
    unittest.main()
