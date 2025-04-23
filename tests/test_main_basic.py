#!/usr/bin/env python3
"""
Basic tests for the main application module.

Tests for:
- Command-line argument parsing
- Configuration loading
- Logging setup
- Error handling
- Intermediate data saving
"""

import os
import sys
import unittest
import tempfile
import shutil
import logging
from unittest.mock import patch, MagicMock
from io import StringIO
import json # Import json for loading saved data

# Ensure src is in path
script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(script_dir, '..', 'src'))
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

from src.main import HaplotypeAnnotationTool, main
from src.config import Config, ConfigurationError
# Import placeholder summaries if needed for testing save_intermediate
from src.analysis.summary_generator import AnalysisSummary, FeatureSummary
from src.analysis.impact_classifier import ImpactType


class TestMain(unittest.TestCase):
    """Test the main application module functionality."""

    def setUp(self):
        """Set up test environment."""
        # Create temporary directory for test files
        self.temp_dir = tempfile.mkdtemp(prefix="main_test_")

        # Create mock input files
        self.gfa_file = os.path.join(self.temp_dir, "test.gfa")
        self.gff_file = os.path.join(self.temp_dir, "test.gff3")
        self.fasta_file = os.path.join(self.temp_dir, "test.fasta")

        # Create simple content for mock files to allow parsing
        with open(self.gfa_file, 'w') as f:
            f.write("H\tVN:Z:1.0\n")
            f.write("S\t1\tACGT\n")
            f.write("P\tpath1\t1+\t*\n")
        with open(self.gff_file, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("seq1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")
            f.write("seq1\t.\texon\t10\t50\t.\t+\t.\tID=exon1;Parent=gene1\n")
        with open(self.fasta_file, 'w') as f:
            f.write(">seq1\n")
            f.write("ACGT" * 25 + "\n")


        # Intermediate data directory
        self.intermediate_dir = os.path.join(self.temp_dir, "intermediate")

        # Create tool instance
        self.tool = HaplotypeAnnotationTool()

        # Store original handlers and level to restore later
        self.original_handlers = logging.getLogger().handlers[:]
        self.original_level = logging.getLogger().level


    def tearDown(self):
        """Clean up after tests."""
        # Remove temporary directory
        shutil.rmtree(self.temp_dir)

        # Restore original logging configuration
        root_logger = logging.getLogger()
        # Remove any handlers added during tests or by the tool
        for handler in root_logger.handlers[:]:
             root_logger.removeHandler(handler)
        # Restore original handlers
        for handler in self.original_handlers:
             root_logger.addHandler(handler)
        # Restore original level
        root_logger.setLevel(self.original_level)


    def test_command_line_args_parsing(self):
        """Test that command-line arguments are correctly parsed."""
        args = [
            "--gfa-file", self.gfa_file,
            "--gff3-file", self.gff_file,
            "--reference-fasta", self.fasta_file,
            "--log-level", "DEBUG",
            "--num-workers", "4",
            "--batch-size", "100",
            "--save-intermediate",
            "--intermediate-dir", self.intermediate_dir
        ]

        # Mock the Config methods used within load_config
        with patch.object(Config, 'load') as mock_load, \
             patch.object(Config, 'get_all') as mock_get_all:

            # Define what get_all should return after load is called
            mock_get_all.return_value = {
                'gfa_file': self.gfa_file,
                'gff3_file': self.gff_file,
                'reference_fasta': self.fasta_file,
                'log_level': 'DEBUG',
                'num_workers': 4,
                'batch_size': 100,
                'save_intermediate': True,
                'intermediate_dir': self.intermediate_dir,
                # Add other defaults expected by the tool if necessary
                'output_file': None,
                'output_format': 'tsv',
                'comparative': False,
            }

            # Call the method that uses Config.load and Config.get_all
            loaded_config = self.tool.load_config(args)

            # Verify load was called correctly
            mock_load.assert_called_once_with(args)

            # Check the returned config dictionary
            self.assertEqual(loaded_config['gfa_file'], self.gfa_file)
            self.assertEqual(loaded_config['log_level'], 'DEBUG')
            self.assertEqual(loaded_config['num_workers'], 4)
            self.assertTrue(loaded_config['save_intermediate'])
            self.assertEqual(loaded_config['intermediate_dir'], self.intermediate_dir)

    def test_config_loading_from_file(self):
        """Test configuration loading from a YAML file."""
        # Create a test config file
        config_file = os.path.join(self.temp_dir, "test_config.yaml")
        with open(config_file, 'w') as f:
            f.write(f"""
gfa_file: {self.gfa_file}
gff3_file: {self.gff_file}
reference_fasta: {self.fasta_file}
log_level: INFO
num_workers: 2
save_intermediate: true
intermediate_dir: {self.intermediate_dir}
output_file: output.tsv
""")

        args = ["--config-file", config_file]

        # Mock get_all to simulate what Config would return after loading the file
        # This avoids needing a full Config integration test here
        with patch.object(Config, 'get_all') as mock_get_all:
             mock_get_all.return_value = {
                 'gfa_file': self.gfa_file,
                 'gff3_file': self.gff_file,
                 'reference_fasta': self.fasta_file,
                 'log_level': 'INFO',
                 'num_workers': 2,
                 'batch_size': 100, # Assuming default if not in file
                 'save_intermediate': True,
                 'intermediate_dir': self.intermediate_dir,
                 'output_file': 'output.tsv',
                 'output_format': 'tsv', # Assuming default
                 'comparative': False, # Assuming default
                 'config_file': config_file # Config usually stores this too
             }

             # Call load_config, which internally calls config.load -> config.get_all
             config = self.tool.load_config(args)

             # Check that parameters were correctly loaded
             self.assertEqual(config['log_level'], 'INFO')
             self.assertEqual(config['num_workers'], 2)
             self.assertTrue(config['save_intermediate'])
             self.assertEqual(config['intermediate_dir'], self.intermediate_dir)
             self.assertEqual(config['output_file'], 'output.tsv')


    def test_logging_configuration(self):
        """Test that logging is properly configured using a custom log handler."""
        # Create a StringIO object to capture log output
        log_capture = StringIO()
        # Create a handler that writes to the StringIO object
        handler = logging.StreamHandler(log_capture)
        formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
        handler.setFormatter(formatter)
        
        # Add the handler to the root logger to capture all logs
        root_logger = logging.getLogger()
        root_logger.addHandler(handler)
        # Save original level to restore later
        original_level = root_logger.level
        # Ensure the root logger level is set to DEBUG to capture all logs
        root_logger.setLevel(logging.DEBUG)
        
        # Test DEBUG level
        self.tool.configure_logging('DEBUG')
        # Log messages using different loggers to test propagation
        logging.getLogger('src.main').debug("Test main debug")
        logging.getLogger('src.parsers').info("Test parser info")
        logging.getLogger().warning("Test root warning")
        
        # Get the captured log output
        log_output = log_capture.getvalue()
        
        # Check captured output
        self.assertIn("INFO:src.main:Logging configured at DEBUG level", log_output)
        self.assertIn("DEBUG:src.main:Test main debug", log_output)
        self.assertIn("INFO:src.parsers:Test parser info", log_output)
        self.assertIn("WARNING:root:Test root warning", log_output)
        
        # Clear the log capture for the next test
        log_capture.truncate(0)
        log_capture.seek(0)
        
        # Clear the log capture for the next test
        log_capture.truncate(0)
        log_capture.seek(0)
        
        # Reset the root logger for the next test
        root_logger = logging.getLogger()
        # Ensure the handler is still attached
        if handler not in root_logger.handlers:
            root_logger.addHandler(handler)
        # Set appropriate level for INFO test
        root_logger.setLevel(logging.INFO)
        
        # Test INFO level
        self.tool.configure_logging('INFO')
        logging.getLogger('src.main').debug("Test main debug") # Should not be captured at INFO level
        logging.getLogger('src.main').info("Test main info")
        logging.getLogger().error("Test root error")
        
        # Get the captured log output
        log_output = log_capture.getvalue()
        
        # Check captured output
        self.assertIn("INFO:src.main:Logging configured at INFO level", log_output)
        self.assertIn("INFO:src.main:Test main info", log_output)
        self.assertIn("ERROR:root:Test root error", log_output)
        # Verify debug message is NOT present (since we're at INFO level)
        self.assertNotIn("DEBUG:src.main:Test main debug", log_output)
        
        # Clean up
        root_logger = logging.getLogger()
        root_logger.removeHandler(handler)
        # Restore original level
        root_logger.setLevel(original_level)


    def test_error_handling_config(self):
        """Test error handling for configuration errors using assertLogs."""
        args = ["--gfa-file", "nonexistent.gfa"] # Missing other required args

        # Configure logging first to ensure logs are captured
        self.tool.configure_logging('ERROR')
        
        # Expect ConfigurationError during config.load()
        # Capture logs from the specific logger used in the tool
        with self.assertLogs('src.main', level='ERROR') as cm:
             # Mock Config.load to raise the expected error
             with patch.object(Config, 'load', side_effect=ConfigurationError("Missing required files")):
                  exit_code = self.tool.run(args)

        self.assertNotEqual(exit_code, 0)
        # Check if the ConfigurationError message was logged
        self.assertTrue(any("Configuration error: Missing required files" in log for log in cm.output))


    def test_error_handling_runtime(self):
        """Test error handling for runtime errors using assertLogs."""
        # Use valid args to pass config loading
        valid_args = [
            "--gfa-file", self.gfa_file,
            "--gff3-file", self.gff_file,
            "--reference-fasta", self.fasta_file,
            "--log-level", "INFO" # Ensure level is sufficient for CRITICAL
        ]

        # Configure logging first to ensure logs are captured
        self.tool.configure_logging('CRITICAL')
        
        # Mock a later step (e.g., load_input_files) to raise an error
        runtime_error_msg = "Test runtime error during file load"
        with patch.object(HaplotypeAnnotationTool, 'load_input_files',
                         side_effect=Exception(runtime_error_msg)):
            # Capture logs from the specific logger used in the tool
            with self.assertLogs('src.main', level='CRITICAL') as cm:
                exit_code = self.tool.run(valid_args)

        self.assertNotEqual(exit_code, 0)
        # Check if the critical error message was logged
        expected_log = f"A critical unexpected error occurred: {runtime_error_msg}"
        self.assertTrue(any(expected_log in log for log in cm.output))


    def test_intermediate_data_saving(self):
        """Test that intermediate data can be saved correctly."""
        # Create directory for intermediate data
        os.makedirs(self.intermediate_dir, exist_ok=True)

        # Create a StringIO object to capture log output
        log_capture = StringIO()
        # Create a handler that writes to the StringIO object
        handler = logging.StreamHandler(log_capture)
        formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
        handler.setFormatter(formatter)
        
        # Add the handler to the logger
        logger = logging.getLogger('src.main')
        logger.addHandler(handler)
        
        # Configure logging
        self.tool.configure_logging('INFO')
        
        # Set up mock intermediate data and analysis summaries
        self.tool.intermediate_data = {
            'input_files': {
                'gfa_file': self.gfa_file,
                'gff3_file': self.gff_file,
                'reference_fasta': self.fasta_file
            },
            'selected_paths': ['path1']
        }
        # Create a mock AnalysisSummary (matching actual structure)
        mock_summary = AnalysisSummary(
             path_id='path1',
             feature_count=1,
             feature_by_impact={ImpactType.PRESENT.value: 1},
             variant_counts={},
             reconciliation_counts={},
             feature_summaries={
                 'gene1': FeatureSummary(
                     feature_id='gene1',
                     feature_type='gene',
                     impact_type=ImpactType.PRESENT,
                     sequence_identity=1.0,
                     coverage=1.0,
                     variants=[],
                     path_id='path1'
                 )
             }
        )
        self.tool.analysis_summaries = {'path1': mock_summary}
        # Mock features dict needed for saving features_info
        self.tool.features = {'gene1': MagicMock(type='gene')}

        # Mock config - provide necessary file paths for saving config itself
        mock_config_data = {
            'gfa_file': self.gfa_file,
            'gff3_file': self.gff_file,
            'reference_fasta': self.fasta_file,
            'intermediate_dir': self.intermediate_dir,
            'save_intermediate': True # Assume this was set
        }
        
        # Patch the config object within the tool instance
        with patch.object(self.tool.config, 'get_all', return_value=mock_config_data):
            # Save the data
            self.tool.save_intermediate_data(self.intermediate_dir)
        
        # Get the captured log output
        log_output = log_capture.getvalue()
        
        # Check log message - verify it contains the expected text
        self.assertIn("Successfully saved intermediate data files to", log_output)
        self.assertIn(self.intermediate_dir, log_output)
        
        # Clean up
        logger.removeHandler(handler)

        # Check that files were created
        input_files_path = os.path.join(self.intermediate_dir, 'input_files.json')
        paths_file_path = os.path.join(self.intermediate_dir, 'selected_paths.json')
        config_file_path = os.path.join(self.intermediate_dir, 'configuration.json')
        features_info_path = os.path.join(self.intermediate_dir, 'parsed_features_info.json')
        summaries_path = os.path.join(self.intermediate_dir, 'analysis_summaries.json')

        self.assertTrue(os.path.exists(input_files_path))
        self.assertTrue(os.path.exists(paths_file_path))
        self.assertTrue(os.path.exists(config_file_path))
        self.assertTrue(os.path.exists(features_info_path))
        self.assertTrue(os.path.exists(summaries_path))

        # Verify content (basic checks)
        with open(paths_file_path) as f:
            paths_data = json.load(f)
            self.assertEqual(paths_data, ['path1'])

        with open(input_files_path) as f:
            input_data = json.load(f)
            self.assertEqual(input_data.get('gfa_file'), self.gfa_file)

        with open(config_file_path) as f:
             config_saved_data = json.load(f)
             self.assertEqual(config_saved_data.get('intermediate_dir'), self.intermediate_dir)

        with open(features_info_path) as f:
             features_data = json.load(f)
             self.assertEqual(features_data, {'gene1': 'gene'})

        with open(summaries_path) as f:
             summaries_data = json.load(f)
             self.assertIn('path1', summaries_data)
             self.assertEqual(summaries_data['path1']['path_id'], 'path1')
             self.assertIn('feature_summaries', summaries_data['path1'])
             self.assertIn('gene1', summaries_data['path1']['feature_summaries'])
             # Use .name for Enum comparison after loading from JSON
             self.assertEqual(summaries_data['path1']['feature_summaries']['gene1']['impact_type'], ImpactType.PRESENT.name)


if __name__ == '__main__':
    unittest.main()
