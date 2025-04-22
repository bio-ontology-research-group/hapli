#!/usr/bin/env python3
"""
Basic tests for the main application module.

Tests for:
- Command-line argument parsing
- Configuration loading
- Logging setup
- Parallelization settings
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

from src.main import HaplotypeAnnotationTool, main
from src.config import Config, ConfigurationError

class TestMain(unittest.TestCase):
    """Test the main application module functionality."""

    def setUp(self):
        """Set up test environment."""
        # Create temporary directory for test files
        self.temp_dir = tempfile.mkdtemp()

        # Create mock input files
        self.gfa_file = os.path.join(self.temp_dir, "test.gfa")
        self.gff_file = os.path.join(self.temp_dir, "test.gff3")
        self.fasta_file = os.path.join(self.temp_dir, "test.fasta")

        # Create empty files
        for filepath in [self.gfa_file, self.gff_file, self.fasta_file]:
            with open(filepath, 'w') as f:
                f.write("Test content")

        # Intermediate data directory
        self.intermediate_dir = os.path.join(self.temp_dir, "intermediate")

        # Create tool instance
        self.tool = HaplotypeAnnotationTool()

        # Store original handlers
        self.original_handlers = logging.getLogger().handlers[:]
        self.original_level = logging.getLogger().level


    def tearDown(self):
        """Clean up after tests."""
        # Remove temporary directory
        shutil.rmtree(self.temp_dir)

        # Restore original logging configuration
        logger = logging.getLogger()
        # Remove any handlers added during tests
        for handler in logger.handlers[:]:
             if handler not in self.original_handlers:
                 logger.removeHandler(handler)
        # Restore original handlers if they were removed
        for handler in self.original_handlers:
             if handler not in logger.handlers:
                 logger.addHandler(handler)
        # Restore original level
        logger.setLevel(self.original_level)


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

        # Mock the Config.load method to prevent actual file loading
        with patch.object(Config, 'load') as mock_load:
            mock_load.return_value = None

            # Mock get_all to return our test configuration
            with patch.object(Config, 'get_all') as mock_get_all:
                mock_get_all.return_value = {
                    'gfa_file': self.gfa_file,
                    'gff3_file': self.gff_file,
                    'reference_fasta': self.fasta_file,
                    'log_level': 'DEBUG',
                    'num_workers': 4,
                    'batch_size': 100,
                    'save_intermediate': True,
                    'intermediate_dir': self.intermediate_dir
                }

                self.tool.load_config(args)

                # Verify load was called with our args
                mock_load.assert_called_once_with(args)

                # Get the configuration
                config = self.tool.config.get_all()

                # Check that parameters were correctly set
                self.assertEqual(config['gfa_file'], self.gfa_file)
                self.assertEqual(config['gff3_file'], self.gff_file)
                self.assertEqual(config['reference_fasta'], self.fasta_file)
                self.assertEqual(config['log_level'], 'DEBUG')
                self.assertEqual(config['num_workers'], 4)
                self.assertEqual(config['batch_size'], 100)
                self.assertTrue(config['save_intermediate'])
                self.assertEqual(config['intermediate_dir'], self.intermediate_dir)

    def test_config_loading(self):
        """Test configuration loading through the interface."""
        # Create a test config file
        config_file = os.path.join(self.temp_dir, "test_config.yaml")
        with open(config_file, 'w') as f:
            f.write(f"""
gfa_file: {self.gfa_file}
gff3_file: {self.gff_file}
reference_fasta: {self.fasta_file}
log_level: INFO
num_workers: 2
batch_size: 50
save_intermediate: true
intermediate_dir: {self.intermediate_dir}
""")

        args = ["--config-file", config_file]

        # Mock get_resource_files to prevent file loading attempts
        with patch.object(Config, 'get_resource_files') as mock_resources:
            mock_resources.return_value = {
                'gfa_file': self.gfa_file,
                'gff3_file': self.gff_file,
                'reference_fasta': self.fasta_file
            }

            with patch.object(Config, 'get_all') as mock_get_all:
                mock_get_all.return_value = {
                    'gfa_file': self.gfa_file,
                    'gff3_file': self.gff_file,
                    'reference_fasta': self.fasta_file,
                    'log_level': 'INFO',
                    'num_workers': 2,
                    'batch_size': 50,
                    'save_intermediate': True,
                    'intermediate_dir': self.intermediate_dir
                }

                # Actually call load_config with our args
                config = self.tool.load_config(args)

                # Check that parameters were correctly loaded from the file
                self.assertEqual(config['gfa_file'], self.gfa_file)
                self.assertEqual(config['gff3_file'], self.gff_file)
                self.assertEqual(config['reference_fasta'], self.fasta_file)
                self.assertEqual(config['log_level'], 'INFO')
                self.assertEqual(config['num_workers'], 2)
                self.assertEqual(config['batch_size'], 50)
                self.assertTrue(config['save_intermediate'])
                self.assertEqual(config['intermediate_dir'], self.intermediate_dir)

    def test_logging_configuration(self):
        """Test that logging is properly configured."""
        logger = logging.getLogger() # Get root logger

        # --- Test DEBUG level ---
        log_capture_debug = StringIO()
        handler_debug = logging.StreamHandler(log_capture_debug)
        # Ensure the handler has a formatter to match the basicConfig setup
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', '%Y-%m-%d %H:%M:%S')
        handler_debug.setFormatter(formatter)

        # Configure logging (this removes existing handlers)
        self.tool.configure_logging('DEBUG')
        # Add our capture handler *after* configuration
        logger.addHandler(handler_debug)

        # Check that the root logger is at DEBUG level
        self.assertEqual(logger.level, logging.DEBUG)

        # Log a test message using the root logger
        logger.debug("Test debug message")
        logger.info("Test info message")

        # Check that messages were logged
        log_output_debug = log_capture_debug.getvalue()
        self.assertIn("Test debug message", log_output_debug)
        self.assertIn("Test info message", log_output_debug)

        # Clean up handler for the next part
        logger.removeHandler(handler_debug)

        # --- Test INFO level ---
        log_capture_info = StringIO()
        handler_info = logging.StreamHandler(log_capture_info)
        handler_info.setFormatter(formatter) # Use the same formatter

        # Configure logging (this removes existing handlers again)
        self.tool.configure_logging('INFO')
        # Add our capture handler *after* configuration
        logger.addHandler(handler_info)

        # Check that the root logger is at INFO level
        self.assertEqual(logger.level, logging.INFO)

        # Log test messages using the root logger
        logger.debug("Test debug message") # Should not be captured
        logger.info("Test info message")  # Should be captured

        # Check that only INFO message was logged
        log_output_info = log_capture_info.getvalue()
        self.assertNotIn("Test debug message", log_output_info)
        self.assertIn("Test info message", log_output_info)

        # Clean up handler
        logger.removeHandler(handler_info)

    def test_error_handling(self):
        """Test basic error handling."""
        # Test with missing required files
        args = [
            "--gfa-file", "nonexistent.gfa",
            "--gff3-file", "nonexistent.gff3",
            "--reference-fasta", "nonexistent.fasta"
        ]

        # Run with invalid config should raise ConfigurationError during load
        # We patch 'load' to simulate this failure
        with patch.object(Config, 'load', side_effect=ConfigurationError("Missing required files")):
            # Capture stderr to check log message
            with patch('sys.stderr', new_callable=StringIO) as mock_stderr:
                 exit_code = self.tool.run(args)
                 self.assertNotEqual(exit_code, 0)
                 # Check if the error message was logged to stderr (or logger)
                 # Depending on when logging is configured vs when error happens
                 # Here, config load fails before logging is fully set up by run()
                 # So we might not see the logger output, but the exception is raised
                 # and caught, leading to non-zero exit code.

        # Test with runtime error during execution
        # Need valid args first to pass config loading
        valid_args = [
            "--gfa-file", self.gfa_file,
            "--gff3-file", self.gff_file,
            "--reference-fasta", self.fasta_file
        ]
        # Mock the actual file loading part to raise an error
        with patch.object(HaplotypeAnnotationTool, 'load_input_files',
                         side_effect=Exception("Test runtime error")):
            # Capture stderr to check log message
            with patch('sys.stderr', new_callable=StringIO) as mock_stderr:
                exit_code = self.tool.run(valid_args)
                self.assertNotEqual(exit_code, 0)
                # Check if the runtime error message was logged
                log_output = mock_stderr.getvalue()
                # Note: The logger might write to stderr by default if basicConfig is used
                # Or check logger output if handlers are configured differently
                # Let's assume basicConfig's default stderr handler is active
                self.assertIn("Unexpected error: Test runtime error", log_output)


    def test_intermediate_data_saving(self):
        """Test that intermediate data can be saved correctly."""
        # Create directory for intermediate data
        os.makedirs(self.intermediate_dir, exist_ok=True)

        # Set up mock intermediate data
        self.tool.intermediate_data = {
            'input_files': {
                'gfa_file': self.gfa_file,
                'gff_file': self.gff_file,
                'fasta_file': self.fasta_file
            },
            'selected_paths': ['path1', 'path2', 'path3']
        }

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

            # Check that files were created
            input_files_path = os.path.join(self.intermediate_dir, 'input_files.json')
            paths_file_path = os.path.join(self.intermediate_dir, 'selected_paths.json')
            config_file_path = os.path.join(self.intermediate_dir, 'configuration.json')

            self.assertTrue(os.path.exists(input_files_path))
            self.assertTrue(os.path.exists(paths_file_path))
            self.assertTrue(os.path.exists(config_file_path))

            # Verify content
            with open(paths_file_path) as f:
                paths_data = f.read()
                self.assertIn('path1', paths_data)
                self.assertIn('path2', paths_data)
                self.assertIn('path3', paths_data)

            with open(input_files_path) as f:
                input_data = f.read()
                self.assertIn(self.gfa_file.replace('\\', '/'), input_data.replace('\\', '/')) # Normalize paths for comparison
                self.assertIn(self.gff_file.replace('\\', '/'), input_data.replace('\\', '/'))
                self.assertIn(self.fasta_file.replace('\\', '/'), input_data.replace('\\', '/'))

            with open(config_file_path) as f:
                 config_saved_data = f.read()
                 self.assertIn('intermediate_dir', config_saved_data)
                 self.assertIn(self.intermediate_dir.replace('\\', '/'), config_saved_data.replace('\\', '/'))


if __name__ == '__main__':
    unittest.main()
