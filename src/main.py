#!/usr/bin/env python3
"""
Main application entry point for the Haplotype Annotation Tool.

This module implements the command-line interface and orchestrates the
entire annotation pipeline, from loading input files to generating reports.
"""

import os
import sys
import logging
import argparse
import json
import time
from typing import Dict, List, Any, Optional, Tuple

from src.config import Config, ConfigurationError
from src.parallel.task_manager import create_worker_pool, execute_parallel
from src.parallel.hierarchical_executor import HierarchicalExecutor

# Configure logger
logger = logging.getLogger(__name__)

class HaplotypeAnnotationTool:
    """
    Main class for the Haplotype Annotation Tool.
    
    Orchestrates the entire annotation pipeline, including:
    - Loading and validating input files
    - Selecting paths/haplotypes to analyze
    - Aligning features to selected paths
    - Analyzing feature impacts
    - Generating reports
    
    Supports parallel execution for performance optimization.
    """
    
    def __init__(self):
        """Initialize the annotation tool."""
        self.config = Config()
        self.path_analyzer = None
        self.feature_graph = None
        self.alignment_processor = None
        self.impact_classifier = None
        self.summary_generator = None
        self.report_generator = None
        self.intermediate_data = {}
        
    def configure_logging(self, log_level: str) -> None:
        """
        Configure logging with the specified verbosity level.
        
        Args:
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        numeric_level = getattr(logging, log_level.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f"Invalid log level: {log_level}")
        
        # Reset existing handlers to avoid duplicate log messages
        root_logger = logging.getLogger()
        for handler in root_logger.handlers[:]:
            root_logger.removeHandler(handler)
            
        # Configure root logger
        logging.basicConfig(
            level=numeric_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # Make sure the level is directly set on the root logger
        # (basicConfig might be ignored if it was previously called)
        root_logger.setLevel(numeric_level)
        
        logger.info(f"Logging configured at {log_level} level")
        
    def load_config(self, args: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Load configuration from command-line arguments and/or config file.
        
        Args:
            args: Command-line arguments (uses sys.argv if None)
            
        Returns:
            Dictionary of configuration settings
            
        Raises:
            ConfigurationError: If configuration is invalid
        """
        try:
            self.config.load(args)
            config_dict = self.config.get_all()
            logger.info(f"Configuration loaded successfully")
            logger.debug(f"Configuration: {config_dict}")
            return config_dict
        except ConfigurationError as e:
            logger.error(f"Configuration error: {e}")
            raise
            
    def initialize_components(self) -> None:
        """
        Initialize all components needed for the annotation pipeline.
        
        Sets up parsers, analyzers, and processors based on the loaded configuration.
        """
        logger.info("Initializing pipeline components")
        
        # Create path analyzer
        self.path_analyzer = None  # Will be instantiated when modules are imported
        
        # Create feature graph
        self.feature_graph = None  # Will be instantiated when modules are imported
        
        # Create alignment processor with parallelization settings
        self.alignment_processor = None  # Will be instantiated when modules are imported
        
        # Create impact classifier
        self.impact_classifier = None  # Will be instantiated when modules are imported
        
        # Create summary generator
        self.summary_generator = None  # Will be instantiated when modules are imported
        
        # Create report generator
        self.report_generator = None  # Will be instantiated when modules are imported
        
        logger.info("Components initialized successfully")
        
    def load_input_files(self) -> Tuple[bool, str]:
        """
        Load all input files specified in the configuration.
        
        Returns:
            Tuple of (success, error_message)
        """
        logger.info("Loading input files")
        
        try:
            # Get file paths from config
            resources = self.config.get_resource_files()
            gfa_file = resources.get('gfa_file')
            gff_file = resources.get('gff3_file')
            fasta_file = resources.get('reference_fasta')
            
            # Verify files exist
            for file_path, file_desc in [
                (gfa_file, "GFA"),
                (gff_file, "GFF3"),
                (fasta_file, "Reference FASTA")
            ]:
                if not os.path.exists(file_path):
                    return False, f"{file_desc} file not found: {file_path}"
            
            # Save file paths in intermediate data
            self.intermediate_data['input_files'] = {
                'gfa_file': gfa_file,
                'gff_file': gff_file,
                'fasta_file': fasta_file
            }
            
            logger.info(f"Input files validated successfully")
            return True, ""
            
        except Exception as e:
            error_msg = f"Error loading input files: {str(e)}"
            logger.error(error_msg, exc_info=True)
            return False, error_msg
            
    def select_paths(self) -> List[str]:
        """
        Select paths to annotate based on configuration settings.
        
        Returns:
            List of selected path IDs
        """
        sample_names = self.config.get('sample_names')
        haplotype_ids = self.config.get('haplotype_ids')
        path_ids = self.config.get('path_ids')
        
        # Mock implementation until path_analyzer is implemented
        selected_paths = []
        
        # Handle explicit path IDs
        if isinstance(path_ids, str) and path_ids.strip():
            selected_paths = [p.strip() for p in path_ids.split(',') if p.strip()]
        elif isinstance(path_ids, list):
            selected_paths = path_ids
        
        # In a real implementation, we would also handle sample_names and haplotype_ids here
        # For now, just log if they were provided
        if sample_names:
            logger.info(f"Sample filtering requested but not yet implemented: {sample_names}")
        if haplotype_ids:
            logger.info(f"Haplotype filtering requested but not yet implemented: {haplotype_ids}")
        
        # If path_ids, sample_names, and haplotype_ids are all None/empty,
        # we would typically use a default set of paths or all available paths
        if not selected_paths and not sample_names and not haplotype_ids:
            logger.info("No specific paths requested, using default path selection strategy")
            # In a real implementation, this would select some default paths
            # For now, just return an empty list
        
        logger.info(f"Selected {len(selected_paths)} paths for annotation")
        if selected_paths:
            logger.debug(f"Selected paths: {selected_paths}")
        
        return selected_paths
        
    def save_intermediate_data(self, output_dir: str) -> None:
        """
        Save intermediate data for later visualization.
        
        Args:
            output_dir: Directory to save the intermediate data
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        # Save input file paths
        input_files = os.path.join(output_dir, 'input_files.json')
        with open(input_files, 'w') as f:
            json.dump(self.intermediate_data.get('input_files', {}), f, indent=2)
        
        # Save path information
        paths_file = os.path.join(output_dir, 'selected_paths.json')
        with open(paths_file, 'w') as f:
            json.dump(self.intermediate_data.get('selected_paths', []), f, indent=2)
        
        # Save configuration
        config_file = os.path.join(output_dir, 'configuration.json')
        with open(config_file, 'w') as f:
            # Filter out sensitive or unprintable data
            config_data = {k: v for k, v in self.config.get_all().items()
                          if isinstance(v, (str, int, float, bool, list, dict))}
            json.dump(config_data, f, indent=2)
        
        logger.info(f"Saved intermediate data to {output_dir}")
        
    def run(self, args: Optional[List[str]] = None) -> int:
        """
        Run the full annotation pipeline.
        
        Args:
            args: Command-line arguments (uses sys.argv if None)
            
        Returns:
            Exit code (0 for success, non-zero for error)
        """
        start_time = time.time()
        
        try:
            # Load configuration
            config_dict = self.load_config(args)
            
            # Configure logging
            self.configure_logging(config_dict.get('log_level', 'INFO'))
            
            # Initialize components
            self.initialize_components()
            
            # Load input files
            success, error_msg = self.load_input_files()
            if not success:
                logger.error(error_msg)
                return 1
                
            # Select paths to annotate
            selected_paths = self.select_paths()
            self.intermediate_data['selected_paths'] = selected_paths
            
            # Only warn and exit if path_ids were explicitly provided but no paths were found
            path_ids = self.config.get('path_ids')
            if not selected_paths and path_ids and isinstance(path_ids, str) and path_ids.strip():
                logger.warning(f"No paths found matching the provided path_ids: {path_ids}")
                logger.warning("Exiting without processing.")
                return 0
                
            # Save intermediate data if requested
            if config_dict.get('save_intermediate'):
                intermediate_dir = config_dict.get('intermediate_dir', 'intermediate_data')
                self.save_intermediate_data(intermediate_dir)
                
            # TODO: Implement the actual feature alignment, impact analysis, and report generation
            # These steps will be added in future implementations
            
            elapsed_time = time.time() - start_time
            logger.info(f"Pipeline completed successfully in {elapsed_time:.2f} seconds")
            return 0
            
        except ConfigurationError as e:
            logger.error(f"Configuration error: {e}")
            return 1
        except Exception as e:
            logger.error(f"Unexpected error: {e}", exc_info=True)
            return 1
            
def main() -> int:
    """
    Command-line entry point.
    
    Returns:
        Exit code (0 for success, non-zero for error)
    """
    tool = HaplotypeAnnotationTool()
    return tool.run()
    
if __name__ == "__main__":
    sys.exit(main())
