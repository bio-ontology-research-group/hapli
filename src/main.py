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
from typing import Dict, List, Any, Optional, Tuple, Set

# Configuration and Core
from src.config import Config, ConfigurationError

# Parsers
from src.parsers.gfa_parser import GFAParser
from src.parsers.gff_parser import GFF3Parser
from src.parsers.fasta_parser import FastaParser
from src.parsers.feature_graph import FeatureGraph
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

# Path Analysis
from src.path_analysis import PathAnalyzer

# Parallelization (assuming needed by components)
from src.parallel.task_manager import create_worker_pool, execute_parallel
from src.parallel.hierarchical_executor import HierarchicalExecutor
# from src.parallel.parallel_alignment import ParallelAligner # Potentially used within AlignmentProcessor

# Alignment
# NOTE: Assuming AlignmentProcessor exists and handles alignment logic
# from src.alignment.alignment_processor import AlignmentProcessor
# from src.alignment.minimap_wrapper import MinimapAligner # Likely used by AlignmentProcessor

# Analysis
# NOTE: Assuming these classes exist and perform their described functions
from src.analysis.impact_classifier import ImpactClassifier, ImpactType
from src.analysis.variant_detector import VariantDetector
from src.analysis.reconciliation import FeatureReconciler
from src.analysis.summary_generator import SummaryGenerator, AnalysisSummary, FeatureSummary

# Reporting
# NOTE: Assuming these classes exist and perform their described functions
from src.reporting.report_generator import ReportGenerator
from src.reporting.comparative_report import ComparativeReportGenerator
from src.reporting.formatters.tsv_formatter import TSVFormatter
from src.reporting.formatters.json_formatter import JSONFormatter
from src.reporting.formatters.rdf_formatter import RDFFormatter


# Configure logger for this module
logger = logging.getLogger(__name__)

# Placeholder classes for components not fully defined in context
# These allow the main structure to be built. Replace with actual imports later.
class PlaceholderAlignmentProcessor:
    def __init__(self, reference_fasta: Dict[str, SeqRecord], feature_graph: FeatureGraph, num_workers: int, batch_size: int):
        self.reference_fasta = reference_fasta
        self.feature_graph = feature_graph
        self.num_workers = num_workers
        self.batch_size = batch_size
        logger.info(f"PlaceholderAlignmentProcessor initialized (workers={num_workers}, batch={batch_size})")

    def align_features_to_path(self, path_id: str, path_sequence: str, features_to_align: List[SeqFeature]) -> Dict[str, Tuple[Optional[SeqFeature], Optional[float], Optional[float], Optional[float]]]:
        """ Placeholder: Returns mock alignment results """
        logger.info(f"Placeholder: Aligning {len(features_to_align)} features to path {path_id}...")
        # Simulate alignment: return some features as aligned, some as not
        results = {}
        for i, feature in enumerate(features_to_align):
            feature_id = feature.id
            if i % 3 != 0: # Simulate some features aligning
                 # Mock aligned feature, score, coverage, identity
                aligned_feature = feature # In reality, location/sequence might change
                results[feature_id] = (aligned_feature, 0.95, 0.98, 0.99)
            else: # Simulate alignment failure
                results[feature_id] = (None, None, None, None)
        time.sleep(0.01) # Simulate work slightly faster
        logger.info(f"Placeholder: Finished aligning features to path {path_id}")
        return results

class PlaceholderImpactClassifier:
    def __init__(self, identity_threshold: float = 0.9, coverage_threshold: float = 0.8):
         self.identity_threshold = identity_threshold
         self.coverage_threshold = coverage_threshold
         logger.info(f"PlaceholderImpactClassifier initialized (id_thresh={identity_threshold}, cov_thresh={coverage_threshold})")

    def classify_feature_impact(self, reference_feature: SeqFeature, aligned_feature: Optional[SeqFeature], alignment_score: Optional[float], coverage: Optional[float], identity: Optional[float]) -> Tuple[ImpactType, Dict]:
        """ Placeholder: Classifies based on simple presence/absence """
        if aligned_feature and identity is not None and coverage is not None:
            if identity >= self.identity_threshold and coverage >= self.coverage_threshold:
                return ImpactType.PRESENT, {"identity": identity, "coverage": coverage}
            else:
                 return ImpactType.MODIFIED, {"identity": identity, "coverage": coverage}
        else:
            return ImpactType.ABSENT, {}

class PlaceholderVariantDetector:
    def __init__(self):
        self.reference_fasta = None # Allow setting reference later
        logger.info("PlaceholderVariantDetector initialized")

    def detect_variants(self, reference_feature: SeqFeature, aligned_feature: SeqFeature) -> List[Any]:
        """ Placeholder: Returns mock variants """
        # Simulate variant detection if features differ slightly (e.g., modified)
        # In a real scenario, compare sequences or use CIGAR
        if not self.reference_fasta:
             logger.warning("Reference FASTA not set for PlaceholderVariantDetector. Cannot compare sequences.")
             return []
        try:
            # Ensure feature.ref points to a valid key in reference_sequences
            if reference_feature.ref not in self.reference_fasta:
                 logger.warning(f"Reference sequence ID '{reference_feature.ref}' not found for feature {reference_feature.id}")
                 return []
            # Assuming aligned_feature also has a .ref pointing to the *path sequence* (which isn't stored here)
            # This placeholder logic is flawed as it needs the actual aligned sequence context.
            # Let's just mock based on ID for simplicity in placeholder.
            if "modified" in reference_feature.id.lower(): # Mock based on ID
                 # Mock variant data
                 return [{"type": "SNP", "pos": 10, "ref": "A", "alt": "T"}]
        except Exception as e:
             logger.warning(f"Error extracting sequence for variant detection placeholder: {e}")
        return []

class PlaceholderFeatureReconciler:
    def __init__(self):
        logger.info("PlaceholderFeatureReconciler initialized")

    def reconcile_feature_hierarchy(self, feature_graph: Dict[str, List[str]], aligned_features: Dict[str, SeqFeature]) -> Dict[str, List[Any]]:
        """ Placeholder: Returns empty reconciliation results """
        logger.info("Placeholder: Skipping feature reconciliation")
        return {} # No conflicts simulated

class PlaceholderSummaryGenerator:
    def __init__(self):
        logger.info("PlaceholderSummaryGenerator initialized")

    def generate_summary(self, path_id: str, analysis_results: Dict[str, Dict]) -> AnalysisSummary:
        """ Placeholder: Creates a basic summary """
        feature_summaries = []
        for feature_id, results in analysis_results.items():
            # Ensure FeatureSummary is instantiated correctly based on its definition
            # Assuming FeatureSummary takes these kwargs
            fs = FeatureSummary(
                feature_id=feature_id,
                impact=results.get('impact_type', ImpactType.UNKNOWN),
                metrics=results.get('impact_details', {}),
                variants=results.get('variants', []),
                reconciliation_status=results.get('reconciliation', 'N/A')
            )
            feature_summaries.append(fs)

        # Assuming AnalysisSummary takes path_id and features list
        return AnalysisSummary(path_id=path_id, features=feature_summaries)

class PlaceholderReportGenerator:
    def __init__(self, formatter):
        self.formatter = formatter
        logger.info(f"PlaceholderReportGenerator initialized with {type(formatter).__name__}")

    def generate(self, summary: AnalysisSummary, output_file: Optional[str]):
        """ Placeholder: Formats and prints/writes summary """
        try:
            # Assume formatter.format exists and works
            formatted_output = self.formatter.format(summary, output_file or "stdout")
            if output_file:
                # Ensure directory exists before writing
                output_dir = os.path.dirname(output_file)
                if output_dir and not os.path.exists(output_dir):
                    os.makedirs(output_dir, exist_ok=True)
                with open(output_file, 'w') as f:
                    f.write(formatted_output)
                logger.info(f"Placeholder: Report written to {output_file}")
            else:
                # Use print for stdout simulation
                print(formatted_output)
                logger.info("Placeholder: Report written to stdout")
        except Exception as e:
             logger.error(f"PlaceholderReportGenerator failed: {e}", exc_info=True)
             raise # Re-raise to indicate failure

class PlaceholderComparativeReportGenerator:
     def __init__(self, formatter):
         self.formatter = formatter
         logger.info(f"PlaceholderComparativeReportGenerator initialized with {type(formatter).__name__}")

     def generate(self, summaries: Dict[str, AnalysisSummary], output_file: Optional[str]):
         """ Placeholder: Formats and prints/writes comparative summary """
         try:
             # Mock comparative data structure
             # Check if AnalysisSummary has a to_dict method, otherwise use a placeholder
             comparative_data = {}
             for path_id, summary in summaries.items():
                 if hasattr(summary, 'to_dict'):
                     comparative_data[path_id] = summary.to_dict()
                 else:
                     # Placeholder if to_dict is missing
                     comparative_data[path_id] = {"path_id": summary.path_id, "feature_count": len(summary.features)}

             # Assume formatter.format_comparative exists
             formatted_output = self.formatter.format_comparative(comparative_data, output_file or "stdout")
             if output_file:
                 # Ensure directory exists before writing
                 output_dir = os.path.dirname(output_file)
                 if output_dir and not os.path.exists(output_dir):
                     os.makedirs(output_dir, exist_ok=True)
                 with open(output_file, 'w') as f:
                     f.write(formatted_output)
                 logger.info(f"Placeholder: Comparative report written to {output_file}")
             else:
                 # Use print for stdout simulation
                 print(formatted_output)
                 logger.info("Placeholder: Comparative report written to stdout")
         except Exception as e:
             logger.error(f"PlaceholderComparativeReportGenerator failed: {e}", exc_info=True)
             raise # Re-raise to indicate failure


# --- Main Tool Class ---

class HaplotypeAnnotationTool:
    """
    Main class for the Haplotype Annotation Tool.

    Orchestrates the entire annotation pipeline.
    """

    def __init__(self):
        """Initialize the annotation tool."""
        self.config = Config()
        self.gfa_parser = GFAParser()
        self.gff_parser = GFF3Parser()
        self.fasta_parser = FastaParser()
        self.feature_graph_builder = FeatureGraph()
        self.path_analyzer = PathAnalyzer() # PathAnalyzer now holds GFA data internally

        # Parsed data storage
        self.gfa_data = None
        self.features: Dict[str, SeqFeature] = {}
        self.reference_sequences: Dict[str, SeqRecord] = {}
        self.feature_graph: Optional[Any] = None # Should be nx.DiGraph

        # Components (initialized later)
        self.alignment_processor = None
        self.impact_classifier = None
        self.variant_detector = None
        self.feature_reconciler = None
        self.summary_generator = None
        self.report_generator = None
        self.comparative_report_generator = None

        # Results storage
        self.analysis_summaries: Dict[str, AnalysisSummary] = {}
        self.intermediate_data: Dict[str, Any] = {}


    def configure_logging(self, log_level: str) -> None:
        """Configure logging using explicit handlers."""
        numeric_level = getattr(logging, log_level.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f"Invalid log level: {log_level}")

        root_logger = logging.getLogger()
        root_logger.setLevel(numeric_level) # Set level on root logger

        # Remove existing handlers attached to the root logger
        # This prevents duplicate messages if called multiple times or in tests
        for handler in root_logger.handlers[:]:
            root_logger.removeHandler(handler)

        # Create a default StreamHandler (writes to stderr)
        handler = logging.StreamHandler(sys.stderr)
        handler.setLevel(numeric_level) # Set level on handler

        # Create formatter and add it to the handler
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)

        # Add the handler to the root logger
        root_logger.addHandler(handler)

        # Log the configuration action using the module-level logger
        logger.info(f"Logging configured at {log_level} level")


    def load_config(self, args: Optional[List[str]] = None) -> Dict[str, Any]:
        """Load configuration."""
        try:
            self.config.load(args)
            config_dict = self.config.get_all()
            logger.info("Configuration loaded successfully")
            logger.debug(f"Configuration: {json.dumps(config_dict, indent=2, default=str)}") # Use json for cleaner debug output
            return config_dict
        except ConfigurationError as e:
            logger.error(f"Configuration error: {e}")
            raise

    def initialize_components(self) -> None:
        """Initialize pipeline components based on configuration."""
        logger.info("Initializing pipeline components")
        config_dict = self.config.get_all()

        # --- Instantiate Real Components (where possible) ---
        # Note: Pass necessary data like reference sequences or feature graph if needed at init
        # For now, pass placeholders or None where data isn't loaded yet

        # Alignment Processor (using placeholder)
        num_workers = config_dict.get('num_workers', os.cpu_count() or 4)
        batch_size = config_dict.get('batch_size', 100)
        # Real AlignmentProcessor would need reference sequences, feature graph etc.
        self.alignment_processor = PlaceholderAlignmentProcessor(
             reference_fasta=self.reference_sequences, # Pass loaded sequences
             feature_graph=self.feature_graph, # Pass loaded graph
             num_workers=num_workers,
             batch_size=batch_size
        )

        # Analysis Components (using placeholders)
        # TODO: Pass actual config thresholds when available
        self.impact_classifier = PlaceholderImpactClassifier()
        self.variant_detector = PlaceholderVariantDetector()
        # Pass reference sequences to variant detector if needed for extraction
        self.variant_detector.reference_fasta = self.reference_sequences
        self.feature_reconciler = PlaceholderFeatureReconciler()
        self.summary_generator = PlaceholderSummaryGenerator()

        # Reporting Components (using placeholders)
        output_format = config_dict.get('output_format', 'tsv')
        rdf_format = config_dict.get('rdf_format', 'turtle') # Needed if output_format is rdf

        if output_format == 'tsv':
            formatter = TSVFormatter()
        elif output_format == 'json':
            formatter = JSONFormatter()
        elif output_format == 'rdf':
            formatter = RDFFormatter() # Real formatter might need rdf_format
        else:
            raise ConfigurationError(f"Unsupported output format: {output_format}")

        self.report_generator = PlaceholderReportGenerator(formatter)
        self.comparative_report_generator = PlaceholderComparativeReportGenerator(formatter)

        logger.info("Components initialized successfully")


    def load_input_files(self) -> Tuple[bool, str]:
        """Load and parse GFA, GFF3, and FASTA files."""
        logger.info("Loading and parsing input files...")
        try:
            resources = self.config.get_resource_files()
            gfa_file = resources.get('gfa_file')
            gff_file = resources.get('gff3_file')
            fasta_file = resources.get('reference_fasta')

            if not all([gfa_file, gff_file, fasta_file]):
                 return False, "Missing one or more required input file paths in configuration."

            # Parse GFA
            logger.info(f"Parsing GFA file: {gfa_file}")
            self.gfa_data = self.gfa_parser.parse(gfa_file)
            logger.info(f"GFA parsing complete. Found {len(self.gfa_parser.get_segments())} segments, {len(self.gfa_parser.get_paths())} paths.")
            # Load GFA data into PathAnalyzer
            self.path_analyzer.load_gfa(self.gfa_data)


            # Parse FASTA
            logger.info(f"Parsing reference FASTA file: {fasta_file}")
            self.reference_sequences = self.fasta_parser.parse(fasta_file)
            logger.info(f"FASTA parsing complete. Found {len(self.reference_sequences)} sequences.")

            # Parse GFF3
            logger.info(f"Parsing GFF3 file: {gff_file}")
            # GFFParser returns list of SeqRecord, each containing features for a contig/chromosome
            gff_records = self.gff_parser.parse(gff_file)
            # Index all features by ID
            self.features = {}
            feature_count = 0
            for record in gff_records:
                 # Check if reference sequence exists for this record
                 if record.id not in self.reference_sequences:
                     logger.warning(f"Reference sequence '{record.id}' not found in FASTA for features in GFF. Skipping features on this sequence.")
                     continue

                 for feature in self.gff_parser._iter_features(record.features): # Use internal iterator helper
                     if feature.id:
                         # Store reference SeqRecord ID with the feature for context
                         feature.ref = record.id
                         self.features[feature.id] = feature
                         feature_count += 1
                     else:
                         logger.warning(f"Feature of type '{feature.type}' lacks an ID. Skipping.")
            logger.info(f"GFF3 parsing complete. Indexed {feature_count} features with IDs.")

            # Build Feature Graph
            logger.info("Building feature relationship graph...")
            self.feature_graph = self.feature_graph_builder.build_from_features(self.features)
            logger.info(f"Feature graph built with {self.feature_graph.number_of_nodes()} nodes and {self.feature_graph.number_of_edges()} edges.")


            # Update intermediate data
            self.intermediate_data['input_files'] = resources
            self.intermediate_data['parsed_features_count'] = len(self.features)
            self.intermediate_data['reference_seq_ids'] = list(self.reference_sequences.keys())

            # Pass loaded data to components that need it now (if not passed at init)
            if hasattr(self.variant_detector, 'reference_fasta'):
                 self.variant_detector.reference_fasta = self.reference_sequences
            # Update AlignmentProcessor if it needs data post-init
            if hasattr(self.alignment_processor, 'reference_fasta'):
                 self.alignment_processor.reference_fasta = self.reference_sequences
            if hasattr(self.alignment_processor, 'feature_graph'):
                 self.alignment_processor.feature_graph = self.feature_graph


            logger.info("Input files loaded and parsed successfully.")
            return True, ""

        except Exception as e:
            error_msg = f"Error loading or parsing input files: {str(e)}"
            logger.error(error_msg, exc_info=True)
            return False, error_msg


    def select_paths(self) -> List[str]:
        """Select paths using PathAnalyzer based on configuration."""
        logger.info("Selecting paths for annotation...")
        sample_names_str = self.config.get('sample_names')
        haplotype_ids_str = self.config.get('haplotype_ids')
        path_ids_str = self.config.get('path_ids')

        sample_names = [s.strip() for s in sample_names_str.split(',')] if sample_names_str else None
        haplotype_ids = [h.strip() for h in haplotype_ids_str.split(',')] if haplotype_ids_str else None
        path_ids = [p.strip() for p in path_ids_str.split(',')] if path_ids_str else None

        selected_paths = self.path_analyzer.select_paths(
            sample_names=sample_names,
            haplotype_ids=haplotype_ids,
            path_ids=path_ids
        )

        if not selected_paths:
             # If specific filters were given but nothing matched, warn.
             if sample_names or haplotype_ids or path_ids:
                 logger.warning("No paths matched the specified selection criteria.")
             # If no filters were given, PathAnalyzer.select_paths should return all paths.
             # If it still returns empty, the GFA might have no paths.
             elif not self.gfa_parser.get_paths():
                  logger.warning("No paths found in the GFA file.")
             else:
                  # This case shouldn't happen if select_paths returns all paths by default
                  logger.warning("Path selection resulted in an empty list unexpectedly.")

        logger.info(f"Selected {len(selected_paths)} paths for annotation.")
        if selected_paths:
            logger.debug(f"Selected paths: {selected_paths}")

        return selected_paths


    def run_annotation_pipeline(self, selected_paths: List[str]) -> None:
        """Run the core annotation pipeline for selected paths."""
        if not selected_paths:
            logger.warning("No paths selected for annotation. Skipping pipeline.")
            return

        logger.info(f"Starting annotation pipeline for {len(selected_paths)} paths...")
        self.analysis_summaries = {} # Reset summaries

        # Get all top-level features (e.g., genes) to align first
        # Assumes FeatureGraph provides a way to get roots or top-level nodes
        top_level_feature_ids = [
             f_id for f_id, data in self.feature_graph.nodes(data=True)
             if self.feature_graph.in_degree(f_id) == 0 and f_id in self.features
        ]
        top_level_features = [self.features[f_id] for f_id in top_level_feature_ids]
        logger.info(f"Identified {len(top_level_features)} top-level features for initial alignment.")


        for path_id in selected_paths:
            logger.info(f"--- Processing path: {path_id} ---")
            path_start_time = time.time()

            try:
                # 1. Extract Path Sequence
                logger.debug(f"Extracting sequence for path {path_id}")
                path_sequence = self.path_analyzer.get_path_sequence(path_id)
                if not path_sequence:
                    logger.warning(f"Could not extract sequence for path {path_id}. Skipping.")
                    continue
                logger.info(f"Extracted path sequence (length: {len(path_sequence)})")

                # 2. Align Features to Path
                # Align top-level features first. AlignmentProcessor should handle hierarchy.
                logger.info(f"Aligning {len(top_level_features)} top-level features to path {path_id}")
                # AlignmentProcessor returns dict: feature_id -> (aligned_feature, score, coverage, identity)
                # It should internally handle aligning children within aligned parent boundaries.
                alignment_results = self.alignment_processor.align_features_to_path(
                    path_id=path_id,
                    path_sequence=path_sequence,
                    features_to_align=top_level_features # Pass only top-level, processor handles children
                )
                logger.info(f"Alignment complete for path {path_id}. Got results for {len(alignment_results)} features.")

                # 3. Analyze Alignments (Impact, Variants, Reconciliation)
                path_analysis_results = {}
                aligned_features_for_reconciliation = {} # Collect successfully aligned features

                # Iterate through all features (not just top-level) to classify impact
                for feature_id, reference_feature in self.features.items():
                    logger.debug(f"Analyzing feature: {feature_id}")
                    analysis_result = {}

                    # Get alignment result for this specific feature
                    # The AlignmentProcessor *should* return results for all relevant features (parents and children)
                    # If it only returns top-level, this needs adjustment. Assuming it returns all.
                    aln_data = alignment_results.get(feature_id)

                    if aln_data:
                        aligned_feature, score, coverage, identity = aln_data
                        analysis_result['alignment_score'] = score
                        analysis_result['coverage'] = coverage
                        analysis_result['identity'] = identity

                        if aligned_feature:
                             aligned_features_for_reconciliation[feature_id] = aligned_feature
                             # Detect Variants
                             analysis_result['variants'] = self.variant_detector.detect_variants(
                                 reference_feature=reference_feature,
                                 aligned_feature=aligned_feature
                             )
                        else: # Alignment attempted but failed quality checks?
                             analysis_result['variants'] = []

                        # Classify Impact
                        impact_type, impact_details = self.impact_classifier.classify_feature_impact(
                            reference_feature=reference_feature,
                            aligned_feature=aligned_feature, # Can be None if alignment failed
                            alignment_score=score,
                            coverage=coverage,
                            identity=identity
                        )
                        analysis_result['impact_type'] = impact_type
                        analysis_result['impact_details'] = impact_details

                    else: # Feature was not aligned (or not returned by aligner)
                        logger.debug(f"Feature {feature_id} not found in alignment results. Classifying as ABSENT.")
                        analysis_result['impact_type'] = ImpactType.ABSENT
                        analysis_result['impact_details'] = {}
                        analysis_result['variants'] = []


                    path_analysis_results[feature_id] = analysis_result

                # 4. Reconcile Feature Hierarchy (if necessary)
                # This step might modify aligned_features or add reconciliation info to path_analysis_results
                logger.info(f"Reconciling feature hierarchy for path {path_id}")
                # Build the parent->child map expected by the reconciler
                feature_hierarchy_graph = {
                    parent: list(self.feature_graph.successors(parent))
                    for parent in self.feature_graph.nodes() if self.feature_graph.out_degree(parent) > 0
                }
                reconciliation_results = self.feature_reconciler.reconcile_feature_hierarchy(
                     feature_graph=feature_hierarchy_graph,
                     aligned_features=aligned_features_for_reconciliation # Pass only successfully aligned ones
                )
                # Integrate reconciliation results back into path_analysis_results
                for feature_id, recon_list in reconciliation_results.items():
                     if feature_id in path_analysis_results:
                         # Store simplified status or details
                         path_analysis_results[feature_id]['reconciliation'] = str(recon_list) # Placeholder representation


                # 5. Generate Summary for Path
                logger.info(f"Generating summary for path {path_id}")
                path_summary = self.summary_generator.generate_summary(
                    path_id=path_id,
                    analysis_results=path_analysis_results
                )
                self.analysis_summaries[path_id] = path_summary

                path_elapsed = time.time() - path_start_time
                logger.info(f"--- Finished processing path {path_id} in {path_elapsed:.2f} seconds ---")

            except Exception as e:
                logger.error(f"Error processing path {path_id}: {e}", exc_info=True)
                # Optionally store error state for this path
                # self.analysis_summaries[path_id] = AnalysisSummary(path_id=path_id, error=str(e))


    def generate_report(self) -> None:
        """Generate the final report based on analysis summaries."""
        logger.info("Generating final report...")
        config_dict = self.config.get_all()
        output_file = config_dict.get('output_file')
        is_comparative = config_dict.get('comparative', False)

        if not self.analysis_summaries:
            logger.warning("No analysis summaries generated. Skipping report generation.")
            return

        try:
            if is_comparative:
                if len(self.analysis_summaries) < 2:
                    logger.warning("Comparative report requested, but fewer than 2 paths were analyzed. Generating standard report instead.")
                    # Fallback to standard report - just use the first summary if available
                    if self.analysis_summaries:
                         first_summary = next(iter(self.analysis_summaries.values()))
                         self.report_generator.generate(first_summary, output_file)
                    else:
                         logger.error("No summaries available to generate report.")

                else:
                    logger.info(f"Generating comparative report for {len(self.analysis_summaries)} paths.")
                    self.comparative_report_generator.generate(self.analysis_summaries, output_file)
            else:
                 # Handle single or multiple summaries for standard report
                 # Current placeholder just processes the first one.
                 # A real implementation might concatenate reports or require a single path.
                 if len(self.analysis_summaries) > 1:
                      logger.warning(f"Multiple paths analyzed ({len(self.analysis_summaries)}) but non-comparative report requested. Generating report for the first path only: {next(iter(self.analysis_summaries.keys()))}")

                 # Ensure there's at least one summary before accessing
                 if not self.analysis_summaries:
                      logger.error("No summaries available to generate report.")
                      return

                 first_summary = next(iter(self.analysis_summaries.values()))
                 logger.info(f"Generating standard report for path: {first_summary.path_id}")
                 self.report_generator.generate(first_summary, output_file)


            logger.info("Report generation complete.")

        except Exception as e:
            logger.error(f"Error generating report: {e}", exc_info=True)


    def save_intermediate_data(self, output_dir: str) -> None:
        """Save intermediate data."""
        try:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                logger.info(f"Created intermediate data directory: {output_dir}")

            # Save configuration
            config_file = os.path.join(output_dir, 'configuration.json')
            with open(config_file, 'w') as f:
                config_data = {k: v for k, v in self.config.get_all().items()
                              if isinstance(v, (str, int, float, bool, list, dict, type(None)))} # Allow None
                json.dump(config_data, f, indent=2, default=str) # Use default=str for non-serializable types

            # Save input file paths
            input_files_path = os.path.join(output_dir, 'input_files.json')
            with open(input_files_path, 'w') as f:
                json.dump(self.intermediate_data.get('input_files', {}), f, indent=2)

            # Save selected paths
            paths_file = os.path.join(output_dir, 'selected_paths.json')
            with open(paths_file, 'w') as f:
                json.dump(self.intermediate_data.get('selected_paths', []), f, indent=2)

            # Save basic feature info (IDs and types)
            features_info_file = os.path.join(output_dir, 'parsed_features_info.json')
            features_info = {fid: f.type for fid, f in self.features.items()}
            with open(features_info_file, 'w') as f:
                 json.dump(features_info, f, indent=2)

            # Optionally save summaries if generated
            summaries_file = os.path.join(output_dir, 'analysis_summaries.json')
            try:
                # Convert summaries to dicts if they have a method for it
                summaries_dict = {}
                for path_id, summary in self.analysis_summaries.items():
                     if hasattr(summary, 'to_dict'):
                         summaries_dict[path_id] = summary.to_dict()
                     else: # Placeholder if no to_dict method
                          summaries_dict[path_id] = {"path_id": summary.path_id, "feature_count": len(summary.features)}

                with open(summaries_file, 'w') as f:
                    json.dump(summaries_dict, f, indent=2, default=str) # Use default=str
            except Exception as e:
                 logger.warning(f"Could not save analysis summaries to intermediate data: {e}")


            logger.info(f"Saved intermediate data to {output_dir}")
        except Exception as e:
             logger.error(f"Failed to save intermediate data to {output_dir}: {e}", exc_info=True)


    def run(self, args: Optional[List[str]] = None) -> int:
        """Run the full annotation pipeline."""
        start_time = time.time()
        exit_code = 0
        log_level_for_init = 'INFO' # Default level before config is loaded

        try:
            # Attempt to parse log level from args early for initial setup
            # This is a bit hacky but allows seeing config load errors if DEBUG is passed
            temp_parser = argparse.ArgumentParser(add_help=False)
            temp_parser.add_argument('--log-level', default='INFO')
            parsed_args, _ = temp_parser.parse_known_args(args)
            log_level_for_init = parsed_args.log_level

            # Configure logging minimally first
            self.configure_logging(log_level_for_init)

            # 1. Load Config
            config_dict = self.load_config(args)

            # 2. Re-Configure Logging based on final config value
            final_log_level = config_dict.get('log_level', 'INFO')
            if final_log_level.upper() != log_level_for_init.upper():
                 self.configure_logging(final_log_level)
            logger.info("--- Starting Haplotype Annotation Tool ---")


            # 3. Load Input Files (includes parsing)
            success, error_msg = self.load_input_files()
            if not success:
                # Error already logged by load_input_files
                return 1

            # 4. Initialize Components (after data is loaded)
            self.initialize_components()

            # 5. Select Paths
            selected_paths = self.select_paths()
            self.intermediate_data['selected_paths'] = selected_paths # Store selected paths

            # 6. Run Annotation Pipeline
            self.run_annotation_pipeline(selected_paths)

            # 7. Generate Report
            self.generate_report()

            # 8. Save Intermediate Data (if requested)
            if config_dict.get('save_intermediate'):
                intermediate_dir = config_dict.get('intermediate_dir', 'intermediate_data')
                self.save_intermediate_data(intermediate_dir)

            logger.info("--- Pipeline finished ---")

        except ConfigurationError as e:
            # Error should have been logged by load_config
            exit_code = 1
        except Exception as e:
            # Log the final unexpected error
            logger.critical(f"A critical unexpected error occurred during the pipeline execution: {e}", exc_info=True)
            exit_code = 1
        finally:
            elapsed_time = time.time() - start_time
            if exit_code == 0:
                logger.info(f"Haplotype Annotation Tool completed successfully in {elapsed_time:.2f} seconds.")
            else:
                logger.error(f"Haplotype Annotation Tool failed after {elapsed_time:.2f} seconds. Exit code: {exit_code}")

        return exit_code


def main() -> int:
    """Command-line entry point."""
    # Basic logging setup before full configuration, captures early errors
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    tool = HaplotypeAnnotationTool()
    return tool.run()

if __name__ == "__main__":
    # Ensure the script can find modules in the 'src' directory
    # This is often needed when running a script directly inside a package
    script_dir = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.dirname(script_dir) # Assumes src is one level up from main.py's dir
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)

    sys.exit(main())
