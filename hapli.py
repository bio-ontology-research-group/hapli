#!/usr/bin/env python3
"""
Hapli: Main script for pangenome variant impact analysis.

This script provides a unified interface to run the complete Hapli analysis workflow:
1. GFF alignment to pangenome graph
2. GAM parsing and grouping
3. Impact detection
4. Diploid analysis

Usage:
    python hapli.py --gff features.gff3 --reference ref.fa --graph graph.gfa --output-dir results/
"""

import argparse
import logging
import sys
import os
from pathlib import Path
from typing import Optional

# Add the current directory to Python path to allow imports
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

# Import the individual modules using direct imports
import hapli.gff_alignment as gff_alignment_module
import hapli.gam_parser as gam_parser_module
import hapli.impact_detector as impact_detector_module
import hapli.diploid_analyzer as diploid_analyzer_module


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def run_gff_alignment(gff_file: Path, reference_file: Path, graph_file: Path, 
                     output_file: Path, reference_path_name: str = "reference",
                     max_workers: int = 4, verbose: bool = False) -> None:
    """Run GFF alignment step."""
    logging.info("Step 1: Running GFF alignment to pangenome graph")
    
    aligner = gff_alignment_module.GFFAligner(
        gff_file=gff_file,
        reference_file=reference_file,
        graph_file=graph_file,
        reference_path_name=reference_path_name
    )
    
    # Align features in parallel
    alignments = aligner.align_features_parallel(max_workers=max_workers)
    
    # Write GAM output
    aligner.write_alignments(alignments, output_file, format_type="gam")
    
    logging.info(f"GFF alignment complete. Output: {output_file}")


def run_gam_parsing(gam_file: Path, output_file: Path, verbose: bool = False) -> None:
    """Run GAM parsing step."""
    logging.info("Step 2: Parsing GAM file and grouping alignments")
    
    parser = gam_parser_module.GAMParser(gam_file)
    parser.load_alignments()
    
    # Group alignments by sample/haplotype
    grouped_data = parser.group_alignments_by_sample_haplotype()
    
    # Save grouped data
    parser.save_grouped_data(output_file, grouped_data)
    
    logging.info(f"GAM parsing complete. Output: {output_file}")


def run_impact_detection(gam_data_file: Path, gfa_file: Path, gff_file: Path,
                        output_file: Path, verbose: bool = False) -> None:
    """Run impact detection step."""
    logging.info("Step 3: Detecting variant impacts on genomic features")
    
    detector = impact_detector_module.ImpactDetector(gfa_file, gff_file)
    
    # Load GAM data
    import json
    with open(gam_data_file, 'r') as f:
        gam_data = json.load(f)
    
    # Analyze impacts
    impact_results = detector.analyze_impacts(gam_data)
    
    # Save results
    detector.save_results(impact_results, output_file)
    
    logging.info(f"Impact detection complete. Output: {output_file}")


def run_diploid_analysis(impact_data_file: Path, output_file: Path, 
                        verbose: bool = False) -> None:
    """Run diploid analysis step."""
    logging.info("Step 4: Running diploid analysis and generating reports")
    
    analyzer = diploid_analyzer_module.DiploidAnalyzer()
    
    # Load impact data
    import json
    with open(impact_data_file, 'r') as f:
        impact_results = json.load(f)
    
    # Analyze diploid impacts
    diploid_results = analyzer.analyze_diploid_impacts(impact_results)
    
    # Save results
    analyzer.save_results(diploid_results, output_file)
    
    logging.info(f"Diploid analysis complete. Output: {output_file}")


def main():
    """Main entry point for Hapli analysis workflow."""
    parser = argparse.ArgumentParser(
        description="Hapli: Pangenome variant impact analysis workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete workflow
  python hapli.py --gff features.gff3 --reference ref.fa --graph graph.gfa --output-dir results/
  
  # Run with custom settings
  python hapli.py --gff features.gff3 --reference ref.fa --graph graph.gfa \\
                  --output-dir results/ --reference-path-name reference \\
                  --max-workers 8 --verbose
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--gff",
        type=Path,
        required=True,
        help="GFF3 file with genomic features"
    )
    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="Reference genome FASTA file"
    )
    parser.add_argument(
        "--graph",
        type=Path,
        required=True,
        help="Pangenome graph file (GFA format)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for analysis results"
    )
    
    # Optional arguments
    parser.add_argument(
        "--reference-path-name",
        type=str,
        default="reference",
        help="Name of reference path in the graph (default: reference)"
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=4,
        help="Maximum number of parallel workers (default: 4)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    parser.add_argument(
        "--step",
        type=str,
        choices=["alignment", "parsing", "impact", "diploid", "all"],
        default="all",
        help="Run specific step only (default: all)"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Validate input files
    if not args.gff.exists():
        logging.error(f"GFF file not found: {args.gff}")
        sys.exit(1)
    if not args.reference.exists():
        logging.error(f"Reference file not found: {args.reference}")
        sys.exit(1)
    if not args.graph.exists():
        logging.error(f"Graph file not found: {args.graph}")
        sys.exit(1)
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define intermediate file paths
    gam_file = args.output_dir / "alignments.gam"
    grouped_file = args.output_dir / "grouped_alignments.json"
    impact_file = args.output_dir / "impact_results.json"
    diploid_file = args.output_dir / "diploid_results.json"
    
    try:
        if args.step in ["alignment", "all"]:
            run_gff_alignment(
                gff_file=args.gff,
                reference_file=args.reference,
                graph_file=args.graph,
                output_file=gam_file,
                reference_path_name=args.reference_path_name,
                max_workers=args.max_workers,
                verbose=args.verbose
            )
        
        if args.step in ["parsing", "all"]:
            if args.step == "parsing" and not gam_file.exists():
                logging.error(f"GAM file not found: {gam_file}. Run alignment step first.")
                sys.exit(1)
            run_gam_parsing(
                gam_file=gam_file,
                output_file=grouped_file,
                verbose=args.verbose
            )
        
        if args.step in ["impact", "all"]:
            if args.step == "impact" and not grouped_file.exists():
                logging.error(f"Grouped alignments file not found: {grouped_file}. Run parsing step first.")
                sys.exit(1)
            run_impact_detection(
                gam_data_file=grouped_file,
                gfa_file=args.graph,
                gff_file=args.gff,
                output_file=impact_file,
                verbose=args.verbose
            )
        
        if args.step in ["diploid", "all"]:
            if args.step == "diploid" and not impact_file.exists():
                logging.error(f"Impact results file not found: {impact_file}. Run impact step first.")
                sys.exit(1)
            run_diploid_analysis(
                impact_data_file=impact_file,
                output_file=diploid_file,
                verbose=args.verbose
            )
        
        logging.info("Hapli analysis workflow completed successfully!")
        logging.info(f"Final results available in: {args.output_dir}")
        
        if args.step == "all":
            logging.info("Output files:")
            logging.info(f"  - GAM alignments: {gam_file}")
            logging.info(f"  - Grouped alignments: {grouped_file}")
            logging.info(f"  - Impact results: {impact_file}")
            logging.info(f"  - Diploid analysis: {diploid_file}")
    
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        if args.verbose:
            logging.exception("Full traceback:")
        sys.exit(1)


if __name__ == "__main__":
    main()
