#!/usr/bin/env python3
"""
Align GFF3 features to genome graph paths with hierarchical constraints.

This script aligns each feature in a GFF3 file to all paths in a genome graph
(GFA/VG/XG format) using minimap2, respecting the hierarchical structure of
features and processing in parallel.
"""

import argparse
import logging
import sys
from pathlib import Path

# Add hapli to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from hapli.gff_alignment import GFFAligner, setup_logging


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Align GFF3 features to genome graph paths",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        "gff_file",
        type=Path,
        help="Path to GFF3 file containing features to align"
    )
    
    parser.add_argument(
        "reference_genome",
        type=Path,
        help="Path to reference genome FASTA file"
    )
    
    parser.add_argument(
        "graph_file",
        type=Path,
        help="Path to genome graph file (GFA/VG/XG format)"
    )
    
    parser.add_argument(
        "output_file",
        type=Path,
        help="Path to output alignment file"
    )
    
    # Optional arguments
    parser.add_argument(
        "--reference-path-name",
        type=str,
        help="Name of reference path in graph (for GFA files)"
    )
    
    parser.add_argument(
        "--output-format",
        choices=["json", "gam", "bam"],
        default="json",
        help="Output format for alignments"
    )
    
    parser.add_argument(
        "--max-workers",
        type=int,
        default=4,
        help="Number of parallel workers for alignment"
    )
    
    parser.add_argument(
        "--max-features",
        type=int,
        help="Maximum number of features to process (for testing)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Validate input files
    if not args.gff_file.exists():
        logging.error(f"GFF3 file not found: {args.gff_file}")
        sys.exit(1)
    
    if not args.reference_genome.exists():
        logging.error(f"Reference genome file not found: {args.reference_genome}")
        sys.exit(1)
    
    if not args.graph_file.exists():
        logging.error(f"Graph file not found: {args.graph_file}")
        sys.exit(1)
    
    # Validate output format and file extension
    if args.output_format == "gam" and not str(args.output_file).endswith('.gam'):
        logging.warning(f"GAM format specified but output file doesn't end with .gam: {args.output_file}")
    elif args.output_format == "bam" and not str(args.output_file).endswith('.bam'):
        logging.warning(f"BAM format specified but output file doesn't end with .bam: {args.output_file}")
    elif args.output_format == "json" and not str(args.output_file).endswith('.json'):
        logging.warning(f"JSON format specified but output file doesn't end with .json: {args.output_file}")
    
    # Create output directory if needed
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Initialize aligner
    logging.info("Initializing GFF aligner")
    aligner = GFFAligner(
        gff_path=args.gff_file,
        reference_path=args.reference_genome,
        graph_file=args.graph_file,
        reference_path_name=args.reference_path_name
    )
    
    try:
        # Process alignment
        logging.info("Starting graph alignment process")
        alignments = aligner.process_graph_alignment(
            output_path=args.output_file,
            output_format=args.output_format,
            max_workers=args.max_workers
        )
        
        # Print summary
        total_features = len(aligner.sorted_features)
        aligned_features = len(set(a.feature_id for a in alignments))
        total_paths = len(aligner.graph_paths)
        used_paths = len(set(a.path_name for a in alignments))
        
        print(f"\nAlignment Summary:")
        print(f"================")
        print(f"Total features: {total_features}")
        print(f"Features with alignments: {aligned_features}")
        print(f"Total graph paths: {total_paths}")
        print(f"Paths with alignments: {used_paths}")
        print(f"Total alignments: {len(alignments)}")
        print(f"Output format: {args.output_format.upper()}")
        print(f"Output file: {args.output_file}")
        
        # Feature type breakdown
        feature_types = {}
        for alignment in alignments:
            feature_types[alignment.feature_type] = feature_types.get(alignment.feature_type, 0) + 1
        
        if feature_types:
            print(f"\nAlignments by feature type:")
            for feature_type, count in sorted(feature_types.items()):
                print(f"  {feature_type}: {count}")
        
        # Hierarchy level breakdown
        hierarchy_levels = {}
        for alignment in alignments:
            level = alignment.hierarchy_level
            hierarchy_levels[level] = hierarchy_levels.get(level, 0) + 1
        
        if hierarchy_levels:
            print(f"\nAlignments by hierarchy level:")
            for level, count in sorted(hierarchy_levels.items()):
                print(f"  Level {level}: {count}")
        
        # Additional info for BAM output
        if args.output_format == "bam":
            output_dir = args.output_file.parent
            bam_files = list(output_dir.glob("*.bam"))
            print(f"\nBAM files created: {len(bam_files)}")
            for bam_file in sorted(bam_files):
                print(f"  {bam_file}")
        
        logging.info("Graph alignment completed successfully")
        
    except Exception as e:
        logging.error(f"Error during graph alignment: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
