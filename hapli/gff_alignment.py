#!/usr/bin/env python3
"""
GFF3 topological sorting and sequence alignment functionality for hapli.

This module provides core functionality to:
1. Topologically sort GFF3 features by hierarchy and size
2. Extract corresponding sequences from reference genome
3. Align features to GFA paths using minimap2
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict, deque
import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


class GFFAligner:
    """Handle GFF3 topological sorting and sequence alignment."""
    
    def __init__(self, gff_path: Path, reference_path: Path):
        """Initialize with GFF3 and reference files."""
        self.gff_path = gff_path
        self.reference_path = reference_path
        self.db = None
        self.reference_genome = {}
        self.sorted_features = []
        
    def load_gff_database(self) -> gffutils.FeatureDB:
        """Load GFF3 file into gffutils database."""
        logging.info(f"Loading GFF3 file: {self.gff_path}")
        
        # Create temporary database
        db_path = str(self.gff_path) + '.db'
        
        try:
            # Try to load existing database
            if Path(db_path).exists():
                self.db = gffutils.FeatureDB(db_path)
                logging.info("Loaded existing GFF3 database")
            else:
                raise FileNotFoundError("Database file does not exist")
        except (gffutils.exceptions.FeatureDBError, ValueError, FileNotFoundError):
            # Create new database
            logging.info("Creating new GFF3 database")
            self.db = gffutils.create_db(
                str(self.gff_path),
                db_path,
                force=True,
                keep_order=True,
                merge_strategy='merge',
                sort_attribute_values=True
            )
        
        return self.db
    
    def load_reference_genome(self) -> Dict[str, SeqRecord]:
        """Load reference genome from FASTA file."""
        logging.info(f"Loading reference genome: {self.reference_path}")
        
        genome = {}
        with open(self.reference_path, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                genome[record.id] = record
        
        logging.info(f"Loaded {len(genome)} chromosomes")
        self.reference_genome = genome
        return genome
    
    def get_feature_hierarchy_level(self, feature: gffutils.Feature) -> int:
        """Determine hierarchy level of a feature (0=top level, higher=nested)."""
        level = 0
        current = feature
        
        # Count parent relationships
        while True:
            try:
                parents = list(self.db.parents(current))
                if not parents:
                    break
                current = parents[0]  # Take first parent
                level += 1
            except (gffutils.exceptions.FeatureNotFoundError, StopIteration):
                break
        
        return level
    
    def get_feature_size(self, feature: gffutils.Feature) -> int:
        """Get the size of a feature."""
        return feature.end - feature.start + 1
    
    def topological_sort_features(self) -> List[gffutils.Feature]:
        """
        Topologically sort features by hierarchy and size.
        
        Returns features sorted by:
        1. Hierarchy level (outermost first)
        2. Size (largest first within same level)
        3. Genomic position (for deterministic ordering)
        """
        logging.info("Performing topological sort of GFF3 features")
        
        if not self.db:
            raise ValueError("GFF3 database not loaded")
        
        # Get all features
        all_features = list(self.db.all_features())
        logging.info(f"Found {len(all_features)} total features")
        
        # Group features by hierarchy level
        features_by_level = defaultdict(list)
        
        for feature in all_features:
            level = self.get_feature_hierarchy_level(feature)
            features_by_level[level].append(feature)
        
        # Sort features within each level by size (largest first) and position
        sorted_features = []
        
        for level in sorted(features_by_level.keys()):
            level_features = features_by_level[level]
            
            # Sort by size (descending) then by chromosome and position
            level_features.sort(key=lambda f: (
                -self.get_feature_size(f),  # Negative for descending order
                f.seqid,
                f.start
            ))
            
            sorted_features.extend(level_features)
            logging.info(f"Level {level}: {len(level_features)} features")
        
        self.sorted_features = sorted_features
        logging.info(f"Topological sort complete: {len(sorted_features)} features")
        
        return sorted_features
    
    def extract_feature_sequence(self, feature: gffutils.Feature) -> str:
        """Extract the DNA sequence for a feature from the reference genome."""
        if not self.reference_genome:
            raise ValueError("Reference genome not loaded")
        
        chromosome = feature.seqid
        if chromosome not in self.reference_genome:
            raise ValueError(f"Chromosome {chromosome} not found in reference genome")
        
        # Extract sequence (GFF3 uses 1-based coordinates)
        start_idx = feature.start - 1  # Convert to 0-based
        end_idx = feature.end
        
        sequence = str(self.reference_genome[chromosome].seq[start_idx:end_idx])
        
        # Handle strand
        if feature.strand == '-':
            from Bio.Seq import Seq
            sequence = str(Seq(sequence).reverse_complement())
        
        return sequence
    
    def get_feature_sequences(self, features: List[gffutils.Feature] = None) -> List[Tuple[gffutils.Feature, str]]:
        """
        Extract sequences for features.
        
        Args:
            features: List of features to process. If None, uses sorted_features.
            
        Returns:
            List of (feature, sequence) tuples
        """
        if features is None:
            if not self.sorted_features:
                raise ValueError("No sorted features available. Run topological_sort_features() first.")
            features = self.sorted_features
        
        logging.info(f"Extracting sequences for {len(features)} features")
        
        feature_sequences = []
        
        for feature in features:
            try:
                sequence = self.extract_feature_sequence(feature)
                feature_sequences.append((feature, sequence))
            except Exception as e:
                logging.warning(f"Could not extract sequence for feature {feature.id}: {e}")
                continue
        
        logging.info(f"Successfully extracted {len(feature_sequences)} sequences")
        return feature_sequences
    
    def print_feature_sequences(self, max_features: int = None, max_seq_length: int = 100) -> None:
        """
        Print topologically sorted features and their sequences.
        
        Args:
            max_features: Maximum number of features to print (None for all)
            max_seq_length: Maximum sequence length to display (truncate longer sequences)
        """
        if not self.sorted_features:
            logging.error("No sorted features available. Run topological_sort_features() first.")
            return
        
        features_to_print = self.sorted_features
        if max_features:
            features_to_print = features_to_print[:max_features]
        
        feature_sequences = self.get_feature_sequences(features_to_print)
        
        print(f"\nTopologically sorted features ({len(feature_sequences)} total):")
        print("=" * 80)
        
        for i, (feature, sequence) in enumerate(feature_sequences, 1):
            # Get hierarchy level for display
            level = self.get_feature_hierarchy_level(feature)
            size = self.get_feature_size(feature)
            
            # Truncate sequence if too long
            display_seq = sequence
            if len(sequence) > max_seq_length:
                display_seq = sequence[:max_seq_length] + f"... ({len(sequence)} bp total)"
            
            # Get feature attributes for display
            feature_id = getattr(feature, 'id', 'N/A')
            feature_name = feature.attributes.get('Name', ['N/A'])[0] if 'Name' in feature.attributes else 'N/A'
            
            print(f"\n{i:3d}. Feature: {feature.featuretype}")
            print(f"     ID: {feature_id}")
            print(f"     Name: {feature_name}")
            print(f"     Location: {feature.seqid}:{feature.start}-{feature.end} ({feature.strand})")
            print(f"     Hierarchy Level: {level}")
            print(f"     Size: {size:,} bp")
            print(f"     Sequence: {display_seq}")
    
    def process_gff_alignment(self, max_features: int = None) -> List[Tuple[gffutils.Feature, str]]:
        """
        Complete workflow: load data, sort features, extract sequences.
        
        Args:
            max_features: Maximum number of features to process (None for all)
            
        Returns:
            List of (feature, sequence) tuples in topological order
        """
        # Load data
        self.load_gff_database()
        self.load_reference_genome()
        
        # Sort features
        self.topological_sort_features()
        
        # Extract sequences
        features_to_process = self.sorted_features
        if max_features:
            features_to_process = features_to_process[:max_features]
        
        return self.get_feature_sequences(features_to_process)
