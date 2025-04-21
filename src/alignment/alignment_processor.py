"""
Processor for the two-phase feature alignment strategy.
"""
import logging
from typing import Dict, List, Optional, Set, Tuple, Union

import networkx as nx
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from src.parsers.feature_graph import FeatureGraph
from src.parsers.gfa_parser import GFAParser
from src.parsers.gff_parser import GFF3Parser
from src.parsers.fasta_parser import FastaParser
from src.alignment.minimap_wrapper import MinimapAligner

logger = logging.getLogger(__name__)

class AlignmentProcessor:
    """
    Processor for aligning genomic features to paths in a variation graph.
    
    This class implements a two-phase alignment strategy:
    1. First phase: Aligns parent features (e.g., genes) to target paths
    2. Second phase: Aligns child features (e.g., exons) within parent boundaries
    
    It handles feature duplications and multiple alignments.
    """
    
    def __init__(self):
        """Initialize the alignment processor."""
        self.gfa_parser = None
        self.gff_parser = None
        self.fasta_parser = None
        self.feature_graph = FeatureGraph()
        self.aligner = MinimapAligner(preset="splice")
        
        # Results storage
        self.path_sequences = {}  # Path ID -> sequence
        self.aligned_features = {}  # Path ID -> {feature_id: [aligned_features]}
        
    def load_data(self, gfa_file: str, gff_file: str, reference_fasta: str):
        """
        Load all required data for alignment.
        
        Args:
            gfa_file: Path to GFA file
            gff_file: Path to GFF3 file
            reference_fasta: Path to reference FASTA file
        """
        # Load GFA data
        self.gfa_parser = GFAParser()
        self.gfa_parser.parse(gfa_file)
        
        # Load GFF3 data
        self.gff_parser = GFF3Parser()
        self.gff_parser.parse(gff_file)
        
        # Load reference sequences
        self.fasta_parser = FastaParser()
        self.fasta_parser.parse(reference_fasta)
        
        # Build feature relationship graph
        # First call index_features to make sure features are indexed
        self.gff_parser._index_features()
        
        # Collect all features using public methods
        all_types = self.gff_parser.get_all_feature_types()
        all_features = []
        for feat_type in all_types:
            all_features.extend(self.gff_parser.get_features_by_type(feat_type))
            
        # Create a dictionary of features by ID
        features = {f.id: f for f in all_features if hasattr(f, 'id') and f.id}
        self.feature_graph.build_from_features(features)
        
    def extract_path_sequences(self, path_ids: List[str]):
        """
        Extract sequences for specified paths.
        
        Args:
            path_ids: List of path IDs to extract sequences for
        """
        for path_id in path_ids:
            # Get the path from the GFA
            paths = self.gfa_parser.get_paths()
            if path_id not in paths:
                logger.warning(f"Path {path_id} not found in GFA")
                continue
                
            # Extract segment IDs from the path
            path = paths[path_id]
            segments = []
            
            # Handle different GFA versions
            try:
                # GFApy has different representations for GFA1 and GFA2
                if hasattr(path, 'segment_names'):
                    # For string segment names
                    segments = []
                    for seg_name in path.segment_names:
                        if isinstance(seg_name, str):
                            # String segment name
                            segments.append(seg_name.strip('+').strip('-'))
                        elif hasattr(seg_name, 'name'):
                            # GFA1 segment object
                            segments.append(str(seg_name.name))
                        else:
                            # Try converting to string
                            seg_str = str(seg_name)
                            # Remove orientation if present
                            segments.append(seg_str.strip('+').strip('-'))
                elif hasattr(path, 'items'):
                    # Direct access to segment items
                    segments = []
                    for item in path.items:
                        if hasattr(item, 'name'):
                            segments.append(str(item.name).strip('+').strip('-'))
                        else:
                            # Try as string
                            segments.append(str(item).strip('+').strip('-'))
                else:
                    # Try with path.path_name if available (some GFA implementations)
                    segments = []
                    path_str = getattr(path, 'path_name', str(path))
                    if isinstance(path_str, str) and ',' in path_str:
                        # Parse comma-separated path
                        for seg in path_str.split(','):
                            # Remove orientation characters if present
                            segments.append(seg.strip('+').strip('-'))
            except Exception as e:
                logger.warning(f"Error extracting segments from path {path_id}: {e}")
                
            if not segments:
                logger.warning(f"Unknown path format for {path_id}, cannot extract segments")
                continue
                
            if not segments:
                logger.warning(f"No segments found for path {path_id}")
                continue
                
            # Concatenate segment sequences
            path_seq = ""
            for seg_id in segments:
                # Get the segment sequence - handle different possible methods
                seg_seq = None
                try:
                    # Try get_segment_sequence if it exists
                    if hasattr(self.gfa_parser, 'get_segment_sequence'):
                        seg_seq = self.gfa_parser.get_segment_sequence(seg_id)
                    # Try getting the segment and then its sequence
                    elif hasattr(self.gfa_parser, 'get_segments'):
                        segments_dict = self.gfa_parser.get_segments()
                        if seg_id in segments_dict:
                            seg = segments_dict[seg_id]
                            # Try different attribute names for sequence
                            for attr in ['sequence', 'seq', 'sequence_str']:
                                if hasattr(seg, attr):
                                    seg_seq = getattr(seg, attr)
                                    if seg_seq:
                                        break
                    
                    # Try direct GFApy segment access if available
                    if not seg_seq and hasattr(self.gfa_parser.gfa, 'segment'):
                        seg = self.gfa_parser.gfa.segment(seg_id)
                        if seg and hasattr(seg, 'sequence'):
                            seg_seq = seg.sequence
                
                except Exception as e:
                    logger.warning(f"Error getting sequence for segment {seg_id}: {e}")
                
                if seg_seq:
                    path_seq += seg_seq
                else:
                    logger.warning(f"No sequence found for segment {seg_id} in path {path_id}")
            
            if path_seq:
                self.path_sequences[path_id] = path_seq
                logger.info(f"Extracted sequence of length {len(path_seq)} for path {path_id}")
            else:
                logger.warning(f"Failed to extract sequence for path {path_id}")
    
    def align_features_to_paths(self, 
                              path_ids: List[str], 
                              feature_types: Optional[List[str]] = None):
        """
        Align features to specified paths using the two-phase strategy.
        
        Args:
            path_ids: List of path IDs to align features to
            feature_types: Optional list of feature types to align (default: genes)
            
        Returns:
            Dictionary mapping path IDs to dictionaries of feature IDs to lists of aligned features
        """
        if not feature_types:
            feature_types = ["gene"]
        
        # Extract path sequences if not already done
        if not self.path_sequences:
            self.extract_path_sequences(path_ids)
        
        # First phase: Align parent features
        self._align_parent_features(path_ids, feature_types)
        
        # Second phase: Align child features within parent boundaries
        self._align_child_features(path_ids)
        
        return self.aligned_features
    
    def _align_parent_features(self, path_ids: List[str], feature_types: List[str]):
        """
        First phase: Align parent features (e.g., genes) to paths.
        
        Args:
            path_ids: List of path IDs to align to
            feature_types: List of feature types to align
        """
        # Get all features of the specified types
        parent_features = []
        for feature_type in feature_types:
            parent_features.extend(self.gff_parser.get_features_by_type(feature_type))
        
        logger.info(f"Aligning {len(parent_features)} parent features of types {feature_types}")
        
        # For each path
        for path_id in path_ids:
            if path_id not in self.path_sequences:
                logger.warning(f"No sequence available for path {path_id}, skipping")
                continue
                
            path_seq = self.path_sequences[path_id]
            self.aligned_features[path_id] = {}
            
            # Load the path sequence as reference for alignment
            self.aligner.load_reference(path_seq)
            
            # For each parent feature
            for feature in parent_features:
                # Get the feature sequence from the reference
                ref_id = None
                
                # Try to get the sequence ID from qualifiers
                if 'seqid' in feature.qualifiers:
                    ref_id = feature.qualifiers['seqid'][0]
                # If not in qualifiers, try getting from GFF columns
                elif hasattr(feature, 'ref'):
                    ref_id = feature.ref
                elif hasattr(feature, 'seq_id'):
                    ref_id = feature.seq_id
                
                # If we still don't have a ref_id, try to get from the first field of GFF
                if not ref_id and hasattr(feature, 'qualifiers'):
                    # Look at standard GFF columns that might be stored in other qualifiers
                    for possible_key in ['chromosome', 'ref', 'contig', 'seq']:
                        if possible_key in feature.qualifiers:
                            ref_id = feature.qualifiers[possible_key][0]
                            break
                
                if not ref_id:
                    logger.debug(f"No sequence ID found for feature {feature.id}")
                    continue
                    
                ref_seq = self.fasta_parser.get_sequence(ref_id)
                if not ref_seq:
                    logger.warning(f"Reference sequence {ref_id} not found")
                    continue
                
                # Extract the feature sequence
                start = int(feature.location.start)
                end = int(feature.location.end)
                if start >= len(ref_seq.seq) or end > len(ref_seq.seq):
                    logger.warning(f"Feature coordinates ({start}, {end}) out of bounds for reference {ref_id} (length {len(ref_seq.seq)})")
                    continue
                    
                feature_seq = str(ref_seq.seq[start:end])
                
                # Align the feature to the path
                alignments = self.aligner.align_sequence(feature_seq)
                
                if not alignments:
                    logger.debug(f"No alignments found for feature {feature.id} on path {path_id}")
                    continue
                
                # Create aligned features for each alignment
                aligned_features = []
                for aln in alignments:
                    # Create a new feature with coordinates mapped to the path
                    new_feature = SeqFeature(
                        id=feature.id,
                        type=feature.type,
                        location=FeatureLocation(aln.r_st, aln.r_en, strand=feature.location.strand),
                        qualifiers=feature.qualifiers.copy()
                    )
                    
                    # Add additional alignment metadata
                    new_feature.qualifiers['alignment_score'] = [str(aln.mapq)]
                    new_feature.qualifiers['alignment_cigar'] = [aln.cigar_str]
                    new_feature.qualifiers['original_location'] = [f"{ref_id}:{start}-{end}"]
                    
                    aligned_features.append(new_feature)
                
                # Store the aligned features
                self.aligned_features[path_id][feature.id] = aligned_features
                logger.info(f"Aligned feature {feature.id} to path {path_id} ({len(alignments)} alignments)")
    
    def _align_child_features(self, path_ids: List[str]):
        """
        Second phase: Align child features within parent boundaries.
        
        Args:
            path_ids: List of path IDs to align to
        """
        # For each path
        for path_id in path_ids:
            if path_id not in self.aligned_features:
                logger.warning(f"No aligned parent features for path {path_id}, skipping")
                continue
                
            path_seq = self.path_sequences[path_id]
            
            # Create a copy of parent IDs to avoid dictionary modification during iteration
            parent_ids = list(self.aligned_features[path_id].keys())
            
            # For each aligned parent feature
            for parent_id in parent_ids:
                aligned_parents = self.aligned_features[path_id][parent_id]
                # Get the children of this parent
                children_ids = self.feature_graph.get_children(parent_id)
                if not children_ids:
                    logger.info(f"No children found for feature {parent_id}")
                    continue
                    
                # For each aligned instance of the parent
                for aligned_parent in aligned_parents:
                    parent_start = int(aligned_parent.location.start)
                    parent_end = int(aligned_parent.location.end)
                    
                    # Extract the parent region from the path
                    parent_region = path_seq[parent_start:parent_end]
                    
                    # Load this region as the reference for alignment
                    self.aligner.load_reference(parent_region)
                    
                    # For each child feature
                    for child_id in children_ids:
                        child_feature = self.gff_parser.get_feature_by_id(child_id)
                        if not child_feature:
                            logger.warning(f"Child feature {child_id} not found")
                            continue
                            
                        # Get the original child sequence from the reference
                        ref_id = None
                        
                        # Try to get the sequence ID from qualifiers
                        if 'seqid' in child_feature.qualifiers:
                            ref_id = child_feature.qualifiers['seqid'][0]
                        # If not in qualifiers, try getting from GFF columns
                        elif hasattr(child_feature, 'ref'):
                            ref_id = child_feature.ref
                        elif hasattr(child_feature, 'seq_id'):
                            ref_id = child_feature.seq_id
                            
                        # If we still don't have a ref_id, try to get from the first field of GFF
                        if not ref_id and hasattr(child_feature, 'qualifiers'):
                            # Look at standard GFF columns that might be stored in other qualifiers
                            for possible_key in ['chromosome', 'ref', 'contig', 'seq']:
                                if possible_key in child_feature.qualifiers:
                                    ref_id = child_feature.qualifiers[possible_key][0]
                                    break
                        
                        if not ref_id:
                            logger.debug(f"No sequence ID found for feature {child_id}")
                            continue
                            
                        ref_seq = self.fasta_parser.get_sequence(ref_id)
                        if not ref_seq:
                            logger.warning(f"Reference sequence {ref_id} not found")
                            continue
                        
                        # Extract the child feature sequence
                        child_start = int(child_feature.location.start)
                        child_end = int(child_feature.location.end)
                        if child_start >= len(ref_seq.seq) or child_end > len(ref_seq.seq):
                            logger.warning(f"Child feature coordinates ({child_start}, {child_end}) out of bounds")
                            continue
                            
                        child_seq = str(ref_seq.seq[child_start:child_end])
                        
                        # Align the child to the parent region
                        alignments = self.aligner.align_sequence(child_seq, min_score=30)
                        
                        if not alignments:
                            logger.warning(f"No alignments found for child {child_id} within parent {parent_id}")
                            continue
                        
                        # Create aligned features for each alignment
                        for aln in alignments:
                            # Map the child coordinates to the parent region, then to the path
                            path_child_start = parent_start + aln.r_st
                            path_child_end = parent_start + aln.r_en
                            
                            # Create a new feature with coordinates mapped to the path
                            new_child = SeqFeature(
                                id=child_id,
                                type=child_feature.type,
                                location=FeatureLocation(path_child_start, path_child_end, 
                                                        strand=child_feature.location.strand),
                                qualifiers=child_feature.qualifiers.copy()
                            )
                            
                            # Add additional alignment metadata
                            new_child.qualifiers['alignment_score'] = [str(aln.mapq)]
                            new_child.qualifiers['alignment_cigar'] = [aln.cigar_str]
                            new_child.qualifiers['original_location'] = [f"{ref_id}:{child_start}-{child_end}"]
                            new_child.qualifiers['parent_feature'] = [parent_id]
                            
                            # Store the aligned child feature
                            if child_id not in self.aligned_features[path_id]:
                                self.aligned_features[path_id][child_id] = []
                            
                            self.aligned_features[path_id][child_id].append(new_child)
                            
                        logger.info(f"Aligned child {child_id} within parent {parent_id} on path {path_id} ({len(alignments)} alignments)")
    
    def get_alignment_results(self) -> Dict[str, Dict[str, List[SeqFeature]]]:
        """
        Get all alignment results.
        
        Returns:
            Dictionary mapping path IDs to dictionaries of feature IDs to lists of aligned features
        """
        return self.aligned_features
    
    def get_path_alignments(self, path_id: str) -> Dict[str, List[SeqFeature]]:
        """
        Get alignment results for a specific path.
        
        Args:
            path_id: Path ID to get alignments for
            
        Returns:
            Dictionary mapping feature IDs to lists of aligned features for the specified path
        """
        return self.aligned_features.get(path_id, {})
    
    def get_feature_alignments(self, path_id: str, feature_id: str) -> List[SeqFeature]:
        """
        Get alignment results for a specific feature on a specific path.
        
        Args:
            path_id: Path ID to get alignments for
            feature_id: Feature ID to get alignments for
            
        Returns:
            List of aligned features for the specified feature on the specified path
        """
        return self.aligned_features.get(path_id, {}).get(feature_id, [])
    
    def export_alignments_tsv(self, output_file: str):
        """
        Export alignment results to a TSV file.
        
        Args:
            output_file: Path to output file
        """
        with open(output_file, 'w') as f:
            # Write header
            f.write("path_id\tfeature_id\ttype\tstart\tend\tstrand\tparent\tscore\tcigar\toriginal_location\n")
            
            # Write alignments
            for path_id, features in self.aligned_features.items():
                for feature_id, alignments in features.items():
                    for aln in alignments:
                        # Get feature info
                        feature_type = aln.type
                        start = int(aln.location.start)
                        end = int(aln.location.end)
                        strand = '+' if aln.location.strand == 1 else '-' if aln.location.strand == -1 else '.'
                        
                        # Get alignment info
                        parent = aln.qualifiers.get('parent_feature', ['.'])[0]
                        score = aln.qualifiers.get('alignment_score', [''])[0]
                        cigar = aln.qualifiers.get('alignment_cigar', [''])[0]
                        orig_loc = aln.qualifiers.get('original_location', [''])[0]
                        
                        # Write line
                        f.write(f"{path_id}\t{feature_id}\t{feature_type}\t{start}\t{end}\t{strand}\t{parent}\t{score}\t{cigar}\t{orig_loc}\n")
            
        logger.info(f"Exported alignment results to {output_file}")
