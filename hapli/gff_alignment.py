#!/usr/bin/env python3
"""
GFF3 topological sorting and sequence alignment functionality for hapli.

This module provides core functionality to:
1. Topologically sort GFF3 features by hierarchy and size
2. Extract corresponding sequences from reference genome
3. Align features to GFA paths using minimap2
4. Handle hierarchical constraints for sub-feature alignment
"""

import logging
import json
import subprocess
import tempfile
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, NamedTuple, Union
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

try:
    import mappy
    MAPPY_AVAILABLE = True
except ImportError:
    MAPPY_AVAILABLE = False
    logging.warning("mappy not available, falling back to minimap2 subprocess")

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False
    logging.warning("pysam not available, BAM output will use samtools subprocess")


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


@dataclass
class FeatureAlignment:
    """Represents an alignment of a GFF3 feature to a genome graph path."""
    feature_id: str
    feature_type: str
    path_name: str
    path_start: int
    path_end: int
    feature_start: int
    feature_end: int
    strand: str
    mapq: int
    alignment_score: int
    sequence_identity: float
    feature_sequence: str
    parent_feature_id: Optional[str] = None
    hierarchy_level: int = 0
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'feature_id': self.feature_id,
            'feature_type': self.feature_type,
            'path_name': self.path_name,
            'path_start': self.path_start,
            'path_end': self.path_end,
            'feature_start': self.feature_start,
            'feature_end': self.feature_end,
            'strand': self.strand,
            'mapq': self.mapq,
            'alignment_score': self.alignment_score,
            'sequence_identity': self.sequence_identity,
            'feature_sequence': self.feature_sequence,
            'parent_feature_id': self.parent_feature_id,
            'hierarchy_level': self.hierarchy_level
        }
    
    def to_sam_line(self, reference_length: int) -> str:
        """Convert alignment to SAM format line."""
        # SAM fields: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
        qname = self.feature_id
        flag = 16 if self.strand == '-' else 0  # 16 = reverse strand
        rname = self.path_name
        pos = self.path_start + 1  # SAM uses 1-based coordinates
        mapq = self.mapq
        
        # Simple CIGAR - assume match for now (could be improved with actual alignment)
        cigar = f"{len(self.feature_sequence)}M"
        
        rnext = "*"
        pnext = 0
        tlen = 0
        seq = self.feature_sequence
        qual = "*"  # No quality scores available
        
        # Add optional tags
        optional_tags = [
            f"AS:i:{self.alignment_score}",
            f"NM:i:0",  # Edit distance (simplified)
            f"XS:f:{self.sequence_identity}",  # Sequence identity
            f"XL:i:{self.hierarchy_level}",  # Hierarchy level (custom tag)
        ]
        
        if self.parent_feature_id:
            optional_tags.append(f"XP:Z:{self.parent_feature_id}")  # Parent feature ID
        
        sam_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\t" + "\t".join(optional_tags)
        return sam_line


@dataclass
class FeatureInfo:
    """Pre-computed information about a feature for thread-safe processing."""
    feature_id: str
    feature_type: str
    sequence: str
    hierarchy_level: int
    parent_feature_id: Optional[str]
    seqid: str
    start: int
    end: int
    strand: str


class GenomeGraphPathExtractor:
    """Extract paths from genome graph files (GFA/VG/XG)."""
    
    def __init__(self, graph_file: Path, reference_path_name: Optional[str] = None, vg_threads: int = 8):
        """Initialize with genome graph file."""
        self.graph_file = graph_file
        self.reference_path_name = reference_path_name
        self.vg_threads = vg_threads
        self.file_type = self._detect_file_type()
    
    def _detect_file_type(self) -> str:
        """Detect the type of genome graph file."""
        suffix = self.graph_file.suffix.lower()
        if suffix == '.gfa':
            return 'gfa'
        elif suffix == '.vg':
            return 'vg'
        elif suffix == '.xg':
            return 'xg'
        else:
            # Try to detect from content
            try:
                with open(self.graph_file, 'r') as f:
                    first_line = f.readline().strip()
                    if first_line.startswith('H\t'):
                        return 'gfa'
            except UnicodeDecodeError:
                # Binary file, likely VG or XG
                pass
            
            logging.warning(f"Could not detect file type for {self.graph_file}, assuming VG")
            return 'vg'
    
    def extract_paths(self) -> Dict[str, str]:
        """Extract all paths from the genome graph."""
        if self.file_type == 'gfa':
            return self._extract_gfa_paths()
        elif self.file_type in ['vg', 'xg']:
            return self._extract_vg_paths()
        else:
            raise ValueError(f"Unsupported file type: {self.file_type}")
    
    def _extract_gfa_paths(self) -> Dict[str, str]:
        """Extract paths from GFA file."""
        logging.info(f"Extracting paths from GFA file: {self.graph_file}")
        
        # Parse GFA to extract segments and paths
        segments = {}
        paths = {}
        
        with open(self.graph_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                fields = line.split('\t')
                if not fields:
                    continue
                
                if fields[0] == 'S':  # Segment
                    if len(fields) >= 3:
                        seg_id = fields[1]
                        sequence = fields[2]
                        segments[seg_id] = sequence
                
                elif fields[0] == 'P':  # Path
                    if len(fields) >= 3:
                        path_name = fields[1]
                        path_segments = fields[2].split(',')
                        
                        # Reconstruct sequence from segments
                        path_sequence = self._reconstruct_path_sequence(path_segments, segments, path_name)
                        if path_sequence:
                            paths[path_name] = path_sequence
                
                elif fields[0] == 'W':  # Walk (alternative path representation)
                    if len(fields) >= 6:
                        # W format: W <sample_id> <haplotype_index> <seq_id> <seq_start> <seq_end> <walk>
                        sample_id = fields[1]
                        haplotype_index = fields[2]
                        seq_id = fields[3]
                        seq_start = fields[4]
                        seq_end = fields[5]
                        walk = fields[6]
                        
                        # Create path name from walk components
                        path_name = f"{sample_id}#{haplotype_index}#{seq_id}#{seq_start}"
                        
                        # Parse walk to get segments
                        path_segments = []
                        if walk:
                            # Walk format can be like: >1<2>3 or 1+,2-,3+
                            walk_segments = self._parse_walk(walk)
                            path_sequence = self._reconstruct_path_sequence(walk_segments, segments, path_name)
                            if path_sequence:
                                paths[path_name] = path_sequence
        
        logging.info(f"Extracted {len(paths)} paths from GFA")
        if len(paths) == 0:
            logging.warning("No paths found in GFA file. Check if the file contains P or W lines.")
            logging.debug(f"Found {len(segments)} segments")
        
        return paths
    
    def _parse_walk(self, walk: str) -> List[str]:
        """Parse a walk string to extract segment IDs with orientations."""
        segments = []
        
        # Handle different walk formats
        if ',' in walk:
            # Format like: 1+,2-,3+
            for seg in walk.split(','):
                seg = seg.strip()
                if seg:
                    segments.append(seg)
        else:
            # Format like: >1<2>3 or other formats
            i = 0
            current_seg = ""
            orientation = "+"
            
            while i < len(walk):
                char = walk[i]
                if char == '>' or char == '<':
                    if current_seg:
                        # Add previous segment
                        segments.append(current_seg + orientation)
                        current_seg = ""
                    orientation = "+" if char == '>' else "-"
                elif char.isdigit() or char.isalpha() or char == '_':
                    current_seg += char
                elif char in [' ', '\t']:
                    if current_seg:
                        segments.append(current_seg + orientation)
                        current_seg = ""
                        orientation = "+"
                i += 1
            
            # Add final segment
            if current_seg:
                segments.append(current_seg + orientation)
        
        return segments
    
    def _reconstruct_path_sequence(self, path_segments: List[str], segments: Dict[str, str], path_name: str) -> str:
        """Reconstruct sequence from path segments."""
        path_sequence = ""
        
        for seg in path_segments:
            if not seg:
                continue
            
            # Handle orientation
            if seg.endswith('+'):
                seg_id = seg[:-1]
                if seg_id in segments:
                    path_sequence += segments[seg_id]
                else:
                    logging.debug(f"Segment {seg_id} not found for path {path_name}")
            elif seg.endswith('-'):
                seg_id = seg[:-1]
                if seg_id in segments:
                    # Reverse complement
                    from Bio.Seq import Seq
                    path_sequence += str(Seq(segments[seg_id]).reverse_complement())
                else:
                    logging.debug(f"Segment {seg_id} not found for path {path_name}")
            else:
                # No orientation specified, assume forward
                if seg in segments:
                    path_sequence += segments[seg]
                else:
                    logging.debug(f"Segment {seg} not found for path {path_name}")
        
        return path_sequence
    
    def _extract_vg_paths(self) -> Dict[str, str]:
        """Extract paths from VG/XG file using vg commands."""
        logging.info(f"Extracting paths from VG/XG file: {self.graph_file}")
        
        try:
            # Get list of paths
            cmd = ['vg', 'paths', '-L', str(self.graph_file)]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            path_names = [line.strip() for line in result.stdout.split('\n') if line.strip()]
            
            logging.info(f"Found {len(path_names)} paths in graph")
            
            # Extract sequences for each path
            paths = {}
            for path_name in path_names:
                try:
                    cmd = ['vg', 'find', '-p', path_name, '-x', str(self.graph_file), '-t', str(self.vg_threads)]
                    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                    
                    # Parse FASTA output
                    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
                        tmp.write(result.stdout)
                        tmp_path = tmp.name
                    
                    try:
                        with open(tmp_path, 'r') as f:
                            for record in SeqIO.parse(f, 'fasta'):
                                paths[path_name] = str(record.seq)
                                break  # Take first sequence
                    finally:
                        Path(tmp_path).unlink(missing_ok=True)
                        
                except subprocess.CalledProcessError as e:
                    logging.warning(f"Could not extract path {path_name}: {e}")
                    continue
            
            return paths
            
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to extract paths from VG file: {e}")
            return {}


class GFFAligner:
    """Handle GFF3 topological sorting and sequence alignment."""
    
    def __init__(self, gff_path: Path, reference_path: Path, graph_file: Optional[Path] = None, 
                 reference_path_name: Optional[str] = None, vg_threads: int = 8):
        """Initialize with GFF3 and reference files."""
        self.gff_path = gff_path
        self.reference_path = reference_path
        self.graph_file = graph_file
        self.reference_path_name = reference_path_name
        self.vg_threads = vg_threads
        self.db = None
        self.reference_genome = {}
        self.sorted_features = []
        self.graph_paths = {}
        self.feature_parent_map = {}
        self.parent_alignments = {}  # Store alignments of parent features
        
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
                return self.db
        except (gffutils.exceptions.FeatureDBError, ValueError) as e:
            logging.debug(f"Could not load existing database: {e}")
        
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
    
    def load_graph_paths(self) -> Dict[str, str]:
        """Load paths from genome graph file."""
        if not self.graph_file:
            return {}
        
        extractor = GenomeGraphPathExtractor(self.graph_file, self.reference_path_name, self.vg_threads)
        self.graph_paths = extractor.extract_paths()
        
        logging.info(f"Loaded {len(self.graph_paths)} graph paths")
        return self.graph_paths
    
    def build_feature_parent_map(self) -> Dict[str, str]:
        """Build mapping from feature ID to parent feature ID."""
        if not self.db:
            raise ValueError("GFF3 database not loaded")
        
        parent_map = {}
        
        for feature in self.db.all_features():
            try:
                parents = list(self.db.parents(feature))
                if parents:
                    parent_map[feature.id] = parents[0].id
            except (gffutils.exceptions.FeatureNotFoundError, StopIteration):
                continue
        
        self.feature_parent_map = parent_map
        return parent_map
    
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
    
    def prepare_feature_info(self, features: List[gffutils.Feature]) -> List[FeatureInfo]:
        """
        Pre-compute all necessary feature information for thread-safe processing.
        
        Args:
            features: List of features to process
            
        Returns:
            List of FeatureInfo objects with pre-computed data
        """
        logging.info(f"Preparing feature information for {len(features)} features")
        
        feature_infos = []
        
        for feature in features:
            try:
                # Extract sequence
                sequence = self.extract_feature_sequence(feature)
                
                # Get hierarchy level
                hierarchy_level = self.get_feature_hierarchy_level(feature)
                
                # Get parent ID
                parent_id = self.feature_parent_map.get(feature.id)
                
                # Create FeatureInfo object
                feature_info = FeatureInfo(
                    feature_id=feature.id,
                    feature_type=feature.featuretype,
                    sequence=sequence,
                    hierarchy_level=hierarchy_level,
                    parent_feature_id=parent_id,
                    seqid=feature.seqid,
                    start=feature.start,
                    end=feature.end,
                    strand=feature.strand
                )
                
                feature_infos.append(feature_info)
                
            except Exception as e:
                logging.warning(f"Could not prepare info for feature {feature.id}: {e}")
                continue
        
        logging.info(f"Successfully prepared {len(feature_infos)} feature infos")
        return feature_infos
    
    def align_sequence_to_path(self, sequence: str, path_sequence: str, path_name: str) -> List[FeatureAlignment]:
        """Align a feature sequence to a graph path using minimap2/mappy."""
        alignments = []
        
        if MAPPY_AVAILABLE:
            # Use mappy for alignment
            try:
                # Create aligner
                aligner = mappy.Aligner(seq=path_sequence, preset="map-ont")
                
                # Perform alignment
                for hit in aligner.map(sequence):
                    alignment = FeatureAlignment(
                        feature_id="",  # Will be filled by caller
                        feature_type="",  # Will be filled by caller
                        path_name=path_name,
                        path_start=hit.r_st,
                        path_end=hit.r_en,
                        feature_start=hit.q_st,
                        feature_end=hit.q_en,
                        strand="+" if hit.strand == 1 else "-",
                        mapq=hit.mapq,
                        alignment_score=hit.score if hasattr(hit, 'score') else 0,
                        sequence_identity=hit.mlen / max(hit.blen, 1),
                        feature_sequence=sequence
                    )
                    alignments.append(alignment)
                    
            except Exception as e:
                logging.warning(f"mappy alignment failed for path {path_name}: {e}")
        
        else:
            # Fall back to minimap2 subprocess
            try:
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as path_tmp:
                    path_tmp.write(f">{path_name}\n{path_sequence}\n")
                    path_tmp_path = path_tmp.name
                
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as query_tmp:
                    query_tmp.write(f">query\n{sequence}\n")
                    query_tmp_path = query_tmp.name
                
                try:
                    # Run minimap2
                    cmd = ['minimap2', '-a', path_tmp_path, query_tmp_path]
                    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                    
                    # Parse SAM output
                    for line in result.stdout.split('\n'):
                        if line.strip() and not line.startswith('@'):
                            fields = line.split('\t')
                            if len(fields) >= 11 and fields[2] != '*':  # Valid alignment
                                alignment = FeatureAlignment(
                                    feature_id="",  # Will be filled by caller
                                    feature_type="",  # Will be filled by caller
                                    path_name=path_name,
                                    path_start=int(fields[3]) - 1,  # Convert to 0-based
                                    path_end=int(fields[3]) - 1 + len(sequence),  # Approximate
                                    feature_start=0,
                                    feature_end=len(sequence),
                                    strand="+" if int(fields[1]) & 16 == 0 else "-",
                                    mapq=int(fields[4]),
                                    alignment_score=0,  # Not easily available from SAM
                                    sequence_identity=1.0,  # Approximate
                                    feature_sequence=sequence
                                )
                                alignments.append(alignment)
                
                finally:
                    Path(path_tmp_path).unlink(missing_ok=True)
                    Path(query_tmp_path).unlink(missing_ok=True)
                    
            except Exception as e:
                logging.warning(f"minimap2 alignment failed for path {path_name}: {e}")
        
        return alignments
    
    def get_parent_constraint_regions(self, feature_info: FeatureInfo) -> Dict[str, List[Tuple[int, int]]]:
        """Get allowed regions for this feature based on parent alignments."""
        if not feature_info.parent_feature_id:
            return {}  # No parent, no constraints
        
        parent_id = feature_info.parent_feature_id
        if parent_id not in self.parent_alignments:
            return {}  # Parent not aligned yet
        
        # Get allowed regions from parent alignments
        constraint_regions = defaultdict(list)
        
        for alignment in self.parent_alignments[parent_id]:
            path_name = alignment.path_name
            # Add some buffer around parent region
            buffer = 100
            start = max(0, alignment.path_start - buffer)
            end = alignment.path_end + buffer
            constraint_regions[path_name].append((start, end))
        
        return constraint_regions
    
    def is_alignment_within_constraints(self, alignment: FeatureAlignment, 
                                      constraint_regions: Dict[str, List[Tuple[int, int]]]) -> bool:
        """Check if alignment falls within constraint regions."""
        if not constraint_regions:
            return True  # No constraints
        
        path_name = alignment.path_name
        if path_name not in constraint_regions:
            return False  # Path not allowed
        
        # Check if alignment overlaps with any allowed region
        for start, end in constraint_regions[path_name]:
            if alignment.path_start < end and alignment.path_end > start:
                return True
        
        return False
    
    def align_feature_to_all_paths(self, feature_info: FeatureInfo, graph_paths: Dict[str, str], 
                                  constraint_regions: Dict[str, List[Tuple[int, int]]]) -> List[FeatureAlignment]:
        """Align a feature to all graph paths with hierarchical constraints (thread-safe)."""
        if not graph_paths:
            return []
        
        all_alignments = []
        
        # Align to each path
        for path_name, path_sequence in graph_paths.items():
            try:
                alignments = self.align_sequence_to_path(feature_info.sequence, path_sequence, path_name)
                
                for alignment in alignments:
                    # Fill in feature information
                    alignment.feature_id = feature_info.feature_id
                    alignment.feature_type = feature_info.feature_type
                    alignment.hierarchy_level = feature_info.hierarchy_level
                    alignment.parent_feature_id = feature_info.parent_feature_id
                    
                    # Check constraints
                    if self.is_alignment_within_constraints(alignment, constraint_regions):
                        all_alignments.append(alignment)
                
            except Exception as e:
                logging.warning(f"Failed to align feature {feature_info.feature_id} to path {path_name}: {e}")
                continue
        
        return all_alignments
    
    def align_features_parallel(self, max_workers: int = 4) -> List[FeatureAlignment]:
        """Align all features to graph paths in parallel, respecting hierarchy."""
        if not self.sorted_features:
            logging.info("Features not sorted. Running pre-alignment setup.")
            if not self.db:
                self.load_gff_database()
            if not self.reference_genome:
                self.load_reference_genome()
            self.build_feature_parent_map()
            self.topological_sort_features()

        if not self.graph_paths:
            logging.info("Graph paths not loaded. Loading now.")
            self.load_graph_paths()
        
        logging.info(f"Aligning {len(self.sorted_features)} features using {max_workers} workers")
        
        all_alignments = []
        
        # Group features by hierarchy level to maintain ordering constraints
        features_by_level = defaultdict(list)
        for feature in self.sorted_features:
            level = self.get_feature_hierarchy_level(feature)
            features_by_level[level].append(feature)
        
        # Process each level sequentially (but features within level in parallel)
        for level in sorted(features_by_level.keys()):
            level_features = features_by_level[level]
            logging.info(f"Processing hierarchy level {level}: {len(level_features)} features")
            
            # Pre-compute all feature information for this level (thread-safe)
            level_feature_infos = self.prepare_feature_info(level_features)
            
            # Process features at this level in parallel
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                # Submit alignment tasks
                future_to_feature_info = {}
                for feature_info in level_feature_infos:
                    # Get parent constraint regions
                    constraint_regions = self.get_parent_constraint_regions(feature_info)
                    
                    # Submit task with all necessary data
                    future = executor.submit(
                        self.align_feature_to_all_paths, 
                        feature_info, 
                        self.graph_paths, 
                        constraint_regions
                    )
                    future_to_feature_info[future] = feature_info
                
                # Collect results
                for future in as_completed(future_to_feature_info):
                    feature_info = future_to_feature_info[future]
                    try:
                        feature_alignments = future.result()
                        all_alignments.extend(feature_alignments)
                        
                        # Update parent alignments for next level
                        if feature_info.feature_id not in self.parent_alignments and feature_alignments:
                            self.parent_alignments[feature_info.feature_id] = feature_alignments
                        
                        logging.debug(f"Feature {feature_info.feature_id}: {len(feature_alignments)} alignments")
                        
                    except Exception as e:
                        logging.error(f"Error aligning feature {feature_info.feature_id}: {e}")
                        continue
        
        logging.info(f"Alignment complete: {len(all_alignments)} total alignments")
        return all_alignments
    
    def write_alignments_json(self, alignments: List[FeatureAlignment], output_path: Path) -> None:
        """Write alignments to JSON file."""
        logging.info(f"Writing {len(alignments)} alignments to JSON: {output_path}")
        
        # Convert to serializable format
        alignment_data = {
            'alignments': [alignment.to_dict() for alignment in alignments],
            'summary': {
                'total_alignments': len(alignments),
                'features_aligned': len(set(a.feature_id for a in alignments)),
                'paths_used': len(set(a.path_name for a in alignments)),
                'hierarchy_levels': len(set(a.hierarchy_level for a in alignments))
            }
        }
        
        with open(output_path, 'w') as f:
            json.dump(alignment_data, f, indent=2)
        
        logging.info(f"JSON alignments written to {output_path}")
    
    def write_alignments_gam(self, alignments: List[FeatureAlignment], output_path: Path, vg_threads: int = 8, keep_temp_files: bool = False) -> None:
        """Write alignments to GAM file using vg tools."""
        logging.info(f"Writing {len(alignments)} alignments to GAM: {output_path}")
        
        if not self.graph_file:
            raise ValueError("Graph file required for GAM output")
        
        # Create temporary SAM file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as sam_tmp:
            sam_tmp_path = sam_tmp.name
            
            # Write SAM header
            sam_tmp.write("@HD\tVN:1.0\tSO:unsorted\n")
            
            # Write reference sequences (paths)
            for path_name, path_sequence in self.graph_paths.items():
                sam_tmp.write(f"@SQ\tSN:{path_name}\tLN:{len(path_sequence)}\n")
            
            # Write alignments
            for alignment in alignments:
                path_sequence = self.graph_paths.get(alignment.path_name, "")
                reference_length = len(path_sequence)
                sam_line = alignment.to_sam_line(reference_length)
                sam_tmp.write(sam_line + "\n")
        
        try:
            logging.debug(f"Created temporary SAM file: {sam_tmp_path}")
            
            # Check SAM file size and content for debugging
            sam_size = Path(sam_tmp_path).stat().st_size
            logging.debug(f"SAM file size: {sam_size} bytes")
            
            # For debugging: show first few lines of SAM file
            if logging.getLogger().level <= logging.DEBUG:
                with open(sam_tmp_path, 'r') as f:
                    sample_lines = []
                    for i, line in enumerate(f):
                        if i < 10:  # First 10 lines
                            sample_lines.append(line.strip())
                        else:
                            break
                    logging.debug(f"Sample SAM content:\n" + "\n".join(sample_lines))
            
            # Convert SAM to GAM using vg
            if self.graph_file.suffix.lower() == '.gfa':
                # For GFA files, convert GFA to VG first, then use VG for GAM conversion
                with tempfile.NamedTemporaryFile(suffix='.vg', delete=False) as vg_tmp:
                    vg_tmp_path = vg_tmp.name
                
                try:
                    # Convert GFA to VG
                    logging.debug(f"Converting GFA to VG: {self.graph_file} -> {vg_tmp_path}")
                    cmd_gfa_to_vg = ['vg', 'view', '-Fv', '--threads', str(vg_threads), str(self.graph_file)]
                    with open(vg_tmp_path, 'wb') as vg_out:
                        result = subprocess.run(cmd_gfa_to_vg, stdout=vg_out, stderr=subprocess.PIPE, check=True)
                    
                    vg_size = Path(vg_tmp_path).stat().st_size
                    logging.debug(f"VG file size: {vg_size} bytes")
                    
                    # Convert SAM to GAM using the VG file
                    logging.debug(f"Converting SAM to GAM using VG file")
                    cmd = ['vg', 'inject', '-x', vg_tmp_path, '-t', str(vg_threads), sam_tmp_path]
                    with open(output_path, 'wb') as gam_out:
                        result = subprocess.run(cmd, stdout=gam_out, stderr=subprocess.PIPE, check=True)
                    
                finally:
                    # Clean up temporary VG file
                    if not keep_temp_files:
                        Path(vg_tmp_path).unlink(missing_ok=True)
                    else:
                        logging.info(f"Kept temporary VG file: {vg_tmp_path}")
            
            else:
                # For VG/XG files, convert SAM to GAM directly
                logging.debug(f"Converting SAM to GAM using VG/XG file")
                cmd = ['vg', 'inject', '-x', str(self.graph_file), '-t', str(vg_threads), sam_tmp_path]
                with open(output_path, 'wb') as gam_out:
                    result = subprocess.run(cmd, stdout=gam_out, stderr=subprocess.PIPE, check=True)
            
            # Check GAM file size
            gam_size = Path(output_path).stat().st_size
            logging.info(f"GAM file size: {gam_size} bytes")
            
            if gam_size < 100:  # Suspiciously small
                logging.warning(f"GAM file seems very small ({gam_size} bytes) for {len(alignments)} alignments")
            
            logging.info(f"GAM alignments written to {output_path}")
            
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to create GAM file: {e}")
            if e.stderr:
                stderr_text = e.stderr.decode() if isinstance(e.stderr, bytes) else str(e.stderr)
                logging.error(f"stderr: {stderr_text}")
            raise
        
        finally:
            # Clean up temporary SAM file
            if not keep_temp_files:
                Path(sam_tmp_path).unlink(missing_ok=True)
            else:
                logging.info(f"Kept temporary SAM file: {sam_tmp_path}")
    
    def write_alignments_bam(self, alignments: List[FeatureAlignment], output_path: Path) -> None:
        """Write alignments to separate BAM files for each path."""
        logging.info(f"Writing {len(alignments)} alignments to BAM files")
        
        # Group alignments by path
        alignments_by_path = defaultdict(list)
        for alignment in alignments:
            alignments_by_path[alignment.path_name].append(alignment)
        
        output_dir = output_path.parent
        output_base = output_path.stem
        
        bam_files_created = []
        
        for path_name, path_alignments in alignments_by_path.items():
            # Sanitize path name for filename
            safe_path_name = path_name.replace('/', '_').replace('#', '_').replace(':', '_')
            bam_file = output_dir / f"{output_base}_{safe_path_name}.bam"
            
            try:
                self._write_single_bam(path_alignments, bam_file, path_name)
                bam_files_created.append(bam_file)
                
            except Exception as e:
                logging.error(f"Failed to create BAM file for path {path_name}: {e}")
                continue
        
        logging.info(f"Created {len(bam_files_created)} BAM files in {output_dir}")
        
        # Create index files for BAM files
        for bam_file in bam_files_created:
            try:
                self._index_bam_file(bam_file)
            except Exception as e:
                logging.warning(f"Failed to index BAM file {bam_file}: {e}")
    
    def _write_single_bam(self, alignments: List[FeatureAlignment], bam_file: Path, path_name: str) -> None:
        """Write alignments for a single path to BAM file."""
        if not alignments:
            return
        
        path_sequence = self.graph_paths.get(path_name, "")
        reference_length = len(path_sequence)
        
        if PYSAM_AVAILABLE:
            # Use pysam to write BAM
            try:
                # Create header
                header = {
                    'HD': {'VN': '1.0', 'SO': 'unsorted'},
                    'SQ': [{'LN': reference_length, 'SN': path_name}]
                }
                
                with pysam.AlignmentFile(str(bam_file), "wb", header=header) as bam_out:
                    for alignment in alignments:
                        # Create pysam AlignedSegment
                        read = pysam.AlignedSegment()
                        read.query_name = alignment.feature_id
                        read.query_sequence = alignment.feature_sequence
                        read.flag = 16 if alignment.strand == '-' else 0
                        read.reference_id = 0  # First (and only) reference
                        read.reference_start = alignment.path_start
                        read.mapping_quality = alignment.mapq
                        read.cigar = [(0, len(alignment.feature_sequence))]  # Match
                        read.next_reference_id = -1
                        read.next_reference_start = -1
                        read.template_length = 0
                        read.query_qualities = None
                        
                        # Add tags
                        read.set_tag('AS', alignment.alignment_score, 'i')
                        read.set_tag('XS', alignment.sequence_identity, 'f')
                        read.set_tag('XL', alignment.hierarchy_level, 'i')
                        if alignment.parent_feature_id:
                            read.set_tag('XP', alignment.parent_feature_id, 'Z')
                        
                        bam_out.write(read)
                
                logging.debug(f"Created BAM file {bam_file} with {len(alignments)} alignments")
                
            except Exception as e:
                logging.error(f"pysam failed to write BAM file {bam_file}: {e}")
                # Fall back to samtools
                self._write_bam_with_samtools(alignments, bam_file, path_name, reference_length)
        
        else:
            # Fall back to samtools
            self._write_bam_with_samtools(alignments, bam_file, path_name, reference_length)
    
    def _write_bam_with_samtools(self, alignments: List[FeatureAlignment], bam_file: Path, 
                                path_name: str, reference_length: int) -> None:
        """Write BAM file using samtools subprocess."""
        # Create temporary SAM file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as sam_tmp:
            sam_tmp_path = sam_tmp.name
            
            # Write SAM header
            sam_tmp.write("@HD\tVN:1.0\tSO:unsorted\n")
            sam_tmp.write(f"@SQ\tSN:{path_name}\tLN:{reference_length}\n")
            
            # Write alignments
            for alignment in alignments:
                sam_line = alignment.to_sam_line(reference_length)
                sam_tmp.write(sam_line + "\n")
        
        try:
            # Convert SAM to BAM using samtools
            cmd = ['samtools', 'view', '-b', '-o', str(bam_file), sam_tmp_path]
            subprocess.run(cmd, check=True, capture_output=True)
            
            logging.debug(f"Created BAM file {bam_file} with {len(alignments)} alignments using samtools")
            
        except subprocess.CalledProcessError as e:
            logging.error(f"samtools failed to create BAM file {bam_file}: {e}")
            raise
        
        finally:
            # Clean up temporary file
            Path(sam_tmp_path).unlink(missing_ok=True)
    
    def _index_bam_file(self, bam_file: Path) -> None:
        """Create index for BAM file."""
        try:
            if PYSAM_AVAILABLE:
                pysam.index(str(bam_file))
            else:
                cmd = ['samtools', 'index', str(bam_file)]
                subprocess.run(cmd, check=True, capture_output=True)
            
            logging.debug(f"Created index for {bam_file}")
            
        except Exception as e:
            logging.warning(f"Failed to index {bam_file}: {e}")
    
    def write_alignments(self, alignments: List[FeatureAlignment], output_path: Path, 
                        output_format: str = "json", vg_threads: int = 8, keep_temp_files: bool = False) -> None:
        """Write alignments in the specified format."""
        if output_format == "json":
            self.write_alignments_json(alignments, output_path)
        elif output_format == "gam":
            self.write_alignments_gam(alignments, output_path, vg_threads=vg_threads, keep_temp_files=keep_temp_files)
        elif output_format == "bam":
            self.write_alignments_bam(alignments, output_path)
        else:
            raise ValueError(f"Unsupported output format: {output_format}")
    
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
    
    def process_graph_alignment(self, output_path: Path, output_format: str = "json", 
                               max_workers: int = 4) -> List[FeatureAlignment]:
        """
        Complete workflow: load data, sort features, align to graph paths.
        
        Args:
            output_path: Path to write alignment results
            output_format: Output format ("json", "gam", "bam")
            max_workers: Number of parallel workers
            
        Returns:
            List of feature alignments
        """
        # Load all data
        self.load_gff_database()
        self.load_reference_genome()
        self.load_graph_paths()
        
        # Build feature hierarchy
        self.build_feature_parent_map()
        
        # Sort features topologically
        self.topological_sort_features()
        
        # Perform alignments
        alignments = self.align_features_parallel(max_workers=max_workers)
        
        # Write results
        self.write_alignments(alignments, output_path, output_format)
        
        return alignments
