#!/usr/bin/env python3
"""
Impact detector for analyzing feature alignments against pangenome graphs.

This module provides functionality to analyze GAM alignment data and determine
the impact of structural variations on genomic features.
"""

import logging
import re
from collections import defaultdict, deque, Counter
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any, Union
import json
import gffutils

logger = logging.getLogger(__name__)


class GFAParser:
    """Parse GFA (Graphical Fragment Assembly) files to extract graph structure."""
    
    def __init__(self, gfa_file: Path):
        """
        Initialize GFA parser.
        
        Args:
            gfa_file: Path to GFA file
        """
        self.gfa_file = Path(gfa_file)
        self.segments = {}  # node_id -> sequence
        self.links = defaultdict(list)  # node_id -> [(target_node, orientation)]
        self.paths = {}  # path_name -> list of nodes
        self.node_lengths = {}  # node_id -> length
        self._validate_input()
    
    def _validate_input(self) -> None:
        """Validate input file exists."""
        if not self.gfa_file.exists():
            raise FileNotFoundError(f"GFA file not found: {self.gfa_file}")
    
    def parse(self) -> None:
        """Parse the GFA file to extract graph structure."""
        logger.info(f"Parsing GFA file: {self.gfa_file}")
        
        with open(self.gfa_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                try:
                    self._parse_line(line)
                except Exception as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue
        
        logger.info(f"Parsed {len(self.segments)} segments, "
                   f"{sum(len(links) for links in self.links.values())} links, "
                   f"{len(self.paths)} paths")
    
    def _parse_line(self, line: str) -> None:
        """Parse a single GFA line."""
        fields = line.split('\t')
        if not fields:
            return
        
        record_type = fields[0]
        
        if record_type == 'S':  # Segment
            self._parse_segment(fields)
        elif record_type == 'L':  # Link
            self._parse_link(fields)
        elif record_type == 'P':  # Path
            self._parse_path(fields)
        elif record_type == 'W':  # Walk (alternative path format)
            self._parse_walk(fields)
    
    def _parse_segment(self, fields: List[str]) -> None:
        """Parse segment line: S <id> <sequence> [tags...]"""
        if len(fields) < 3:
            return
        
        node_id = fields[1]
        sequence = fields[2]
        
        self.segments[node_id] = sequence
        self.node_lengths[node_id] = len(sequence)
    
    def _parse_link(self, fields: List[str]) -> None:
        """Parse link line: L <from> <from_orient> <to> <to_orient> <overlap>"""
        if len(fields) < 6:
            return
        
        from_node = fields[1]
        from_orient = fields[2]
        to_node = fields[3]
        to_orient = fields[4]
        
        # Store bidirectional links
        self.links[from_node].append((to_node, from_orient + to_orient))
        self.links[to_node].append((from_node, to_orient + from_orient))
    
    def _parse_path(self, fields: List[str]) -> None:
        """Parse path line: P <path_name> <node_list> <cigar>"""
        if len(fields) < 3:
            return
        
        path_name = fields[1]
        node_string = fields[2]
        
        # Parse node list (format: node1+,node2-,node3+)
        nodes = []
        for node_spec in node_string.split(','):
            if node_spec.endswith('+') or node_spec.endswith('-'):
                node_id = node_spec[:-1]
                orientation = node_spec[-1]
                nodes.append((node_id, orientation))
        
        self.paths[path_name] = nodes
    
    def _parse_walk(self, fields: List[str]) -> None:
        """Parse walk line: W <sample> <haplotype> <seqid> <start> <end> <walk>"""
        if len(fields) < 7:
            return
        
        sample = fields[1]
        haplotype = fields[2]
        seqid = fields[3]
        walk = fields[6]
        
        # Construct path name
        path_name = f"{sample}#{haplotype}#{seqid}"
        
        # Parse walk (similar to path but different format)
        nodes = []
        for step in walk.split('>'):
            step = step.strip('<')
            if step:
                if step.endswith('+') or step.endswith('-'):
                    node_id = step[:-1]
                    orientation = step[-1]
                    nodes.append((node_id, orientation))
        
        self.paths[path_name] = nodes
    
    def get_connected_components(self) -> List[Set[str]]:
        """Find connected components in the graph."""
        visited = set()
        components = []
        
        for node in self.segments:
            if node not in visited:
                component = self._dfs_component(node, visited)
                components.append(component)
        
        return components
    
    def _dfs_component(self, start_node: str, visited: Set[str]) -> Set[str]:
        """DFS to find connected component starting from a node."""
        component = set()
        stack = [start_node]
        
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            
            visited.add(node)
            component.add(node)
            
            # Add neighbors
            for neighbor, _ in self.links.get(node, []):
                if neighbor not in visited:
                    stack.append(neighbor)
        
        return component
    
    def get_path_nodes(self, path_name: str) -> List[str]:
        """Get list of node IDs for a path."""
        if path_name not in self.paths:
            return []
        
        return [node_id for node_id, _ in self.paths[path_name]]
    
    def are_nodes_adjacent(self, node1: str, node2: str) -> bool:
        """Check if two nodes are adjacent in the graph."""
        if node1 not in self.links:
            return False
        
        for neighbor, _ in self.links[node1]:
            if neighbor == node2:
                return True
        return False


class GFF3Parser:
    """Parse GFF3 files to extract feature information."""
    
    def __init__(self, gff3_file: Path):
        """
        Initialize GFF3 parser.
        
        Args:
            gff3_file: Path to GFF3 file
        """
        self.gff3_file = Path(gff3_file)
        self.features = {}  # feature_id -> feature_info
        self.feature_hierarchy = {}  # child_id -> parent_id
        self.children_map = defaultdict(list)  # parent_id -> [child_ids]
        self._validate_input()
    
    def _validate_input(self) -> None:
        """Validate input file exists."""
        if not self.gff3_file.exists():
            raise FileNotFoundError(f"GFF3 file not found: {self.gff3_file}")
    
    def parse(self) -> None:
        """Parse the GFF3 file to extract feature information."""
        logger.info(f"Parsing GFF3 file: {self.gff3_file}")
        
        try:
            # Create an in-memory database
            db = gffutils.create_db(
                str(self.gff3_file),
                ':memory:',
                merge_strategy='create_unique',
                verbose=False,
                disable_infer_genes=True,
                disable_infer_transcripts=True
            )
            
            # Extract feature information
            for feature in db.all_features():
                feature_info = self._extract_feature_info(feature)
                self.features[feature_info['id']] = feature_info
                
                # Build hierarchy
                if feature_info['parent_ids']:
                    for parent_id in feature_info['parent_ids']:
                        self.feature_hierarchy[feature_info['id']] = parent_id
                        self.children_map[parent_id].append(feature_info['id'])
            
            logger.info(f"Parsed {len(self.features)} features from GFF3")
            
        except Exception as e:
            logger.error(f"Error parsing GFF3 file: {e}")
            # Fall back to manual parsing if gffutils fails
            self._manual_parse()
    
    def _extract_feature_info(self, feature) -> Dict[str, Any]:
        """Extract relevant information from a gffutils feature."""
        # Get attributes
        attributes = dict(feature.attributes)
        
        # Extract ID
        feature_id = attributes.get('ID', [feature.id])[0] if 'ID' in attributes else feature.id
        
        # Extract parent IDs
        parent_ids = attributes.get('Parent', [])
        
        # Extract name
        name = attributes.get('Name', attributes.get('gene_name', [feature_id]))[0]
        
        # Extract biotype/gene_biotype
        biotype = attributes.get('biotype', attributes.get('gene_biotype', ['unknown']))[0]
        
        return {
            'id': feature_id,
            'name': name,
            'type': feature.featuretype,
            'seqid': feature.seqid,
            'start': feature.start,
            'end': feature.end,
            'strand': feature.strand,
            'biotype': biotype,
            'parent_ids': parent_ids,
            'attributes': attributes
        }
    
    def _manual_parse(self) -> None:
        """Manual parsing fallback if gffutils fails."""
        logger.info("Attempting manual GFF3 parsing...")
        
        with open(self.gff3_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                try:
                    feature_info = self._parse_gff3_line(line)
                    if feature_info:
                        self.features[feature_info['id']] = feature_info
                        
                        # Build hierarchy
                        if feature_info['parent_ids']:
                            for parent_id in feature_info['parent_ids']:
                                self.feature_hierarchy[feature_info['id']] = parent_id
                                self.children_map[parent_id].append(feature_info['id'])
                                
                except Exception as e:
                    logger.warning(f"Error parsing GFF3 line {line_num}: {e}")
                    continue
        
        logger.info(f"Manually parsed {len(self.features)} features from GFF3")
    
    def _parse_gff3_line(self, line: str) -> Optional[Dict[str, Any]]:
        """Parse a single GFF3 line."""
        fields = line.split('\t')
        if len(fields) != 9:
            return None
        
        seqid, source, feature_type, start, end, score, strand, phase, attributes_str = fields
        
        # Parse attributes
        attributes = {}
        for attr in attributes_str.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attributes[key] = [value]  # Store as list for consistency
        
        # Extract ID
        feature_id = attributes.get('ID', [f"feature_{start}_{end}"])[0]
        
        # Extract parent IDs
        parent_ids = attributes.get('Parent', [])
        if isinstance(parent_ids, str):
            parent_ids = parent_ids.split(',')
        
        # Extract name
        name = attributes.get('Name', attributes.get('gene_name', [feature_id]))[0]
        
        # Extract biotype
        biotype = attributes.get('biotype', attributes.get('gene_biotype', ['unknown']))[0]
        
        return {
            'id': feature_id,
            'name': name,
            'type': feature_type.lower(),
            'seqid': seqid,
            'start': int(start),
            'end': int(end),
            'strand': strand,
            'biotype': biotype,
            'parent_ids': parent_ids,
            'attributes': attributes
        }
    
    def get_feature(self, feature_id: str) -> Optional[Dict[str, Any]]:
        """Get feature information by ID."""
        return self.features.get(feature_id)
    
    def get_feature_type(self, feature_id: str) -> str:
        """Get feature type by ID."""
        feature = self.features.get(feature_id)
        return feature['type'] if feature else 'unknown'
    
    def get_parent_feature(self, feature_id: str) -> Optional[Dict[str, Any]]:
        """Get parent feature information."""
        parent_id = self.feature_hierarchy.get(feature_id)
        return self.features.get(parent_id) if parent_id else None
    
    def get_children_features(self, feature_id: str) -> List[Dict[str, Any]]:
        """Get child features."""
        child_ids = self.children_map.get(feature_id, [])
        return [self.features[child_id] for child_id in child_ids if child_id in self.features]
    
    def get_gene_structure(self, feature_id: str) -> Dict[str, Any]:
        """Get complete gene structure for a feature."""
        feature = self.get_feature(feature_id)
        if not feature:
            return {}
        
        # Find the gene (top-level parent)
        gene = feature
        while gene and gene['type'] not in ['gene', 'pseudogene']:
            parent = self.get_parent_feature(gene['id'])
            if parent:
                gene = parent
            else:
                break
        
        # Get all descendants
        children = self._get_all_descendants(gene['id'])
        
        return {
            'gene': gene,
            'children': children,
            'transcripts': [c for c in children if c['type'] in ['mRNA', 'transcript', 'lnc_RNA', 'miRNA', 'tRNA', 'rRNA']],
            'exons': [c for c in children if c['type'] == 'exon'],
            'cds': [c for c in children if c['type'] == 'CDS'],
            'utrs': [c for c in children if c['type'] in ['five_prime_UTR', 'three_prime_UTR', 'UTR']]
        }
    
    def _get_all_descendants(self, feature_id: str) -> List[Dict[str, Any]]:
        """Get all descendant features recursively."""
        descendants = []
        
        def _collect_children(parent_id):
            children = self.get_children_features(parent_id)
            for child in children:
                descendants.append(child)
                _collect_children(child['id'])
        
        _collect_children(feature_id)
        return descendants


class ImpactDetector:
    """
    Detect impact of structural variations on genomic features.
    
    Analyzes GAM alignment data against GFA graph structure to determine
    how structural variations affect feature integrity.
    """
    
    # General impact types
    GENERAL_IMPACT_TYPES = {
        'INTACT': 'Feature is completely preserved',
        'TRUNCATED': 'Feature is partially missing at ends',
        'SPLIT': 'Feature is broken into multiple fragments',
        'MISSING': 'Feature is completely absent or unmappable'
    }
    
    # Structural impact types
    STRUCTURAL_IMPACT_TYPES = {
        'INVERSION': 'Feature sequence is reversed in orientation',
        'DUPLICATION': 'Feature appears multiple times in the path',
        'DELETION': 'Feature is missing from the expected path location',
        'TRANSLOCATION': 'Feature is split across non-adjacent graph nodes',
        'COPY_NUMBER_CHANGE': 'Feature copy count differs from reference',
        'COMPLEX_REARRANGEMENT': 'Combination of multiple structural variations'
    }
    
    # Feature-specific consequences
    CDS_CONSEQUENCES = {
        'FRAMESHIFT': 'Reading frame is disrupted by indel not divisible by 3',
        'START_LOST': 'Start codon is missing or disrupted',
        'STOP_LOST': 'Stop codon is missing or disrupted', 
        'INFRAME_INDEL': 'In-frame insertion or deletion (divisible by 3)',
        'SYNONYMOUS': 'No amino acid change expected',
        'NONSENSE': 'Premature stop codon introduced'
    }
    
    PROMOTER_CONSEQUENCES = {
        'COMPLETE': 'Promoter region is fully preserved',
        'PARTIAL': 'Promoter region is partially preserved',
        'MISSING': 'Promoter region is completely missing',
        'CORE_DISRUPTED': 'Core promoter elements are disrupted',
        'TSS_SHIFTED': 'Transcription start site position is altered'
    }
    
    SPLICE_SITE_CONSEQUENCES = {
        'EXACT': 'Splice site sequence exactly preserved',
        'SHIFTED': 'Splice site position shifted but potentially functional',
        'LOST': 'Splice site sequence completely lost',
        'WEAKENED': 'Splice site sequence altered, reducing efficiency',
        'CRYPTIC_CREATED': 'New cryptic splice site created nearby'
    }
    
    UTR_CONSEQUENCES = {
        'COMPLETE': 'UTR is fully preserved',
        'PARTIAL': 'UTR is partially preserved',
        'MISSING': 'UTR is completely missing',
        'REGULATORY_LOST': 'Regulatory elements within UTR are disrupted'
    }
    
    EXON_CONSEQUENCES = {
        'COMPLETE': 'Exon is fully preserved',
        'PARTIAL': 'Exon is partially preserved',
        'MISSING': 'Exon is completely missing',
        'SKIPPED': 'Exon appears to be skipped in splicing'
    }
    
    def __init__(self, gfa_file: Path, gff3_file: Optional[Path] = None,
                 min_alignment_coverage: float = 0.8,
                 min_identity_threshold: float = 0.9):
        """
        Initialize impact detector.
        
        Args:
            gfa_file: Path to GFA graph file
            gff3_file: Path to GFF3 annotation file (optional)
            min_alignment_coverage: Minimum coverage to consider feature intact
            min_identity_threshold: Minimum identity to consider alignment valid
        """
        self.gfa_file = Path(gfa_file)
        self.gff3_file = Path(gff3_file) if gff3_file else None
        self.min_alignment_coverage = min_alignment_coverage
        self.min_identity_threshold = min_identity_threshold
        
        self.gfa_parser = GFAParser(gfa_file)
        self.gff3_parser = GFF3Parser(gff3_file) if gff3_file else None
        self.connected_components = []
        self.node_to_component = {}
        
        self._parse_graph()
        if self.gff3_parser:
            self._parse_annotations()
    
    def _parse_graph(self) -> None:
        """Parse the GFA graph and identify connected components."""
        logger.info("Parsing GFA graph structure...")
        self.gfa_parser.parse()
        
        # Find connected components
        self.connected_components = self.gfa_parser.get_connected_components()
        
        # Create node to component mapping
        for i, component in enumerate(self.connected_components):
            for node in component:
                self.node_to_component[node] = i
        
        logger.info(f"Found {len(self.connected_components)} connected components")
    
    def _parse_annotations(self) -> None:
        """Parse GFF3 annotations."""
        logger.info("Parsing GFF3 annotations...")
        self.gff3_parser.parse()
    
    def analyze_impacts(self, gam_data: Dict[str, Dict[str, Dict[str, List[Dict[str, Any]]]]]) -> Dict[str, Dict[str, Dict[str, Dict[str, Any]]]]:
        """
        Analyze impacts for all features in GAM data.
        
        Args:
            gam_data: Parsed GAM data from GAMParser.group_alignments_by_sample_haplotype()
            
        Returns:
            Impact analysis: {sample: {haplotype: {feature: {'type': type, 'consequence': consequence, 'structural_impacts': [...], 'details': {...}}}}}
        """
        logger.info("Analyzing feature impacts...")
        
        impact_results = {}
        
        for sample in gam_data:
            impact_results[sample] = {}
            
            for haplotype in gam_data[sample]:
                impact_results[sample][haplotype] = {}
                
                for feature_type in gam_data[sample][haplotype]:
                    features = gam_data[sample][haplotype][feature_type]
                    
                    for feature in features:
                        feature_name = feature['read_name']
                        impact = self._analyze_single_feature_impact(feature, sample, haplotype)
                        impact_results[sample][haplotype][feature_name] = impact
        
        return impact_results
    
    def _analyze_single_feature_impact(self, feature: Dict[str, Any], 
                                     sample: str, haplotype: str) -> Dict[str, Any]:
        """
        Analyze impact for a single feature with feature-type specific analysis.
        
        Args:
            feature: Feature alignment data
            sample: Sample name
            haplotype: Haplotype identifier
            
        Returns:
            Impact analysis for the feature with type, consequence, and structural impacts
        """
        feature_name = feature['read_name']
        
        # Get feature type from GFF3 if available, otherwise from GAM data
        gff3_feature = None
        if self.gff3_parser:
            gff3_feature = self.gff3_parser.get_feature(feature_name)
            feature_type = gff3_feature['type'] if gff3_feature else feature.get('feature_type', 'unknown')
        else:
            feature_type = feature.get('feature_type', 'unknown')
        
        sequence_length = len(feature['sequence'])
        identity = feature.get('identity', 0.0)
        score = feature.get('score', 0)
        
        # Initialize impact analysis with new structure
        impact_analysis = {
            'type': 'MISSING',
            'consequence': 'MISSING',
            'structural_impacts': [],
            'details': {
                'feature_type': feature_type,
                'sequence_length': sequence_length,
                'identity': identity,
                'score': score,
                'coverage': 0.0,
                'components_spanned': [],
                'is_split': False,
                'is_truncated': False,
                'fragments': [],
                'structural_variation_detected': False,
                'orientation_changes': [],
                'path_discontinuities': [],
                'copy_number': 0,
                'expected_copy_number': 1,
                'is_inverted': False,
                'translocation_detected': False,
                'gff3_info': gff3_feature
            }
        }
        
        # Check if we have path positions
        path_positions = feature.get('path_positions', [])
        if not path_positions:
            impact_analysis['details']['reason'] = 'No alignment positions found'
            return impact_analysis
        
        # Check identity threshold
        if identity < self.min_identity_threshold:
            impact_analysis['details']['reason'] = f'Identity {identity:.3f} below threshold {self.min_identity_threshold}'
            return impact_analysis
        
        # Analyze alignment positions
        aligned_nodes = self._extract_aligned_nodes(path_positions)
        if not aligned_nodes:
            impact_analysis['details']['reason'] = 'No valid aligned nodes found'
            return impact_analysis
        
        # Calculate coverage
        coverage = self._calculate_coverage(aligned_nodes, sequence_length)
        impact_analysis['details']['coverage'] = coverage
        
        # Determine which components the feature spans
        components_spanned = self._get_components_spanned(aligned_nodes)
        impact_analysis['details']['components_spanned'] = components_spanned
        
        # Check for structural variation (spans multiple components)
        if len(components_spanned) > 1:
            impact_analysis['details']['structural_variation_detected'] = True
        
        # Analyze fragments and gaps
        fragments = self._analyze_fragments(aligned_nodes, sequence_length)
        impact_analysis['details']['fragments'] = fragments
        
        # Analyze structural impacts
        structural_impacts = self._detect_structural_impacts(feature, aligned_nodes, path_positions)
        impact_analysis['structural_impacts'] = structural_impacts
        
        # Update details with structural analysis
        impact_analysis['details'].update(self._analyze_structural_details(
            aligned_nodes, path_positions, fragments
        ))
        
        # Determine general impact type
        general_impact = self._determine_general_impact_type(coverage, fragments, components_spanned)
        impact_analysis['type'] = general_impact
        
        # Set flags based on impact type
        if general_impact == 'SPLIT':
            impact_analysis['details']['is_split'] = True
        elif general_impact == 'TRUNCATED':
            impact_analysis['details']['is_truncated'] = True
        
        # Dispatch to feature-specific analysis using GFF3 information
        consequence = self._analyze_feature_specific_impact(
            feature, general_impact, coverage, fragments, components_spanned, gff3_feature
        )
        impact_analysis['consequence'] = consequence
        
        return impact_analysis
    
    def _detect_structural_impacts(self, feature: Dict[str, Any], 
                                 aligned_nodes: List[str], 
                                 path_positions: List[Dict[str, Any]]) -> List[str]:
        """
        Detect structural impacts affecting the feature.
        
        Args:
            feature: Feature data
            aligned_nodes: List of aligned node IDs
            path_positions: Path position information
            
        Returns:
            List of detected structural impact types
        """
        structural_impacts = []
        
        # Check for inversion
        if self._detect_inversion(path_positions):
            structural_impacts.append('INVERSION')
        
        # Check for duplication
        if self._detect_duplication(aligned_nodes):
            structural_impacts.append('DUPLICATION')
        
        # Check for deletion
        if self._detect_deletion(feature, aligned_nodes):
            structural_impacts.append('DELETION')
        
        # Check for translocation
        if self._detect_translocation(aligned_nodes):
            structural_impacts.append('TRANSLOCATION')
        
        # Check for copy number change
        if self._detect_copy_number_change(feature, aligned_nodes):
            structural_impacts.append('COPY_NUMBER_CHANGE')
        
        # Check for complex rearrangement
        if len(structural_impacts) >= 2:
            structural_impacts.append('COMPLEX_REARRANGEMENT')
        
        return structural_impacts
    
    def _detect_inversion(self, path_positions: List[Dict[str, Any]]) -> bool:
        """Detect if feature sequence is inverted."""
        if len(path_positions) < 2:
            return False
        
        # Count reverse orientations
        reverse_count = sum(1 for pos in path_positions 
                          if pos.get('is_reverse', False))
        total_positions = len(path_positions)
        
        # If majority of positions are reverse, it's likely inverted
        return reverse_count > total_positions * 0.6
    
    def _detect_duplication(self, aligned_nodes: List[str]) -> bool:
        """Detect if feature appears multiple times (duplication)."""
        if not aligned_nodes:
            return False
        
        # Count node occurrences
        node_counts = Counter(aligned_nodes)
        
        # If any node appears more than twice, likely duplication
        return any(count > 2 for count in node_counts.values())
    
    def _detect_deletion(self, feature: Dict[str, Any], aligned_nodes: List[str]) -> bool:
        """Detect if feature is deleted (missing expected alignment)."""
        sequence_length = len(feature.get('sequence', ''))
        
        # If we have very few aligned nodes for a long sequence, it's likely deleted
        if sequence_length > 100 and len(aligned_nodes) < 3:
            return True
        
        # If coverage is very low, it might be deleted
        coverage = self._calculate_coverage(aligned_nodes, sequence_length)
        return coverage < 0.2
    
    def _detect_translocation(self, aligned_nodes: List[str]) -> bool:
        """Detect if feature is split across non-adjacent nodes."""
        if len(aligned_nodes) < 2:
            return False
        
        # Check if consecutive nodes in alignment are adjacent in graph
        non_adjacent_count = 0
        for i in range(len(aligned_nodes) - 1):
            node1 = aligned_nodes[i]
            node2 = aligned_nodes[i + 1]
            
            if not self.gfa_parser.are_nodes_adjacent(node1, node2):
                non_adjacent_count += 1
        
        # If more than 30% of consecutive pairs are non-adjacent, it's translocation
        return non_adjacent_count > len(aligned_nodes) * 0.3
    
    def _detect_copy_number_change(self, feature: Dict[str, Any], 
                                 aligned_nodes: List[str]) -> bool:
        """Detect if feature copy number differs from expected."""
        if not aligned_nodes:
            return False
        
        # Count unique nodes vs duplicated nodes
        node_counts = Counter(aligned_nodes)
        unique_nodes = len(node_counts)
        total_occurrences = sum(node_counts.values())
        
        # If we have significantly more occurrences than unique nodes, 
        # it suggests duplication/copy number change
        duplication_ratio = total_occurrences / unique_nodes if unique_nodes > 0 else 1
        
        # Only consider it copy number change if ratio is significantly > 1
        # and we have enough nodes to make a meaningful assessment
        return duplication_ratio > 1.8 and unique_nodes >= 2
    
    def _analyze_structural_details(self, aligned_nodes: List[str], 
                                  path_positions: List[Dict[str, Any]],
                                  fragments: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Analyze detailed structural information.
        
        Args:
            aligned_nodes: List of aligned node IDs
            path_positions: Path position information
            fragments: Fragment analysis
            
        Returns:
            Dictionary with structural analysis details
        """
        details = {}
        
        # Analyze orientation changes
        details['orientation_changes'] = self._analyze_orientation_changes(path_positions)
        
        # Analyze path discontinuities
        details['path_discontinuities'] = self._analyze_path_discontinuities(aligned_nodes)
        
        # Calculate copy number
        details['copy_number'] = self._calculate_copy_number(aligned_nodes)
        details['expected_copy_number'] = 1  # Default assumption
        
        # Check if inverted
        details['is_inverted'] = self._detect_inversion(path_positions)
        
        # Check for translocation
        details['translocation_detected'] = self._detect_translocation(aligned_nodes)
        
        return details
    
    def _analyze_orientation_changes(self, path_positions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Analyze orientation changes in the alignment."""
        changes = []
        
        if len(path_positions) < 2:
            return changes
        
        current_orientation = path_positions[0].get('is_reverse', False)
        
        for i, pos in enumerate(path_positions[1:], 1):
            pos_orientation = pos.get('is_reverse', False)
            
            if pos_orientation != current_orientation:
                changes.append({
                    'position': i,
                    'from_orientation': '+' if not current_orientation else '-',
                    'to_orientation': '+' if not pos_orientation else '-',
                    'node_id': pos.get('node_id', 'unknown')
                })
                current_orientation = pos_orientation
        
        return changes
    
    def _analyze_path_discontinuities(self, aligned_nodes: List[str]) -> List[Dict[str, Any]]:
        """Analyze discontinuities in the alignment path."""
        discontinuities = []
        
        if len(aligned_nodes) < 2:
            return discontinuities
        
        for i in range(len(aligned_nodes) - 1):
            node1 = aligned_nodes[i]
            node2 = aligned_nodes[i + 1]
            
            if not self.gfa_parser.are_nodes_adjacent(node1, node2):
                # Calculate component difference
                comp1 = self.node_to_component.get(node1, -1)
                comp2 = self.node_to_component.get(node2, -1)
                
                discontinuities.append({
                    'position': i,
                    'from_node': node1,
                    'to_node': node2,
                    'from_component': comp1,
                    'to_component': comp2,
                    'cross_component': comp1 != comp2
                })
        
        return discontinuities
    
    def _calculate_copy_number(self, aligned_nodes: List[str]) -> int:
        """Calculate estimated copy number based on alignment patterns."""
        if not aligned_nodes:
            return 0
        
        # Count unique nodes
        unique_nodes = len(set(aligned_nodes))
        
        # For simple cases, copy number is just the number of unique nodes
        # In more complex cases, this would need more sophisticated analysis
        return max(1, unique_nodes)
    
    def _analyze_feature_specific_impact(self, feature: Dict[str, Any], 
                                       general_impact: str, coverage: float,
                                       fragments: List[Dict[str, Any]], 
                                       components_spanned: List[int],
                                       gff3_feature: Optional[Dict[str, Any]] = None) -> str:
        """
        Analyze feature-specific consequences based on feature type.
        
        Args:
            feature: Feature data
            general_impact: General impact type (INTACT, TRUNCATED, SPLIT, MISSING)
            coverage: Alignment coverage
            fragments: Fragment analysis
            components_spanned: Components spanned by alignment
            gff3_feature: GFF3 feature information if available
            
        Returns:
            Feature-specific consequence
        """
        # Use GFF3 feature type if available, otherwise fall back to GAM data
        if gff3_feature:
            feature_type = gff3_feature['type'].lower()
            # Also consider parent feature types for context
            if self.gff3_parser:
                parent = self.gff3_parser.get_parent_feature(gff3_feature['id'])
                parent_type = parent['type'].lower() if parent else None
        else:
            feature_type = feature.get('feature_type', '').lower()
            parent_type = None
        
        # Dispatch to appropriate analyzer
        if feature_type == 'cds':
            return self._analyze_cds_impact(feature, general_impact, coverage, fragments, gff3_feature)
        elif feature_type in ['promoter', 'regulatory']:
            return self._analyze_promoter_impact(feature, general_impact, coverage, fragments, gff3_feature)
        elif 'splice' in feature_type or feature_type == 'splice_site':
            return self._analyze_splice_site_impact(feature, general_impact, coverage, fragments, gff3_feature)
        elif 'utr' in feature_type or feature_type in ['five_prime_utr', 'three_prime_utr', 'utr']:
            return self._analyze_utr_impact(feature, general_impact, coverage, fragments, gff3_feature)
        elif feature_type == 'exon':
            return self._analyze_exon_impact(feature, general_impact, coverage, fragments, gff3_feature)
        elif feature_type in ['gene', 'pseudogene']:
            return self._analyze_gene_impact(feature, general_impact, coverage, fragments, gff3_feature)
        elif feature_type in ['mrna', 'transcript', 'lnc_rna', 'mirna', 'trna', 'rrna']:
            return self._analyze_transcript_impact(feature, general_impact, coverage, fragments, gff3_feature)
        else:
            # Default to general impact for unknown feature types
            return general_impact
    
    def _analyze_cds_impact(self, feature: Dict[str, Any], general_impact: str,
                          coverage: float, fragments: List[Dict[str, Any]],
                          gff3_feature: Optional[Dict[str, Any]] = None) -> str:
        """Analyze CDS-specific impact."""
        sequence = feature.get('sequence', '')
        sequence_length = len(sequence)
        
        if general_impact == 'MISSING':
            return 'MISSING'
        
        # Check for frameshift based on fragment analysis
        if len(fragments) > 1:
            # Multiple fragments suggest potential frameshift
            total_gap_length = 0
            for i in range(len(fragments) - 1):
                # Estimate gap between fragments (simplified)
                gap_estimate = abs(len(fragments[i]['nodes']) - len(fragments[i+1]['nodes']))
                total_gap_length += gap_estimate
            
            if total_gap_length % 3 != 0:
                return 'FRAMESHIFT'
            else:
                return 'INFRAME_INDEL'
        
        # Check for start/stop codon loss based on coverage at ends
        if coverage < self.min_alignment_coverage:
            if coverage < 0.5:
                # Low coverage might indicate start or stop loss
                # This is simplified - in reality we'd need position information
                if len(sequence) >= 3:
                    # Check if likely start codon region affected
                    if 'ATG' in sequence[:9]:  # Start region
                        return 'START_LOST'
                    # Check if likely stop codon region affected  
                    elif any(stop in sequence[-9:] for stop in ['TAA', 'TAG', 'TGA']):
                        return 'STOP_LOST'
            
            return 'TRUNCATED'
        
        # High coverage, likely intact or synonymous
        if coverage >= 0.95:
            return 'SYNONYMOUS'
        else:
            return 'INFRAME_INDEL'
    
    def _analyze_promoter_impact(self, feature: Dict[str, Any], general_impact: str,
                               coverage: float, fragments: List[Dict[str, Any]],
                               gff3_feature: Optional[Dict[str, Any]] = None) -> str:
        """Analyze promoter/regulatory element impact."""
        if general_impact == 'MISSING':
            return 'MISSING'
        
        if coverage >= 0.9:
            return 'COMPLETE'
        elif coverage >= 0.5:
            # Check if core promoter region likely affected
            if len(fragments) > 1:
                return 'CORE_DISRUPTED'
            else:
                return 'PARTIAL'
        elif coverage >= 0.2:
            return 'TSS_SHIFTED'
        else:
            return 'MISSING'
    
    def _analyze_splice_site_impact(self, feature: Dict[str, Any], general_impact: str,
                                  coverage: float, fragments: List[Dict[str, Any]],
                                  gff3_feature: Optional[Dict[str, Any]] = None) -> str:
        """Analyze splice site impact."""
        if general_impact == 'MISSING':
            return 'LOST'
        
        sequence = feature.get('sequence', '')
        
        # Check for canonical splice site sequences
        has_gt_ag = 'GT' in sequence and 'AG' in sequence
        has_gc_ag = 'GC' in sequence and 'AG' in sequence
        
        if coverage >= 0.95:
            if has_gt_ag or has_gc_ag:
                return 'EXACT'
            else:
                return 'WEAKENED'
        elif coverage >= 0.7:
            return 'SHIFTED'
        elif coverage >= 0.3:
            if has_gt_ag or has_gc_ag:
                return 'WEAKENED'
            else:
                return 'CRYPTIC_CREATED'
        else:
            return 'LOST'
    
    def _analyze_utr_impact(self, feature: Dict[str, Any], general_impact: str,
                          coverage: float, fragments: List[Dict[str, Any]],
                          gff3_feature: Optional[Dict[str, Any]] = None) -> str:
        """Analyze UTR impact."""
        if general_impact == 'MISSING':
            return 'MISSING'
        
        if coverage >= 0.9:
            return 'COMPLETE'
        elif coverage >= 0.5:
            # Check if regulatory elements likely affected
            if len(fragments) > 1:
                return 'REGULATORY_LOST'
            else:
                return 'PARTIAL'
        else:
            return 'MISSING'
    
    def _analyze_exon_impact(self, feature: Dict[str, Any], general_impact: str,
                           coverage: float, fragments: List[Dict[str, Any]],
                           gff3_feature: Optional[Dict[str, Any]] = None) -> str:
        """Analyze exon impact."""
        if general_impact == 'MISSING':
            return 'MISSING'
        
        if coverage >= 0.9:
            return 'COMPLETE'
        elif coverage >= 0.3:
            if len(fragments) > 1:
                return 'SKIPPED'
            else:
                return 'PARTIAL'
        else:
            return 'MISSING'
    
    def _analyze_gene_impact(self, feature: Dict[str, Any], general_impact: str,
                           coverage: float, fragments: List[Dict[str, Any]],
                           gff3_feature: Optional[Dict[str, Any]] = None) -> str:
        """Analyze gene-level impact."""
        if general_impact == 'MISSING':
            return 'MISSING'
        
        if coverage >= 0.9:
            return 'COMPLETE'
        elif coverage >= 0.5:
            if len(fragments) > 1:
                return 'DISRUPTED'
            else:
                return 'PARTIAL'
        else:
            return 'MISSING'
    
    def _analyze_transcript_impact(self, feature: Dict[str, Any], general_impact: str,
                                 coverage: float, fragments: List[Dict[str, Any]],
                                 gff3_feature: Optional[Dict[str, Any]] = None) -> str:
        """Analyze transcript-level impact."""
        if general_impact == 'MISSING':
            return 'MISSING'
        
        if coverage >= 0.9:
            return 'COMPLETE'
        elif coverage >= 0.5:
            if len(fragments) > 1:
                return 'SPLICING_DISRUPTED'
            else:
                return 'PARTIAL'
        else:
            return 'MISSING'
    
    def _extract_aligned_nodes(self, path_positions: List[Dict[str, Any]]) -> List[str]:
        """Extract list of aligned node IDs from path positions."""
        nodes = []
        for pos in path_positions:
            node_id = str(pos.get('node_id', ''))
            if node_id and node_id in self.gfa_parser.segments:
                nodes.append(node_id)
        return nodes
    
    def _calculate_coverage(self, aligned_nodes: List[str], sequence_length: int) -> float:
        """Calculate approximate coverage based on aligned nodes."""
        if not aligned_nodes or sequence_length == 0:
            return 0.0
        
        # Sum lengths of aligned nodes (approximate)
        total_aligned_length = 0
        for node_id in aligned_nodes:
            total_aligned_length += self.gfa_parser.node_lengths.get(node_id, 0)
        
        # Coverage is aligned length divided by sequence length
        # This is approximate since we don't have exact alignment boundaries
        coverage = min(1.0, total_aligned_length / sequence_length)
        return coverage
    
    def _get_components_spanned(self, aligned_nodes: List[str]) -> List[int]:
        """Get list of connected components spanned by the alignment."""
        components = set()
        for node_id in aligned_nodes:
            component_id = self.node_to_component.get(node_id)
            if component_id is not None:
                components.add(component_id)
        return sorted(list(components))
    
    def _analyze_fragments(self, aligned_nodes: List[str], sequence_length: int) -> List[Dict[str, Any]]:
        """Analyze alignment fragments and gaps."""
        if not aligned_nodes:
            return []
        
        fragments = []
        current_fragment = {
            'nodes': [],
            'start_node': aligned_nodes[0],
            'end_node': aligned_nodes[0],
            'component': self.node_to_component.get(aligned_nodes[0])
        }
        
        for i, node_id in enumerate(aligned_nodes):
            node_component = self.node_to_component.get(node_id)
            
            # If we're in the same component as current fragment, extend it
            if node_component == current_fragment['component']:
                current_fragment['nodes'].append(node_id)
                current_fragment['end_node'] = node_id
            else:
                # New component - finish current fragment and start new one
                if current_fragment['nodes']:
                    fragments.append(current_fragment.copy())
                
                current_fragment = {
                    'nodes': [node_id],
                    'start_node': node_id,
                    'end_node': node_id,
                    'component': node_component
                }
        
        # Add the last fragment
        if current_fragment['nodes']:
            fragments.append(current_fragment)
        
        return fragments
    
    def _determine_general_impact_type(self, coverage: float, fragments: List[Dict[str, Any]], 
                                     components_spanned: List[int]) -> str:
        """Determine the general impact type based on analysis results."""
        # Missing: very low coverage or no fragments
        if coverage < 0.1 or not fragments:
            return 'MISSING'
        
        # Split: multiple fragments in different components
        if len(fragments) > 1 and len(components_spanned) > 1:
            return 'SPLIT'
        
        # Intact: high coverage in single component
        if coverage >= self.min_alignment_coverage and len(components_spanned) == 1:
            return 'INTACT'
        
        # Truncated: moderate coverage but not complete
        if coverage >= 0.3:
            return 'TRUNCATED'
        
        # Default to missing
        return 'MISSING'
    
    def get_impact_summary(self, impact_results: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]]) -> Dict[str, Any]:
        """
        Generate summary statistics for impact analysis.
        
        Args:
            impact_results: Results from analyze_impacts()
            
        Returns:
            Summary statistics
        """
        summary = {
            'total_features': 0,
            'impact_type_counts': {impact_type: 0 for impact_type in self.GENERAL_IMPACT_TYPES},
            'structural_impact_counts': {impact_type: 0 for impact_type in self.STRUCTURAL_IMPACT_TYPES},
            'consequence_counts': defaultdict(int),
            'feature_type_analysis': defaultdict(lambda: {
                'total': 0,
                'impact_types': {impact_type: 0 for impact_type in self.GENERAL_IMPACT_TYPES},
                'structural_impacts': {impact_type: 0 for impact_type in self.STRUCTURAL_IMPACT_TYPES},
                'consequences': defaultdict(int)
            }),
            'structural_variations_detected': 0,
            'features_spanning_multiple_components': 0,
            'orientation_changes_detected': 0,
            'path_discontinuities_detected': 0,
            'avg_coverage': 0.0,
            'avg_identity': 0.0,
            'samples_analyzed': len(impact_results),
            'sample_summaries': {},
            'gff3_annotations_used': self.gff3_parser is not None
        }
        
        total_coverage = 0.0
        total_identity = 0.0
        
        for sample in impact_results:
            sample_summary = {
                'haplotypes': len(impact_results[sample]),
                'total_features': 0,
                'impact_type_counts': {impact_type: 0 for impact_type in self.GENERAL_IMPACT_TYPES},
                'structural_impact_counts': {impact_type: 0 for impact_type in self.STRUCTURAL_IMPACT_TYPES},
                'consequence_counts': defaultdict(int)
            }
            
            for haplotype in impact_results[sample]:
                for feature_name, analysis in impact_results[sample][haplotype].items():
                    summary['total_features'] += 1
                    sample_summary['total_features'] += 1
                    
                    impact_type = analysis['type']
                    consequence = analysis['consequence']
                    structural_impacts = analysis.get('structural_impacts', [])
                    feature_type = analysis['details'].get('feature_type', 'unknown')
                    
                    # Count impact types
                    summary['impact_type_counts'][impact_type] += 1
                    sample_summary['impact_type_counts'][impact_type] += 1
                    
                    # Count structural impacts
                    for struct_impact in structural_impacts:
                        if struct_impact in self.STRUCTURAL_IMPACT_TYPES:
                            summary['structural_impact_counts'][struct_impact] += 1
                            sample_summary['structural_impact_counts'][struct_impact] += 1
                    
                    # Count consequences
                    summary['consequence_counts'][consequence] += 1
                    sample_summary['consequence_counts'][consequence] += 1
                    
                    # Feature type specific analysis
                    summary['feature_type_analysis'][feature_type]['total'] += 1
                    summary['feature_type_analysis'][feature_type]['impact_types'][impact_type] += 1
                    summary['feature_type_analysis'][feature_type]['consequences'][consequence] += 1
                    
                    for struct_impact in structural_impacts:
                        if struct_impact in self.STRUCTURAL_IMPACT_TYPES:
                            summary['feature_type_analysis'][feature_type]['structural_impacts'][struct_impact] += 1
                    
                    details = analysis['details']
                    total_coverage += details.get('coverage', 0.0)
                    total_identity += details.get('identity', 0.0)
                    
                    if details.get('structural_variation_detected', False):
                        summary['structural_variations_detected'] += 1
                    
                    if len(details.get('components_spanned', [])) > 1:
                        summary['features_spanning_multiple_components'] += 1
                    
                    if details.get('orientation_changes', []):
                        summary['orientation_changes_detected'] += 1
                    
                    if details.get('path_discontinuities', []):
                        summary['path_discontinuities_detected'] += 1
            
            # Convert defaultdicts to regular dicts for JSON serialization
            sample_summary['consequence_counts'] = dict(sample_summary['consequence_counts'])
            summary['sample_summaries'][sample] = sample_summary
        
        # Calculate averages
        if summary['total_features'] > 0:
            summary['avg_coverage'] = total_coverage / summary['total_features']
            summary['avg_identity'] = total_identity / summary['total_features']
        
        # Convert defaultdicts to regular dicts for JSON serialization
        summary['consequence_counts'] = dict(summary['consequence_counts'])
        summary['feature_type_analysis'] = {
            ft: {
                'total': data['total'],
                'impact_types': data['impact_types'],
                'structural_impacts': data['structural_impacts'],
                'consequences': dict(data['consequences'])
            }
            for ft, data in summary['feature_type_analysis'].items()
        }
        
        return summary
    
    def save_results(self, impact_results: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]], 
                    output_file: Path) -> None:
        """
        Save impact analysis results to JSON file.
        
        Args:
            impact_results: Results from analyze_impacts()
            output_file: Output file path
        """
        output_data = {
            'impact_analysis': impact_results,
            'summary': self.get_impact_summary(impact_results),
            'parameters': {
                'min_alignment_coverage': self.min_alignment_coverage,
                'min_identity_threshold': self.min_identity_threshold,
                'gfa_file': str(self.gfa_file),
                'gff3_file': str(self.gff3_file) if self.gff3_file else None
            },
            'impact_type_descriptions': {
                'general_types': self.GENERAL_IMPACT_TYPES,
                'structural_types': self.STRUCTURAL_IMPACT_TYPES,
                'cds_consequences': self.CDS_CONSEQUENCES,
                'promoter_consequences': self.PROMOTER_CONSEQUENCES,
                'splice_site_consequences': self.SPLICE_SITE_CONSEQUENCES,
                'utr_consequences': self.UTR_CONSEQUENCES,
                'exon_consequences': self.EXON_CONSEQUENCES
            }
        }
        
        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        logger.info(f"Saved impact analysis results to {output_file}")


def main():
    """Command-line interface for impact detector."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze feature impacts using GAM and GFA data")
    parser.add_argument("gfa_file", help="Input GFA graph file")
    parser.add_argument("gam_json", help="Input GAM data (JSON from GAMParser)")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file for impact analysis")
    parser.add_argument("--gff3", help="Input GFF3 annotation file for feature type information")
    parser.add_argument("--min-coverage", type=float, default=0.8, 
                       help="Minimum alignment coverage for intact features")
    parser.add_argument("--min-identity", type=float, default=0.9,
                       help="Minimum identity threshold for valid alignments")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Load GAM data
    logger.info(f"Loading GAM data from {args.gam_json}")
    with open(args.gam_json, 'r') as f:
        gam_data = json.load(f)
    
    # Initialize impact detector
    detector = ImpactDetector(
        Path(args.gfa_file),
        gff3_file=Path(args.gff3) if args.gff3 else None,
        min_alignment_coverage=args.min_coverage,
        min_identity_threshold=args.min_identity
    )
    
    # Analyze impacts
    impact_results = detector.analyze_impacts(gam_data)
    
    # Save results
    detector.save_results(impact_results, Path(args.output))
    
    # Print summary
    summary = detector.get_impact_summary(impact_results)
    print(f"\nImpact Analysis Summary:")
    print(f"Total features analyzed: {summary['total_features']}")
    if summary['gff3_annotations_used']:
        print("✓ Using GFF3 annotations for feature type information")
    else:
        print("ℹ No GFF3 annotations provided - using feature types from GAM data")
    
    print(f"\nImpact type distribution:")
    for impact_type, count in summary['impact_type_counts'].items():
        percentage = (count / summary['total_features'] * 100) if summary['total_features'] > 0 else 0
        print(f"  {impact_type}: {count} ({percentage:.1f}%)")
    
    print(f"\nStructural impact distribution:")
    for impact_type, count in summary['structural_impact_counts'].items():
        percentage = (count / summary['total_features'] * 100) if summary['total_features'] > 0 else 0
        print(f"  {impact_type}: {count} ({percentage:.1f}%)")
    
    print(f"\nTop consequences:")
    sorted_consequences = sorted(summary['consequence_counts'].items(), 
                               key=lambda x: x[1], reverse=True)[:10]
    for consequence, count in sorted_consequences:
        percentage = (count / summary['total_features'] * 100) if summary['total_features'] > 0 else 0
        print(f"  {consequence}: {count} ({percentage:.1f}%)")
    
    print(f"\nFeature type analysis:")
    for feature_type, type_data in summary['feature_type_analysis'].items():
        if type_data['total'] > 0:
            print(f"  {feature_type}: {type_data['total']} features")
            top_consequences = sorted(type_data['consequences'].items(), 
                                    key=lambda x: x[1], reverse=True)[:3]
            for consequence, count in top_consequences:
                print(f"    {consequence}: {count}")
    
    print(f"\nStructural variations detected: {summary['structural_variations_detected']}")
    print(f"Features spanning multiple components: {summary['features_spanning_multiple_components']}")
    print(f"Orientation changes detected: {summary['orientation_changes_detected']}")
    print(f"Path discontinuities detected: {summary['path_discontinuities_detected']}")
    print(f"Average coverage: {summary['avg_coverage']:.3f}")
    print(f"Average identity: {summary['avg_identity']:.3f}")


if __name__ == "__main__":
    main()
