#!/usr/bin/env python3
"""
Impact detector for analyzing feature alignments against pangenome graphs.

This module provides functionality to analyze GAM alignment data and determine
the impact of structural variations on genomic features.
"""

import logging
import re
from collections import defaultdict, deque
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any, Union
import json

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
    
    def __init__(self, gfa_file: Path, min_alignment_coverage: float = 0.8,
                 min_identity_threshold: float = 0.9):
        """
        Initialize impact detector.
        
        Args:
            gfa_file: Path to GFA graph file
            min_alignment_coverage: Minimum coverage to consider feature intact
            min_identity_threshold: Minimum identity to consider alignment valid
        """
        self.gfa_file = Path(gfa_file)
        self.min_alignment_coverage = min_alignment_coverage
        self.min_identity_threshold = min_identity_threshold
        
        self.gfa_parser = GFAParser(gfa_file)
        self.connected_components = []
        self.node_to_component = {}
        
        self._parse_graph()
    
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
    
    def analyze_impacts(self, gam_data: Dict[str, Dict[str, Dict[str, List[Dict[str, Any]]]]]) -> Dict[str, Dict[str, Dict[str, Dict[str, Any]]]]:
        """
        Analyze impacts for all features in GAM data.
        
        Args:
            gam_data: Parsed GAM data from GAMParser.group_alignments_by_sample_haplotype()
            
        Returns:
            Impact analysis: {sample: {haplotype: {feature: {'type': type, 'consequence': consequence, 'details': {...}}}}}
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
            Impact analysis for the feature with type and consequence
        """
        feature_name = feature['read_name']
        feature_type = feature.get('feature_type', 'unknown')
        sequence_length = len(feature['sequence'])
        identity = feature.get('identity', 0.0)
        score = feature.get('score', 0)
        
        # Initialize impact analysis with new structure
        impact_analysis = {
            'type': 'MISSING',
            'consequence': 'MISSING',
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
                'structural_variation_detected': False
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
        
        # Determine general impact type
        general_impact = self._determine_general_impact_type(coverage, fragments, components_spanned)
        impact_analysis['type'] = general_impact
        
        # Set flags based on impact type
        if general_impact == 'SPLIT':
            impact_analysis['details']['is_split'] = True
        elif general_impact == 'TRUNCATED':
            impact_analysis['details']['is_truncated'] = True
        
        # Dispatch to feature-specific analysis
        consequence = self._analyze_feature_specific_impact(
            feature, general_impact, coverage, fragments, components_spanned
        )
        impact_analysis['consequence'] = consequence
        
        return impact_analysis
    
    def _analyze_feature_specific_impact(self, feature: Dict[str, Any], 
                                       general_impact: str, coverage: float,
                                       fragments: List[Dict[str, Any]], 
                                       components_spanned: List[int]) -> str:
        """
        Analyze feature-specific consequences based on feature type.
        
        Args:
            feature: Feature data
            general_impact: General impact type (INTACT, TRUNCATED, SPLIT, MISSING)
            coverage: Alignment coverage
            fragments: Fragment analysis
            components_spanned: Components spanned by alignment
            
        Returns:
            Feature-specific consequence
        """
        feature_type = feature.get('feature_type', '').lower()
        
        # Dispatch to appropriate analyzer
        if feature_type == 'cds':
            return self._analyze_cds_impact(feature, general_impact, coverage, fragments)
        elif feature_type in ['promoter', 'regulatory']:
            return self._analyze_promoter_impact(feature, general_impact, coverage, fragments)
        elif 'splice' in feature_type or feature_type == 'splice_site':
            return self._analyze_splice_site_impact(feature, general_impact, coverage, fragments)
        elif 'utr' in feature_type or feature_type in ['5_prime_utr', '3_prime_utr']:
            return self._analyze_utr_impact(feature, general_impact, coverage, fragments)
        elif feature_type == 'exon':
            return self._analyze_exon_impact(feature, general_impact, coverage, fragments)
        else:
            # Default to general impact for unknown feature types
            return general_impact
    
    def _analyze_cds_impact(self, feature: Dict[str, Any], general_impact: str,
                          coverage: float, fragments: List[Dict[str, Any]]) -> str:
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
                               coverage: float, fragments: List[Dict[str, Any]]) -> str:
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
                                  coverage: float, fragments: List[Dict[str, Any]]) -> str:
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
                          coverage: float, fragments: List[Dict[str, Any]]) -> str:
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
                           coverage: float, fragments: List[Dict[str, Any]]) -> str:
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
            'consequence_counts': defaultdict(int),
            'feature_type_analysis': defaultdict(lambda: {
                'total': 0,
                'impact_types': {impact_type: 0 for impact_type in self.GENERAL_IMPACT_TYPES},
                'consequences': defaultdict(int)
            }),
            'structural_variations_detected': 0,
            'features_spanning_multiple_components': 0,
            'avg_coverage': 0.0,
            'avg_identity': 0.0,
            'samples_analyzed': len(impact_results),
            'sample_summaries': {}
        }
        
        total_coverage = 0.0
        total_identity = 0.0
        
        for sample in impact_results:
            sample_summary = {
                'haplotypes': len(impact_results[sample]),
                'total_features': 0,
                'impact_type_counts': {impact_type: 0 for impact_type in self.GENERAL_IMPACT_TYPES},
                'consequence_counts': defaultdict(int)
            }
            
            for haplotype in impact_results[sample]:
                for feature_name, analysis in impact_results[sample][haplotype].items():
                    summary['total_features'] += 1
                    sample_summary['total_features'] += 1
                    
                    impact_type = analysis['type']
                    consequence = analysis['consequence']
                    feature_type = analysis['details'].get('feature_type', 'unknown')
                    
                    # Count impact types
                    summary['impact_type_counts'][impact_type] += 1
                    sample_summary['impact_type_counts'][impact_type] += 1
                    
                    # Count consequences
                    summary['consequence_counts'][consequence] += 1
                    sample_summary['consequence_counts'][consequence] += 1
                    
                    # Feature type specific analysis
                    summary['feature_type_analysis'][feature_type]['total'] += 1
                    summary['feature_type_analysis'][feature_type]['impact_types'][impact_type] += 1
                    summary['feature_type_analysis'][feature_type]['consequences'][consequence] += 1
                    
                    details = analysis['details']
                    total_coverage += details.get('coverage', 0.0)
                    total_identity += details.get('identity', 0.0)
                    
                    if details.get('structural_variation_detected', False):
                        summary['structural_variations_detected'] += 1
                    
                    if len(details.get('components_spanned', [])) > 1:
                        summary['features_spanning_multiple_components'] += 1
            
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
                'gfa_file': str(self.gfa_file)
            },
            'impact_type_descriptions': {
                'general_types': self.GENERAL_IMPACT_TYPES,
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
    print(f"Impact type distribution:")
    for impact_type, count in summary['impact_type_counts'].items():
        percentage = (count / summary['total_features'] * 100) if summary['total_features'] > 0 else 0
        print(f"  {impact_type}: {count} ({percentage:.1f}%)")
    
    print(f"\nTop consequences:")
    sorted_consequences = sorted(summary['consequence_counts'].items(), 
                               key=lambda x: x[1], reverse=True)[:10]
    for consequence, count in sorted_consequences:
        percentage = (count / summary['total_features'] * 100) if summary['total_features'] > 0 else 0
        print(f"  {consequence}: {count} ({percentage:.1f}%)")
    
    print(f"\nStructural variations detected: {summary['structural_variations_detected']}")
    print(f"Features spanning multiple components: {summary['features_spanning_multiple_components']}")
    print(f"Average coverage: {summary['avg_coverage']:.3f}")
    print(f"Average identity: {summary['avg_identity']:.3f}")


if __name__ == "__main__":
    main()
