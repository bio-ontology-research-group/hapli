"""
Comparative report generator for multi-haplotype analysis.

This module provides the ComparativeReportGenerator class that creates 
reports comparing features across multiple haplotypes.
"""

import os
import logging
from typing import Dict, List, Optional, Set, Tuple, Any
from datetime import datetime

from src.analysis.impact_classifier import ImpactType
from src.analysis.summary_generator import AnalysisSummary
from .formatters.tsv_formatter import TSVFormatter
from .formatters.json_formatter import JSONFormatter
from .formatters.rdf_formatter import RDFFormatter

logger = logging.getLogger(__name__)

class ComparativeReportGenerator:
    """
    Generates reports comparing features across multiple haplotypes.
    
    This class provides functionality to generate comparative reports
    that highlight similarities and differences in feature annotations
    across multiple haplotypes or samples.
    """
    
    def __init__(self):
        """Initialize the comparative report generator."""
        self.formatters = {
            'tsv': TSVFormatter(),
            'json': JSONFormatter(),
            'rdf': RDFFormatter(),
        }
        
    def generate_comparative_report(self,
                                  summaries: Dict[str, AnalysisSummary],
                                  output_file: str,
                                  format_type: str = 'tsv',
                                  impact_types: Optional[List[ImpactType]] = None,
                                  region: Optional[Dict[str, int]] = None,
                                  rdf_format: str = 'turtle') -> str:
        """
        Generate a comparative report from multiple analysis summaries.
        
        Args:
            summaries: Dict mapping path_id to AnalysisSummary
            output_file: Path to output file
            format_type: Format type ('tsv', 'json', 'rdf')
            impact_types: List of impact types to include (None = all)
            region: Dict with 'start' and 'end' keys specifying a region to filter on
            rdf_format: RDF serialization format if format_type is 'rdf'
                       ('turtle', 'xml', 'json-ld', 'ntriples')
        
        Returns:
            Path to the generated report file
            
        Raises:
            ValueError: If format_type is not supported or if summaries is empty
        """
        if not summaries:
            raise ValueError("No summaries provided for comparison")
            
        if format_type not in self.formatters:
            raise ValueError(f"Unsupported format type: {format_type}")
        
        # Create the comparative data structure
        comparative_data = self._create_comparative_data(summaries, impact_types, region)
        
        # Generate the report
        if format_type == 'rdf':
            output_file = self.formatters[format_type].format_comparative(
                comparative_data, output_file, rdf_format=rdf_format)
        else:
            output_file = self.formatters[format_type].format_comparative(
                comparative_data, output_file)
            
        logger.info(f"Generated comparative {format_type.upper()} report: {output_file}")
        return output_file
    
    def _create_comparative_data(self,
                               summaries: Dict[str, AnalysisSummary],
                               impact_types: Optional[List[ImpactType]] = None,
                               region: Optional[Dict[str, int]] = None) -> Dict[str, Any]:
        """
        Create a comparative data structure from multiple summaries.
        
        Args:
            summaries: Dict mapping path_id to AnalysisSummary
            impact_types: List of impact types to include
            region: Dict with 'start' and 'end' keys specifying a region
            
        Returns:
            Comparative data structure for reporting
        """
        # Find all unique feature IDs across all paths
        all_feature_ids = set()
        for summary in summaries.values():
            all_feature_ids.update(summary.feature_summaries.keys())
            
        # Create comparative data structure
        comparative_data = {
            'metadata': {
                'timestamp': datetime.now().isoformat(),
                'paths': list(summaries.keys()),
                'reference_id': next(iter(summaries.values())).path_id.split('_')[0] if '_' in next(iter(summaries.values())).path_id else "reference",
                'feature_count': len(all_feature_ids)
            },
            'features': {}
        }
        
        # Populate feature data for each path
        for feature_id in sorted(all_feature_ids):
            feature_data = {
                'feature_id': feature_id,
                'paths': {}
            }
            
            # Get data for this feature in each path
            for path_id, summary in summaries.items():
                if feature_id in summary.feature_summaries:
                    feature_summary = summary.feature_summaries[feature_id]
                    
                    # Apply filters
                    if impact_types and feature_summary.impact_type not in impact_types:
                        continue
                        
                    if region:
                        feature_start = feature_summary.location[0] if feature_summary.location else 0
                        feature_end = feature_summary.location[1] if feature_summary.location else 0
                        
                        if (feature_end < region['start'] or
                            feature_start > region['end']):
                            continue
                    
                    # Add to comparative data
                    path_data = {
                        'impact_type': feature_summary.impact_type.name if feature_summary.impact_type else 'UNKNOWN',
                        'type': feature_summary.feature_type if hasattr(feature_summary, 'feature_type') else '',
                        'start': feature_summary.location[0] if feature_summary.location else None,
                        'end': feature_summary.location[1] if feature_summary.location else None,
                        'strand': feature_summary.location[2] if feature_summary.location and len(feature_summary.location) > 2 else None,
                        'sequence_identity': feature_summary.sequence_identity if hasattr(feature_summary, 'sequence_identity') else None,
                        'coverage': feature_summary.coverage if hasattr(feature_summary, 'coverage') else None,
                    }
                    
                    # Add variant data if available
                    if hasattr(feature_summary, 'variants'):
                        variants_by_type = {}
                        for variant in feature_summary.variants:
                            variant_type = variant.variant_type.value
                            variants_by_type[variant_type] = variants_by_type.get(variant_type, 0) + 1
                        path_data['variants'] = variants_by_type
                    
                    feature_data['paths'][path_id] = path_data
                else:
                    feature_data['paths'][path_id] = {
                        'impact_type': 'ABSENT',
                        'type': '',
                        'start': None,
                        'end': None,
                        'strand': None,
                        'variants': None,
                        'sequence_identity': None,
                        'coverage': None
                    }
            
            # Only include if the feature appears in at least one path after filtering
            if any(feature_data['paths'].values()):
                comparative_data['features'][feature_id] = feature_data
                
        return comparative_data
    
    def identify_consensus_features(self,
                                  summaries: Dict[str, AnalysisSummary],
                                  threshold: float = 0.9) -> Set[str]:
        """
        Identify features present in most paths (consensus features).
        
        Args:
            summaries: Dict mapping path_id to AnalysisSummary
            threshold: Proportion of paths that must have the feature (0.0-1.0)
            
        Returns:
            Set of feature IDs that are present in at least threshold proportion of paths
        """
        if not summaries:
            return set()
            
        # Count occurrences of each feature
        feature_counts = {}
        for summary in summaries.values():
            for feature_id, feature_summary in summary.feature_summaries.items():
                if feature_summary.impact_type != ImpactType.ABSENT:
                    feature_counts[feature_id] = feature_counts.get(feature_id, 0) + 1
        
        # Find features meeting threshold
        path_count = len(summaries)
        min_paths = int(path_count * threshold)
        consensus_features = {
            feature_id for feature_id, count in feature_counts.items()
            if count >= min_paths
        }
        
        return consensus_features
    
    def identify_discriminating_features(self,
                                       summaries: Dict[str, AnalysisSummary]) -> Dict[str, List[str]]:
        """
        Identify features that discriminate between paths.
        
        Args:
            summaries: Dict mapping path_id to AnalysisSummary
            
        Returns:
            Dict mapping feature IDs to lists of path IDs where they are present
        """
        if not summaries:
            return {}
            
        # Group paths by feature presence
        feature_to_paths = {}
        for path_id, summary in summaries.items():
            for feature_id, feature_summary in summary.feature_summaries.items():
                if feature_summary.impact_type != ImpactType.ABSENT:
                    if feature_id not in feature_to_paths:
                        feature_to_paths[feature_id] = []
                    feature_to_paths[feature_id].append(path_id)
        
        # Keep only features not present in all paths
        discriminating_features = {
            feature_id: paths for feature_id, paths in feature_to_paths.items()
            if len(paths) < len(summaries) and len(paths) > 0
        }
        
        return discriminating_features
