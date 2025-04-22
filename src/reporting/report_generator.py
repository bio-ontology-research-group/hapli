"""
Report generator for analysis results.

This module provides the ReportGenerator class that creates reports
from feature annotation analysis results.
"""

import os
import logging
from typing import Dict, List, Optional, Union, Any, Set
from datetime import datetime

from src.analysis.impact_classifier import ImpactType
from src.analysis.summary_generator import AnalysisSummary, FeatureSummary
from .formatters.tsv_formatter import TSVFormatter
from .formatters.json_formatter import JSONFormatter
from .formatters.rdf_formatter import RDFFormatter

logger = logging.getLogger(__name__)

class ReportGenerator:
    """
    Generates formatted reports from feature annotation analysis.
    
    This class provides functionality to generate reports in various
    formats (TSV, JSON, RDF) from analysis results, with filtering
    capabilities to focus on specific impact types or genomic regions.
    """
    
    def __init__(self):
        """Initialize the report generator."""
        self.formatters = {
            'tsv': TSVFormatter(),
            'json': JSONFormatter(),
            'rdf': RDFFormatter(),
        }
        
    def generate_report(self, 
                       summary: AnalysisSummary,
                       output_file: str,
                       format_type: str = 'tsv',
                       impact_types: Optional[List[ImpactType]] = None,
                       region: Optional[Dict[str, Union[int, int]]] = None,
                       rdf_format: str = 'turtle') -> str:
        """
        Generate a report from analysis summary.
        
        Args:
            summary: Analysis summary to report
            output_file: Path to output file
            format_type: Format type ('tsv', 'json', 'rdf')
            impact_types: List of impact types to include (None = all)
            region: Dict with 'start' and 'end' keys specifying a region to filter on
            rdf_format: RDF serialization format if format_type is 'rdf'
                       ('turtle', 'xml', 'json-ld', 'ntriples')
        
        Returns:
            Path to the generated report file
        
        Raises:
            ValueError: If format_type is not supported
        """
        if format_type not in self.formatters:
            raise ValueError(f"Unsupported format type: {format_type}")
        
        # Apply filters if provided
        filtered_summary = self._filter_summary(summary, impact_types, region)
        
        # Generate the report
        if format_type == 'rdf':
            output_file = self.formatters[format_type].format(
                filtered_summary, output_file, rdf_format=rdf_format)
        else:
            output_file = self.formatters[format_type].format(
                filtered_summary, output_file)
            
        logger.info(f"Generated {format_type.upper()} report: {output_file}")
        return output_file
    
    def _filter_summary(self, 
                      summary: AnalysisSummary,
                      impact_types: Optional[List[ImpactType]] = None,
                      region: Optional[Dict[str, int]] = None) -> AnalysisSummary:
        """
        Filter an analysis summary by impact types and/or region.
        
        Args:
            summary: Analysis summary to filter
            impact_types: List of impact types to include
            region: Dict with 'start' and 'end' keys specifying a region
        
        Returns:
            Filtered analysis summary
        """
        if not impact_types and not region:
            return summary
        
        # Create a new summary with the same metadata but without features
        filtered_feature_summaries = {}
        
        # Filter features
        for feature_id, feature_summary in summary.feature_summaries.items():
            # Filter by impact type if specified
            if impact_types and feature_summary.impact_type not in impact_types:
                continue
                
            # Filter by region if specified
            if region:
                feature_start = feature_summary.location[0] if feature_summary.location else 0
                feature_end = feature_summary.location[1] if feature_summary.location else 0
                
                if (feature_end < region['start'] or
                    feature_start > region['end']):
                    continue
            
            # Include this feature in the filtered summary
            filtered_feature_summaries[feature_id] = feature_summary
        
        # Create new analysis summary with filtered features
        filtered = AnalysisSummary(
            path_id=summary.path_id,
            feature_count=len(filtered_feature_summaries),
            feature_by_impact={},  # Will be recalculated based on filtered features
            variant_counts={},     # Will be recalculated based on filtered features
            reconciliation_counts={}, # Will be recalculated based on filtered features
            feature_summaries=filtered_feature_summaries
        )
        
        # Recalculate impact counts
        impact_counts = {}
        for fs in filtered_feature_summaries.values():
            if fs.impact_type:
                impact_type = fs.impact_type.value
                impact_counts[impact_type] = impact_counts.get(impact_type, 0) + 1
        
        filtered.feature_by_impact = impact_counts
            
        return filtered
        
    def register_formatter(self, format_type: str, formatter: Any) -> None:
        """
        Register a custom formatter.
        
        Args:
            format_type: Format type identifier
            formatter: Formatter instance
            
        Raises:
            ValueError: If formatter doesn't have a format method
        """
        if not hasattr(formatter, 'format'):
            raise ValueError("Formatter must have a 'format' method")
        
        self.formatters[format_type] = formatter
        
    def validate_rdf_report(self, report_file: str, shex_file: str) -> bool:
        """
        Validate an RDF report against a ShEx schema.
        
        Args:
            report_file: Path to RDF report file
            shex_file: Path to ShEx schema file
            
        Returns:
            True if validation succeeded, False otherwise
        """
        return self.formatters['rdf'].validate(report_file, shex_file)
