"""
JSON formatter for annotation reports.

This module provides the JSONFormatter class for generating
JavaScript Object Notation (JSON) reports.
"""

import os
import json
import logging
from typing import Dict, Any, List, Optional
from datetime import datetime

from src.analysis.summary_generator import AnalysisSummary

logger = logging.getLogger(__name__)

class JSONFormatter:
    """
    Formats analysis results as JSON.
    
    This class converts AnalysisSummary objects into JSON format
    for both single-path and comparative reports.
    """
    
    def __init__(self):
        """Initialize the JSON formatter."""
        pass
        
    def format(self, summary: AnalysisSummary, output_file: str) -> str:
        """
        Format an analysis summary as JSON.
        
        Args:
            summary: Analysis summary to format
            output_file: Path to output file
            
        Returns:
            Path to the generated file
        """
        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        
        # Create JSON structure
        json_data = {
            'metadata': {
                'timestamp': datetime.now().isoformat(),
                'path_id': summary.path_id,
                'reference_id': summary.path_id.split('_')[0] if '_' in summary.path_id else "reference",
                'feature_count': summary.feature_count,
                'features_present': summary.feature_by_impact.get('present', 0),
                'features_absent': summary.feature_by_impact.get('absent', 0),
                'features_modified': summary.feature_by_impact.get('modified', 0)
            },
            'features': {}
        }
        
        # Add feature data
        for feature_id, feature_summary in summary.feature_summaries.items():
            feature_data = {
                'type': feature_summary.feature_type,
                'impact_type': feature_summary.impact_type.value if feature_summary.impact_type else 'UNKNOWN',
                'sequence_identity': feature_summary.sequence_identity if hasattr(feature_summary, 'sequence_identity') else None,
                'coverage': feature_summary.coverage if hasattr(feature_summary, 'coverage') else None,
                'parent_features': feature_summary.parent_features if hasattr(feature_summary, 'parent_features') else []
            }
            
            # Add location if available
            if feature_summary.location:
                feature_data['start'] = feature_summary.location[0]
                feature_data['end'] = feature_summary.location[1]
                feature_data['strand'] = '+' if feature_summary.location[2] == 1 else '-'
            
            # Add variant data
            if hasattr(feature_summary, 'variants') and feature_summary.variants:
                # Summarize variants by type
                variants_summary = {}
                for variant in feature_summary.variants:
                    vtype = variant.variant_type.value
                    variants_summary[vtype] = variants_summary.get(vtype, 0) + 1
                
                feature_data['variants'] = variants_summary
                
                # Add detailed variants
                feature_data['variant_details'] = [
                    {
                        'type': v.variant_type.value,
                        'position': v.position,
                        'reference': v.reference,
                        'alternate': v.alternate,
                        'length': v.length,
                        'quality': v.quality
                    }
                    for v in feature_summary.variants
                ]
            
            json_data['features'][feature_id] = feature_data
            
        # Write to file
        with open(output_file, 'w') as f:
            json.dump(json_data, f, indent=2)
            
        logger.info(f"Generated JSON report at {output_file}")
        return output_file
        
    def format_comparative(self, comparative_data: Dict[str, Any], output_file: str) -> str:
        """
        Format comparative analysis as JSON.
        
        Args:
            comparative_data: Comparative data structure
            output_file: Path to output file
            
        Returns:
            Path to the generated file
        """
        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        
        # Write directly to file
        with open(output_file, 'w') as f:
            json.dump(comparative_data, f, indent=2)
            
        logger.info(f"Generated comparative JSON report at {output_file}")
        return output_file
