"""
TSV formatter for annotation reports.

This module provides the TSVFormatter class for generating
Tab-Separated Values (TSV) reports.
"""

import os
import csv
import logging
from typing import Dict, Any, List, Optional
from datetime import datetime

from src.analysis.summary_generator import AnalysisSummary

logger = logging.getLogger(__name__)

class TSVFormatter:
    """
    Formats analysis results as Tab-Separated Values (TSV).
    
    This class converts AnalysisSummary objects into TSV format
    for both single-path and comparative reports.
    """
    
    def __init__(self):
        """Initialize the TSV formatter."""
        pass
        
    def format(self, summary: AnalysisSummary, output_file: str) -> str:
        """
        Format an analysis summary as TSV.
        
        Args:
            summary: Analysis summary to format
            output_file: Path to output file
            
        Returns:
            Path to the generated file
        """
        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        
        # Define fieldnames
        fieldnames = [
            'feature_id', 'type', 'start', 'end', 'strand', 
            'impact_type', 'sequence_identity', 'coverage',
            'variants_snps', 'variants_insertions', 'variants_deletions', 
            'variants_complex', 'parent_features'
        ]
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            # Write metadata as a comment
            f.write(f"# Report generated: {datetime.now().isoformat()}\n")
            f.write(f"# Path: {summary.path_id}\n")
            f.write(f"# Features: {summary.feature_count}\n")
            f.write(f"# Present: {summary.feature_by_impact.get('present', 0)}\n")
            f.write(f"# Absent: {summary.feature_by_impact.get('absent', 0)}\n")
            f.write(f"# Modified: {summary.feature_by_impact.get('modified', 0)}\n\n")
            
            # Write each feature
            for feature_id, feature_summary in summary.feature_summaries.items():
                row = {
                    'feature_id': feature_id,
                    'type': feature_summary.feature_type,
                    'impact_type': feature_summary.impact_type.value if feature_summary.impact_type else 'UNKNOWN',
                    'sequence_identity': feature_summary.sequence_identity if hasattr(feature_summary, 'sequence_identity') else '',
                    'coverage': feature_summary.coverage if hasattr(feature_summary, 'coverage') else '',
                    'parent_features': ','.join(feature_summary.parent_features) if hasattr(feature_summary, 'parent_features') and feature_summary.parent_features else ''
                }
                
                # Add location if available
                if feature_summary.location:
                    row['start'] = feature_summary.location[0]
                    row['end'] = feature_summary.location[1]
                    row['strand'] = '+' if feature_summary.location[2] == 1 else '-'
                else:
                    row['start'] = ''
                    row['end'] = ''
                    row['strand'] = ''
                
                # Add variant counts
                if hasattr(feature_summary, 'variants'):
                    variant_types = {'SNP': 0, 'INSERTION': 0, 'DELETION': 0, 'COMPLEX': 0}
                    for variant in feature_summary.variants:
                        vtype = variant.variant_type.value.upper()
                        if vtype in variant_types:
                            variant_types[vtype] += 1
                        else:
                            variant_types['COMPLEX'] += 1
                    
                    row['variants_snps'] = variant_types['SNP']
                    row['variants_insertions'] = variant_types['INSERTION']
                    row['variants_deletions'] = variant_types['DELETION']
                    row['variants_complex'] = variant_types['COMPLEX']
                else:
                    row['variants_snps'] = ''
                    row['variants_insertions'] = ''
                    row['variants_deletions'] = ''
                    row['variants_complex'] = ''
                
                writer.writerow(row)
                
        logger.info(f"Generated TSV report at {output_file}")
        return output_file
        
    def format_comparative(self, comparative_data: Dict[str, Any], output_file: str) -> str:
        """
        Format comparative analysis as TSV.
        
        Args:
            comparative_data: Comparative data structure
            output_file: Path to output file
            
        Returns:
            Path to the generated file
        """
        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        
        # Extract path IDs
        path_ids = comparative_data['metadata']['paths']
        
        # Define fieldnames
        fieldnames = ['feature_id', 'type']
        for path_id in path_ids:
            path_prefix = f"{path_id}_"
            fieldnames.extend([
                f"{path_prefix}impact",
                f"{path_prefix}start",
                f"{path_prefix}end",
                f"{path_prefix}identity",
                f"{path_prefix}snps",
                f"{path_prefix}insertions",
                f"{path_prefix}deletions"
            ])
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            # Write metadata as a comment
            f.write(f"# Comparative report generated: {datetime.now().isoformat()}\n")
            f.write(f"# Paths: {', '.join(path_ids)}\n")
            f.write(f"# Reference: {comparative_data['metadata']['reference_id']}\n")
            f.write(f"# Features: {comparative_data['metadata']['feature_count']}\n\n")
            
            # Write each feature
            for feature_id, feature_data in comparative_data['features'].items():
                row = {'feature_id': feature_id, 'type': ''}  # Type will be filled from first available path
                
                for path_id in path_ids:
                    path_prefix = f"{path_id}_"
                    path_data = feature_data['paths'].get(path_id, {})
                    
                    # Set feature type from first available path that has this feature
                    if not row['type'] and 'type' in path_data and path_data['type']:
                        row['type'] = path_data['type']
                    
                    # Add path-specific data
                    row[f"{path_prefix}impact"] = path_data.get('impact_type', 'ABSENT')
                    row[f"{path_prefix}start"] = path_data.get('start', '')
                    row[f"{path_prefix}end"] = path_data.get('end', '')
                    row[f"{path_prefix}identity"] = path_data.get('sequence_identity', '')
                    
                    # Add variant counts
                    variants = path_data.get('variants', {}) or {}
                    row[f"{path_prefix}snps"] = variants.get('SNP', '')
                    row[f"{path_prefix}insertions"] = variants.get('INSERTION', '')
                    row[f"{path_prefix}deletions"] = variants.get('DELETION', '')
                
                writer.writerow(row)
                
        logger.info(f"Generated comparative TSV report at {output_file}")
        return output_file
