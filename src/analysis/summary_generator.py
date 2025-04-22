"""
Summary generation for feature annotation analysis.

This module provides functionality to aggregate and summarize annotation
analysis results, combining data from impact classification, variant detection,
and reconciliation processes.
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from Bio.SeqFeature import SeqFeature

from src.analysis.impact_classifier import ImpactType
from src.analysis.variant_detector import Variant, VariantType
from src.analysis.reconciliation import ReconciliationResult, ReconciliationStrategy

logger = logging.getLogger(__name__)

@dataclass
class FeatureSummary:
    """Summary of analysis results for a single feature."""
    feature_id: str
    feature_type: str
    impact_type: Optional[ImpactType] = None
    variants: List[Variant] = None
    reconciliation: Optional[ReconciliationResult] = None
    parent_features: List[str] = None
    child_features: List[str] = None
    sequence_identity: Optional[float] = None
    coverage: Optional[float] = None
    path_id: Optional[str] = None
    location: Optional[Tuple[int, int, int]] = None  # start, end, strand
    
    def __post_init__(self):
        if self.variants is None:
            self.variants = []
        if self.parent_features is None:
            self.parent_features = []
        if self.child_features is None:
            self.child_features = []

@dataclass
class AnalysisSummary:
    """Overall summary of feature annotation analysis."""
    path_id: str
    feature_count: int
    feature_by_impact: Dict[str, int]
    variant_counts: Dict[str, int]
    reconciliation_counts: Dict[str, int]
    feature_summaries: Dict[str, FeatureSummary]
    
    @property
    def features_present(self) -> int:
        """Number of features classified as present."""
        return self.feature_by_impact.get(ImpactType.PRESENT.value, 0)
    
    @property
    def features_absent(self) -> int:
        """Number of features classified as absent."""
        return self.feature_by_impact.get(ImpactType.ABSENT.value, 0)
    
    @property
    def features_modified(self) -> int:
        """Number of features classified as modified."""
        return self.feature_by_impact.get(ImpactType.MODIFIED.value, 0)

class SummaryGenerator:
    """
    Generates summaries of feature annotation analysis.
    
    This class aggregates results from impact classification, variant detection,
    and reconciliation to provide comprehensive summaries of how features are
    represented in target paths compared to the reference.
    """
    
    def __init__(self):
        """Initialize the SummaryGenerator."""
        pass
    
    def generate_feature_summary(self,
                               feature_id: str,
                               feature: SeqFeature,
                               path_id: str,
                               impact_result: Optional[Tuple[ImpactType, Dict]] = None,
                               variants: Optional[List[Variant]] = None,
                               reconciliation: Optional[ReconciliationResult] = None,
                               parent_features: Optional[List[str]] = None,
                               child_features: Optional[List[str]] = None) -> FeatureSummary:
        """
        Generate a summary for a single feature.
        
        Args:
            feature_id: ID of the feature
            feature: The feature object
            path_id: ID of the path
            impact_result: Result from impact classification (optional)
            variants: List of variants detected in the feature (optional)
            reconciliation: Result from feature reconciliation (optional)
            parent_features: List of parent feature IDs (optional)
            child_features: List of child feature IDs (optional)
            
        Returns:
            FeatureSummary object
        """
        # Extract feature type
        feature_type = feature.type
        
        # Extract location information
        location = (
            int(feature.location.start),
            int(feature.location.end),
            feature.location.strand
        )
        
        # Extract impact information if available
        impact_type = None
        sequence_identity = None
        coverage = None
        
        if impact_result:
            impact_type, metadata = impact_result
            sequence_identity = metadata.get("identity")
            coverage = metadata.get("coverage")
        
        return FeatureSummary(
            feature_id=feature_id,
            feature_type=feature_type,
            impact_type=impact_type,
            variants=variants or [],
            reconciliation=reconciliation,
            parent_features=parent_features or [],
            child_features=child_features or [],
            sequence_identity=sequence_identity,
            coverage=coverage,
            path_id=path_id,
            location=location
        )
    
    def generate_path_summary(self, 
                            path_id: str,
                            feature_summaries: Dict[str, FeatureSummary]) -> AnalysisSummary:
        """
        Generate a summary for all features in a path.
        
        Args:
            path_id: ID of the path
            feature_summaries: Dictionary of feature summaries by ID
            
        Returns:
            AnalysisSummary object
        """
        # Count features by impact type
        impact_counts = {}
        for summary in feature_summaries.values():
            if summary.impact_type:
                impact_type = summary.impact_type.value
                impact_counts[impact_type] = impact_counts.get(impact_type, 0) + 1
        
        # Count variants by type
        variant_counts = {}
        for summary in feature_summaries.values():
            for variant in summary.variants:
                var_type = variant.variant_type.value
                variant_counts[var_type] = variant_counts.get(var_type, 0) + 1
        
        # Count reconciliation strategies
        reconciliation_counts = {}
        for summary in feature_summaries.values():
            if summary.reconciliation:
                strategy = summary.reconciliation.strategy.value
                reconciliation_counts[strategy] = reconciliation_counts.get(strategy, 0) + 1
        
        return AnalysisSummary(
            path_id=path_id,
            feature_count=len(feature_summaries),
            feature_by_impact=impact_counts,
            variant_counts=variant_counts,
            reconciliation_counts=reconciliation_counts,
            feature_summaries=feature_summaries
        )
    
    def export_summary_tsv(self, 
                         summary: AnalysisSummary, 
                         output_file: str) -> None:
        """
        Export a summary to a TSV file.
        
        Args:
            summary: The AnalysisSummary to export
            output_file: Path to the output file
        """
        with open(output_file, 'w') as f:
            # Write header
            f.write("feature_id\tfeature_type\tpath_id\timpact\tidentity\tcoverage\tstart\tend\t" +
                    "strand\tvariants\treconciliation\tparents\tchildren\n")
            
            # Write data for each feature
            for feature_id, fs in summary.feature_summaries.items():
                impact = fs.impact_type.value if fs.impact_type else "unknown"
                identity = f"{fs.sequence_identity:.2f}" if fs.sequence_identity is not None else "NA"
                coverage = f"{fs.coverage:.2f}" if fs.coverage is not None else "NA"
                
                if fs.location:
                    start, end, strand = fs.location
                else:
                    start, end, strand = "NA", "NA", "NA"
                
                # Format variants
                variant_str = ";".join([f"{v.variant_type.value}:{v.position}" for v in fs.variants])
                variant_str = variant_str or "none"
                
                # Format reconciliation
                if fs.reconciliation:
                    reconciliation_str = fs.reconciliation.strategy.value
                else:
                    reconciliation_str = "none"
                
                # Format parents and children
                parents_str = ";".join(fs.parent_features) or "none"
                children_str = ";".join(fs.child_features) or "none"
                
                # Write the line
                f.write(f"{feature_id}\t{fs.feature_type}\t{fs.path_id}\t{impact}\t{identity}\t" +
                        f"{coverage}\t{start}\t{end}\t{strand}\t{variant_str}\t" +
                        f"{reconciliation_str}\t{parents_str}\t{children_str}\n")
    
    def export_summary_json(self, 
                          summary: AnalysisSummary, 
                          output_file: str) -> None:
        """
        Export a summary to a JSON file.
        
        Args:
            summary: The AnalysisSummary to export
            output_file: Path to the output file
        """
        import json
        
        # Convert summary to a dictionary
        summary_dict = {
            "path_id": summary.path_id,
            "feature_count": summary.feature_count,
            "impact_counts": summary.feature_by_impact,
            "variant_counts": summary.variant_counts,
            "reconciliation_counts": summary.reconciliation_counts,
            "features": {}
        }
        
        # Convert each feature summary
        for feature_id, fs in summary.feature_summaries.items():
            feature_dict = {
                "id": fs.feature_id,
                "type": fs.feature_type,
                "path_id": fs.path_id,
                "impact": fs.impact_type.value if fs.impact_type else None,
                "identity": fs.sequence_identity,
                "coverage": fs.coverage,
                "location": fs.location,
                "parents": fs.parent_features,
                "children": fs.child_features,
                "variants": [
                    {
                        "type": v.variant_type.value,
                        "position": v.position,
                        "reference": v.reference,
                        "alternate": v.alternate,
                        "length": v.length,
                        "quality": v.quality
                    } for v in fs.variants
                ]
            }
            
            # Add reconciliation info if available
            if fs.reconciliation:
                feature_dict["reconciliation"] = {
                    "strategy": fs.reconciliation.strategy.value,
                    "description": fs.reconciliation.description,
                    "confidence": fs.reconciliation.confidence
                }
            
            summary_dict["features"][feature_id] = feature_dict
        
        # Write to file
        with open(output_file, 'w') as f:
            json.dump(summary_dict, f, indent=2)
