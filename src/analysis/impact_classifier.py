"""
Impact classification for feature alignments.

This module provides functionality to classify the impact of sequence
variations on genomic features by analyzing alignment results.
"""

import logging
from enum import Enum
from typing import Dict, List, Optional, Tuple

from Bio.SeqFeature import SeqFeature

logger = logging.getLogger(__name__)

class ImpactType(Enum):
    """Types of impacts that can occur in feature alignments."""
    PRESENT = "present"           # Feature is present with no significant changes
    ABSENT = "absent"             # Feature is completely missing
    MODIFIED = "modified"         # Feature is present but with sequence modifications
    TRUNCATED = "truncated"       # Feature is present but with start/end truncation
    EXPANDED = "expanded"         # Feature is present but with expanded boundaries
    FRAGMENTED = "fragmented"     # Feature is split across multiple locations
    INVERTED = "inverted"         # Feature is present but inverted
    DUPLICATED = "duplicated"     # Feature is present with duplications
    REARRANGED = "rearranged"     # Feature is extensively rearranged
    UNCERTAIN = "uncertain"       # Cannot determine impact with confidence

class ImpactClassifier:
    """
    Classifies the impact of variations on genomic features.
    
    This class analyzes alignments of genomic features to determine
    how sequence variations affect the feature (e.g., whether it's
    present, absent, modified, etc.).
    """
    
    def __init__(self, sequence_identity_threshold: float = 0.9, 
                 coverage_threshold: float = 0.8):
        """
        Initialize the ImpactClassifier.
        
        Args:
            sequence_identity_threshold: Minimum identity to consider a feature as "present"
            coverage_threshold: Minimum coverage to consider a feature as "present"
        """
        self.sequence_identity_threshold = sequence_identity_threshold
        self.coverage_threshold = coverage_threshold
    
    def classify_feature_impact(self, 
                               reference_feature: SeqFeature,
                               aligned_feature: Optional[SeqFeature],
                               alignment_score: Optional[float] = None,
                               coverage: Optional[float] = None,
                               identity: Optional[float] = None,
                               variants: Optional[List] = None) -> Tuple[ImpactType, Dict]:
        """
        Classify the impact on a single feature.
        
        Args:
            reference_feature: The original feature from the reference
            aligned_feature: The feature as aligned to the target path (or None if not aligned)
            alignment_score: Optional explicit alignment score
            coverage: Optional explicit coverage value
            identity: Optional explicit sequence identity value
            variants: Optional list of detected variants
            
        Returns:
            A tuple containing:
            - ImpactType: The classified impact type
            - Dict: Additional metadata about the impact
        """
        # If the feature didn't align at all
        if aligned_feature is None:
            return ImpactType.ABSENT, {"reason": "No alignment found"}
        
        # Extract alignment metrics if not provided
        if coverage is None and hasattr(aligned_feature, "qualifiers"):
            coverage = float(aligned_feature.qualifiers.get("coverage", [0])[0])
        
        if identity is None and hasattr(aligned_feature, "qualifiers"):
            identity = float(aligned_feature.qualifiers.get("identity", [0])[0])
        
        # Default values if still not available
        coverage = coverage or 0.0
        identity = identity or 0.0
        
        # First check for complex rearrangements if variants are provided
        if variants:
            complex_impact = self._classify_from_variants(variants)
            if complex_impact:
                impact_type, metadata = complex_impact
                metadata.update({
                    "coverage": coverage,
                    "identity": identity
                })
                return impact_type, metadata
        
        # Get feature lengths for structural analysis
        ref_length = reference_feature.location.end - reference_feature.location.start
        aligned_length = aligned_feature.location.end - aligned_feature.location.start
        
        # Check for truncation or expansion first
        if abs(ref_length - aligned_length) / ref_length > 0.1:
            if aligned_length < ref_length:
                return ImpactType.TRUNCATED, {
                    "coverage": coverage,
                    "identity": identity,
                    "ref_length": ref_length,
                    "aligned_length": aligned_length,
                    "length_ratio": aligned_length / ref_length
                }
            else:
                return ImpactType.EXPANDED, {
                    "coverage": coverage,
                    "identity": identity,
                    "ref_length": ref_length,
                    "aligned_length": aligned_length,
                    "length_ratio": aligned_length / ref_length
                }
                
        # Then check if the feature appears intact
        if (coverage >= self.coverage_threshold and 
            identity >= self.sequence_identity_threshold):
            return ImpactType.PRESENT, {
                "coverage": coverage,
                "identity": identity
            }
        
        # Check if the feature is fragmented (multiple partial alignments)
        if hasattr(aligned_feature, "sub_features") and aligned_feature.sub_features:
            return ImpactType.FRAGMENTED, {
                "coverage": coverage,
                "identity": identity,
                "fragments": len(aligned_feature.sub_features)
            }
        
        # If we reach here, the feature has sequence modifications but is roughly the same length
        return ImpactType.MODIFIED, {
            "coverage": coverage,
            "identity": identity
        }
    
    def _classify_from_variants(self, variants: List) -> Optional[Tuple[ImpactType, Dict]]:
        """
        Classify impact based on detected variants.
        
        Args:
            variants: List of detected variants
            
        Returns:
            Optional tuple of (ImpactType, metadata) if a clear impact is found,
            None otherwise
        """
        from src.analysis.variant_detector import VariantType
        
        # Check for specific complex variant types
        for variant in variants:
            if variant.variant_type == VariantType.INVERSION:
                return ImpactType.INVERTED, {
                    "inversion_length": variant.length,
                    "inversion_start": variant.position,
                    "inversion_end": variant.end_position
                }
                
            elif variant.variant_type == VariantType.DUPLICATION:
                return ImpactType.DUPLICATED, {
                    "duplication_length": variant.length,
                    "duplication_count": variant.metadata.get("dup_count", 2),
                    "duplication_start": variant.position
                }
                
            elif variant.variant_type == VariantType.TRANSLOCATION:
                return ImpactType.REARRANGED, {
                    "translocation_length": variant.length,
                    "original_position": variant.position,
                    "new_position": variant.metadata.get("new_position", 0)
                }
                
            elif variant.variant_type == VariantType.TANDEM_REPEAT:
                # For tandem repeats, decide based on the change in repeat count
                ref_count = variant.metadata.get("ref_count", 1)
                aln_count = variant.metadata.get("aln_count", 1)
                
                if aln_count > ref_count:
                    return ImpactType.EXPANDED, {
                        "repeat_pattern": variant.metadata.get("repeat_pattern", ""),
                        "ref_repeat_count": ref_count,
                        "aln_repeat_count": aln_count,
                        "expansion_factor": aln_count / ref_count
                    }
                else:
                    return ImpactType.TRUNCATED, {
                        "repeat_pattern": variant.metadata.get("repeat_pattern", ""),
                        "ref_repeat_count": ref_count,
                        "aln_repeat_count": aln_count,
                        "contraction_factor": ref_count / aln_count
                    }
        
        # If no specific complex variant was found
        return None
    
    def classify_feature_set(self, 
                           reference_features: Dict[str, SeqFeature],
                           aligned_features: Dict[str, List[SeqFeature]]) -> Dict[str, Tuple[ImpactType, Dict]]:
        """
        Classify impacts for a set of features.
        
        Args:
            reference_features: Dictionary of reference features keyed by ID
            aligned_features: Dictionary of aligned features keyed by ID (values are lists 
                              to handle multiple alignments)
            
        Returns:
            Dictionary mapping feature IDs to their impact classifications (type and metadata)
        """
        results = {}
        
        for feature_id, ref_feature in reference_features.items():
            # Get aligned version of this feature, if available
            aligned_list = aligned_features.get(feature_id, [])
            
            if not aligned_list:
                # Feature is completely absent
                results[feature_id] = (ImpactType.ABSENT, {"reason": "No alignment found"})
                continue
            
            # If multiple alignments exist, use the best one
            if len(aligned_list) > 1:
                # Sort by coverage and identity (if available)
                best_aligned = sorted(
                    aligned_list,
                    key=lambda f: (
                        float(f.qualifiers.get("coverage", [0])[0]) if "coverage" in f.qualifiers else 0,
                        float(f.qualifiers.get("identity", [0])[0]) if "identity" in f.qualifiers else 0
                    ),
                    reverse=True
                )[0]
            else:
                best_aligned = aligned_list[0]
            
            # Classify the impact using the best alignment
            impact_type, metadata = self.classify_feature_impact(ref_feature, best_aligned)
            
            # Add information about multiple alignments if present
            if len(aligned_list) > 1:
                metadata["multiple_alignments"] = len(aligned_list)
            
            results[feature_id] = (impact_type, metadata)
        
        return results
