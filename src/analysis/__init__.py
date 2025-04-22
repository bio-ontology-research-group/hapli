"""
Annotation analysis module.

This package provides functionality to analyze feature alignments, detect
variants, reconcile feature hierarchies, and generate comprehensive summaries.
"""

from src.analysis.impact_classifier import ImpactClassifier, ImpactType
from src.analysis.variant_detector import VariantDetector, Variant, VariantType
from src.analysis.reconciliation import FeatureReconciler, ReconciliationStrategy
from src.analysis.summary_generator import SummaryGenerator, FeatureSummary, AnalysisSummary

__all__ = [
    'ImpactClassifier', 'ImpactType',
    'VariantDetector', 'Variant', 'VariantType',
    'FeatureReconciler', 'ReconciliationStrategy',
    'SummaryGenerator', 'FeatureSummary', 'AnalysisSummary'
]
