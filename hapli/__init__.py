"""
Hapli: Pangenome variant impact analysis toolkit.

This package provides tools for analyzing the impact of structural variations
on genomic features using pangenome graphs.
"""

__version__ = "0.1.0"
__author__ = "Hapli Development Team"

# Import main classes for convenience
from .gff_alignment import GFFAligner
from .gam_parser import GAMParser
from .impact_detector import ImpactDetector
from .diploid_analyzer import DiploidAnalyzer

__all__ = [
    "GFFAligner",
    "GAMParser", 
    "ImpactDetector",
    "DiploidAnalyzer"
]
