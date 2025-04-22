"""
Reporting module for generating and formatting annotation reports.

This module provides functionality for generating reports about feature
annotations and their alignments to paths in a variation graph.
"""

from .report_generator import ReportGenerator
from .comparative_report import ComparativeReportGenerator

__all__ = ['ReportGenerator', 'ComparativeReportGenerator']
