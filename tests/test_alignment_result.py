"""
Unit tests for the alignment result module.

This module tests the functionality of the AlignmentResult class
and related utilities for processing alignment results.
"""

import unittest
import json
import os
import sys
from typing import Dict, Any, Tuple

# Add the src directory to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.alignment.alignment_result import AlignmentResult, AlignmentStatistics, CigarOperation, AlignmentType
from src.alignment.terminal_display import AlignmentDisplay, display_alignment_result
from tests.data.alignment_test_data import AlignmentTestData


class TestAlignmentResult(unittest.TestCase):
    """Test cases for the AlignmentResult class."""
    
    def setUp(self):
        """Set up test data."""
        self.test_scenarios = AlignmentTestData.create_all_test_scenarios()
        self.expected_results = AlignmentTestData.get_expected_results()
    
    def test_from_minimap2(self):
        """Test creation of AlignmentResult from minimap2 alignment."""
        # Test with perfect match alignment
        aln, query_seq, target_seq = self.test_scenarios["perfect_match"]
        
        result = AlignmentResult.from_minimap2(aln, query_seq, target_seq)
        
        # Check basic properties
        self.assertEqual(result.query_name, "query1")
        self.assertEqual(result.target_name, "target1")
        self.assertEqual(result.query_sequence, query_seq)
        self.assertEqual(result.target_sequence, target_seq)
        self.assertEqual(result.query_start, 0)
        self.assertEqual(result.query_end, len(query_seq))
        self.assertEqual(result.target_start, 0)
        self.assertEqual(result.target_end, len(target_seq))
        self.assertEqual(result.score, 60)
        self.assertEqual(result.cigar_string, f"{len(query_seq)}M")
        self.assertFalse(result.is_reverse)
        
        # Check statistics
        self.assertAlmostEqual(result.statistics.identity, 1.0, places=2)
        self.assertAlmostEqual(result.statistics.coverage, 1.0, places=2)
        
        # Check alignment type
        self.assertEqual(result.statistics.alignment_type, AlignmentType.PERFECT)
    
    def test_cigar_parsing(self):
        """Test parsing of CIGAR strings."""
        # Test with a complex CIGAR string
        cigar = "10M2I5M3D"
        operations = AlignmentResult._parse_cigar(cigar)
        
        self.assertEqual(len(operations), 4)
        self.assertEqual(operations[0].operation, "M")
        self.assertEqual(operations[0].length, 10)
        self.assertEqual(operations[1].operation, "I")
        self.assertEqual(operations[1].length, 2)
        self.assertEqual(operations[2].operation, "M")
        self.assertEqual(operations[2].length, 5)
        self.assertEqual(operations[3].operation, "D")
        self.assertEqual(operations[3].length, 3)
    
    def test_alignment_statistics(self):
        """Test calculation of alignment statistics."""
        # Test with different alignment scenarios
        for scenario, (aln, query_seq, target_seq) in self.test_scenarios.items():
            if aln is None:  # Skip no_alignment scenario
                continue
                
            result = AlignmentResult.from_minimap2(aln, query_seq, target_seq)
            expected = self.expected_results[scenario]
            
            # Check identity
            if "identity" in expected:
                self.assertAlmostEqual(result.statistics.identity, expected["identity"], places=2)
            elif "identity_range" in expected:
                min_id, max_id = expected["identity_range"]
                self.assertTrue(min_id <= result.statistics.identity <= max_id)
            
            # Check coverage
            if "coverage" in expected:
                self.assertAlmostEqual(result.statistics.coverage, expected["coverage"], places=2)
            elif "coverage_range" in expected:
                min_cov, max_cov = expected["coverage_range"]
                self.assertTrue(min_cov <= result.statistics.coverage <= max_cov)
            
            # Check reverse flag
            if "is_reverse" in expected:
                self.assertEqual(result.is_reverse, expected["is_reverse"])
    
    def test_alignment_visualization(self):
        """Test generation of alignment visualization."""
        # Test with SNP alignment
        aln, query_seq, target_seq = self.test_scenarios["snp"]
        result = AlignmentResult.from_minimap2(aln, query_seq, target_seq)
        
        # Check that visualization strings were generated
        self.assertIsNotNone(result.aligned_query)
        self.assertIsNotNone(result.aligned_target)
        self.assertIsNotNone(result.alignment_indicator)
        
        # Check lengths match
        self.assertEqual(len(result.aligned_query), len(result.aligned_target))
        self.assertEqual(len(result.aligned_query), len(result.alignment_indicator))
        
        # Check for mismatch indicators
        self.assertIn(".", result.alignment_indicator)
    
    def test_serialization(self):
        """Test serialization and deserialization of AlignmentResult."""
        # Test with insertion alignment
        aln, query_seq, target_seq = self.test_scenarios["insertion"]
        original = AlignmentResult.from_minimap2(aln, query_seq, target_seq)
        
        # Convert to JSON and back
        json_str = original.to_json()
        data = json.loads(json_str)
        restored = AlignmentResult.from_dict(data)
        
        # Check that key properties are preserved
        self.assertEqual(restored.query_name, original.query_name)
        self.assertEqual(restored.target_name, original.target_name)
        self.assertEqual(restored.query_sequence, original.query_sequence)
        self.assertEqual(restored.target_sequence, original.target_sequence)
        self.assertEqual(restored.query_start, original.query_start)
        self.assertEqual(restored.query_end, original.query_end)
        self.assertEqual(restored.target_start, original.target_start)
        self.assertEqual(restored.target_end, original.target_end)
        self.assertEqual(restored.score, original.score)
        self.assertEqual(restored.cigar_string, original.cigar_string)
        self.assertEqual(restored.is_reverse, original.is_reverse)
        
        # Check statistics
        self.assertAlmostEqual(restored.statistics.identity, original.statistics.identity, places=2)
        self.assertAlmostEqual(restored.statistics.coverage, original.statistics.coverage, places=2)
    
    def test_summary_generation(self):
        """Test generation of human-readable summary."""
        # Test with perfect match alignment
        aln, query_seq, target_seq = self.test_scenarios["perfect_match"]
        result = AlignmentResult.from_minimap2(aln, query_seq, target_seq)
        
        summary = result.get_summary()
        
        # Check that summary contains key information
        self.assertIn("query1", summary)
        self.assertIn("target1", summary)
        self.assertIn("Identity: 100.00%", summary)
        self.assertIn("Coverage: 100.00%", summary)
        self.assertIn("perfect", summary.lower())
    
    def test_no_alignment(self):
        """Test handling of no alignment case."""
        # Get the no alignment scenario
        _, query_seq, target_seq = self.test_scenarios["no_alignment"]
        
        # Create a minimal AlignmentResult with no alignment
        result = AlignmentResult(
            query_name="query1",
            target_name="target1",
            query_sequence=query_seq,
            target_sequence=target_seq
        )
        
        # Check statistics
        self.assertEqual(result.statistics.identity, 0.0)
        self.assertEqual(result.statistics.coverage, 0.0)
        self.assertEqual(result.statistics.alignment_type, AlignmentType.NO_ALIGNMENT)
        
        # Check summary
        summary = result.get_summary()
        self.assertIn("Identity: 0.00%", summary)
        self.assertIn("Coverage: 0.00%", summary)
        self.assertIn("no_alignment", summary.lower())


class TestTerminalDisplay(unittest.TestCase):
    """Test cases for the terminal display functionality."""
    
    def setUp(self):
        """Set up test data."""
        self.test_scenarios = AlignmentTestData.create_all_test_scenarios()
        self.display = AlignmentDisplay(use_color=False)  # Disable color for testing
    
    def test_display_alignment(self):
        """Test display of a single alignment."""
        # Test with SNP alignment
        aln, query_seq, target_seq = self.test_scenarios["snp"]
        result = AlignmentResult.from_minimap2(aln, query_seq, target_seq)
        
        # Test detailed display
        detailed = self.display.display_alignment(result, detailed=True)
        
        # Check that output contains key elements
        self.assertIn("Alignment: query1 to target1", detailed)
        self.assertIn("Identity:", detailed)
        self.assertIn("Coverage:", detailed)
        
        # Check that alignment visualization is included
        self.assertTrue(any(line.strip() and line[0].isdigit() for line in detailed.split('\n')))
        
        # Test summary display
        summary = self.display.display_alignment(result, detailed=False)
        
        # Check that summary is shorter but still contains key information
        self.assertIn("Alignment: query1 to target1", summary)
        self.assertIn("Identity:", summary)
        self.assertLess(len(summary), len(detailed))
    
    def test_display_multiple_alignments(self):
        """Test display of multiple alignments."""
        # Create a list of alignment results
        alignments = []
        for scenario, (aln, query_seq, target_seq) in self.test_scenarios.items():
            if aln is not None:
                alignments.append(AlignmentResult.from_minimap2(aln, query_seq, target_seq))
        
        # Test multiple alignment display
        output = self.display.display_multiple_alignments(alignments, detailed=False)
        
        # Check that output contains information about each alignment
        self.assertIn(f"Displaying {len(alignments)} alignments", output)
        for i in range(1, len(alignments) + 1):
            self.assertIn(f"Alignment {i}:", output)
    
    def test_compact_summary(self):
        """Test compact summary of multiple alignments."""
        # Create a list of alignment results
        alignments = []
        for scenario, (aln, query_seq, target_seq) in self.test_scenarios.items():
            if aln is not None:
                alignments.append(AlignmentResult.from_minimap2(aln, query_seq, target_seq))
        
        # Test compact summary
        output = self.display.display_compact_summary(alignments)
        
        # Check that output contains a table header and rows
        self.assertIn("Summary of", output)
        self.assertIn("Query", output)
        self.assertIn("Target", output)
        self.assertIn("Identity", output)
        self.assertIn("Coverage", output)
        self.assertIn("Score", output)
        self.assertIn("Type", output)
        
        # Check that there's a row for each alignment
        lines = output.split('\n')
        # Find lines that start with a digit (row numbers)
        data_rows = [line for line in lines if line.strip() and line.strip()[0].isdigit()]
        self.assertEqual(len(data_rows), len(alignments))


if __name__ == '__main__':
    unittest.main()
