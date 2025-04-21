import unittest
import os
import tempfile
from unittest.mock import patch, MagicMock

from src.path_analysis import PathAnalyzer

class TestPathAnalysis(unittest.TestCase):
    """Tests for the path analysis module."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.TemporaryDirectory()
        self.analyzer = PathAnalyzer()
        
        # Create a mock GFA object
        self.mock_gfa = MagicMock()
        
        # Create mock paths with different naming patterns
        # Using naming patterns that match our sample data generator
        self.mock_paths = {
            'sample1_h1': MagicMock(),
            'sample1_h2': MagicMock(),
            'sample2_h1': MagicMock(),
            'sample2_h2': MagicMock(),
            'sample3.1': MagicMock(),
            'sample3.2': MagicMock(),
            'hap1_sample4': MagicMock(),
            'hap2_sample4': MagicMock(),
            'ref_chr1': MagicMock(),  # Reference sequence path (matches REF_SEQ_ID in the generator)
            'complex_name_with_no_pattern': MagicMock()
        }
        
        # Set up path segment names
        for path_id, path in self.mock_paths.items():
            path.segment_names = [f"seg1_{path_id}", f"seg2_{path_id}", f"seg3_{path_id}"]
        
        # Assign the mock paths to the mock GFA object
        self.mock_gfa.paths = self.mock_paths
    
    def tearDown(self):
        """Clean up test fixtures."""
        self.test_dir.cleanup()
    
    def test_extract_paths(self):
        """Test that paths are correctly extracted from a GFA object."""
        self.analyzer.load_gfa(self.mock_gfa)
        
        # Check that all paths were extracted
        self.assertEqual(len(self.analyzer.paths), len(self.mock_paths))
        for path_id in self.mock_paths:
            self.assertIn(path_id, self.analyzer.paths)
    
    def test_group_paths_by_sample(self):
        """Test that paths are correctly grouped by sample name."""
        self.analyzer.load_gfa(self.mock_gfa)
        sample_groups = self.analyzer.group_paths_by_sample()
        
        # Check that we have the expected number of groups
        expected_samples = {'sample1', 'sample2', 'sample3', 'sample4', 'ref_chr1', 'complex_name_with_no_pattern'}
        self.assertEqual(len(sample_groups), len(expected_samples))
        
        # Check that each group contains the correct paths
        self.assertEqual(set(sample_groups['sample1']), {'sample1_h1', 'sample1_h2'})
        self.assertEqual(set(sample_groups['sample2']), {'sample2_h1', 'sample2_h2'})
        self.assertEqual(set(sample_groups['sample3']), {'sample3.1', 'sample3.2'})
        self.assertEqual(set(sample_groups['sample4']), {'hap1_sample4', 'hap2_sample4'})
        
        # Paths without a recognizable pattern should be in their own group
        self.assertEqual(sample_groups['ref_chr1'], ['ref_chr1'])
        self.assertEqual(sample_groups['complex_name_with_no_pattern'], ['complex_name_with_no_pattern'])
    
    def test_identify_haplotypes(self):
        """Test that haplotype relationships are correctly identified."""
        self.analyzer.load_gfa(self.mock_gfa)
        self.analyzer.group_paths_by_sample()
        haplotype_groups = self.analyzer.identify_haplotypes()
        
        # Check that haplotypes are correctly identified for each sample
        sample1_haplotypes = dict(haplotype_groups['sample1'])
        self.assertEqual(sample1_haplotypes['sample1_h1'], '1')
        self.assertEqual(sample1_haplotypes['sample1_h2'], '2')
        
        sample2_haplotypes = dict(haplotype_groups['sample2'])
        self.assertEqual(sample2_haplotypes['sample2_h1'], '1')
        self.assertEqual(sample2_haplotypes['sample2_h2'], '2')
        
        sample3_haplotypes = dict(haplotype_groups['sample3'])
        self.assertEqual(sample3_haplotypes['sample3.1'], '1')
        self.assertEqual(sample3_haplotypes['sample3.2'], '2')
        
        sample4_haplotypes = dict(haplotype_groups['sample4'])
        self.assertEqual(sample4_haplotypes['hap1_sample4'], '1')
        self.assertEqual(sample4_haplotypes['hap2_sample4'], '2')
        
        # Paths without a recognizable haplotype pattern should default to haplotype "1"
        reference_haplotypes = dict(haplotype_groups['ref_chr1'])
        self.assertEqual(reference_haplotypes['ref_chr1'], '1')
        
        complex_haplotypes = dict(haplotype_groups['complex_name_with_no_pattern'])
        self.assertEqual(complex_haplotypes['complex_name_with_no_pattern'], '1')
    
    def test_select_paths_by_sample(self):
        """Test selecting paths by sample name."""
        self.analyzer.load_gfa(self.mock_gfa)
        self.analyzer.group_paths_by_sample()
        self.analyzer.identify_haplotypes()
        
        # Select paths for specific samples
        selected_paths = self.analyzer.select_paths(sample_names=['sample1', 'sample3'])
        
        # Check that only paths from the specified samples are selected
        expected_paths = {'sample1_hap1', 'sample1_hap2', 'sample3.1', 'sample3.2'}
        self.assertEqual(set(selected_paths), expected_paths)
    
    def test_select_paths_by_haplotype(self):
        """Test selecting paths by haplotype ID."""
        self.analyzer.load_gfa(self.mock_gfa)
        self.analyzer.group_paths_by_sample()
        self.analyzer.identify_haplotypes()
        
        # Select paths for specific haplotypes
        selected_paths = self.analyzer.select_paths(haplotype_ids=['1'])
        
        # Check that only haplotype 1 paths are selected
        expected_paths = {'sample1_h1', 'sample2_h1', 'sample3.1', 'hap1_sample4', 'ref_chr1', 'complex_name_with_no_pattern'}
        self.assertEqual(set(selected_paths), expected_paths)
    
    def test_select_paths_by_sample_and_haplotype(self):
        """Test selecting paths by both sample name and haplotype ID."""
        self.analyzer.load_gfa(self.mock_gfa)
        self.analyzer.group_paths_by_sample()
        self.analyzer.identify_haplotypes()
        
        # Select paths for specific samples and haplotypes
        selected_paths = self.analyzer.select_paths(
            sample_names=['sample1', 'sample2'],
            haplotype_ids=['2']
        )
        
        # Check that only the specified sample+haplotype combinations are selected
        expected_paths = {'sample1_hap2', 'sample2_h2'}
        self.assertEqual(set(selected_paths), expected_paths)
    
    def test_select_specific_paths(self):
        """Test selecting specific paths by ID."""
        self.analyzer.load_gfa(self.mock_gfa)
        
        # Select specific paths by ID
        selected_paths = self.analyzer.select_paths(
            path_ids=['sample1_hap1', 'reference', 'nonexistent_path']
        )
        
        # Check that only the existing specified paths are selected
        expected_paths = {'sample1_hap1', 'reference'}
        self.assertEqual(set(selected_paths), expected_paths)
    
    def test_get_path_segments(self):
        """Test retrieving segments for a specific path."""
        self.analyzer.load_gfa(self.mock_gfa)
        
        # Get segments for a specific path
        segments = self.analyzer.get_path_segments('sample1_hap1')
        
        # Check that the correct segments are returned
        expected_segments = ['seg1_sample1_hap1', 'seg2_sample1_hap1', 'seg3_sample1_hap1']
        self.assertEqual(segments, expected_segments)
    
    def test_edge_case_single_path(self):
        """Test with a GFA that has only a single path."""
        # Create a mock GFA with a single path
        single_path_gfa = MagicMock()
        single_path_gfa.paths = {'ref_chr1': MagicMock()}
        single_path_gfa.paths['ref_chr1'].segment_names = ['seg1', 'seg2', 'seg3']
        
        # Analyze the paths
        self.analyzer.load_gfa(single_path_gfa)
        sample_groups = self.analyzer.group_paths_by_sample()
        
        # Check that we have only one group
        self.assertEqual(len(sample_groups), 1)
        self.assertEqual(sample_groups['ref_chr1'], ['ref_chr1'])
    
    def test_edge_case_no_paths(self):
        """Test with a GFA that has no paths."""
        # Create a mock GFA with no paths
        no_paths_gfa = MagicMock()
        no_paths_gfa.paths = {}
        
        # Analyze the paths
        self.analyzer.load_gfa(no_paths_gfa)
        sample_groups = self.analyzer.group_paths_by_sample()
        
        # Check that we have no groups
        self.assertEqual(len(sample_groups), 0)
        
        # Test selection with empty path list
        selected_paths = self.analyzer.select_paths(
            sample_names=['sample1'],
            haplotype_ids=['1']
        )
        self.assertEqual(len(selected_paths), 0)
