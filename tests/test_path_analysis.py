import unittest
import os
import sys
import tempfile
import importlib.util
from unittest.mock import patch, MagicMock

from src.path_analysis import PathAnalyzer
from src.parsers.gfa_parser import GFAParser

class TestPathAnalysis(unittest.TestCase):
    """Tests for the path analysis module."""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test class with generated data."""
        # Load the data generator script as a module
        script_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 
                                  "scripts", "generate_example_data.py")
        spec = importlib.util.spec_from_file_location("data_generator", script_path)
        cls.data_generator = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.data_generator)
    
    def setUp(self):
        """Set up test fixtures with real, generated test data."""
        # Create temporary directory for test files
        self.test_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.test_dir.name
        
        # Generate test data in temporary directory
        self.data_generator.generate_data(
            ref_length=1000,
            num_features=5,
            paired_haplotypes=True,
            output_dir=self.output_dir,
            basename="test"
        )
        
        # Define paths to generated files
        self.gfa_path = os.path.join(self.output_dir, "test.gfa")
        self.gff_path = os.path.join(self.output_dir, "test.gff3")
        self.fasta_path = os.path.join(self.output_dir, "test.fasta")
        
        # Verify files were generated
        self.assertTrue(os.path.exists(self.gfa_path), f"GFA file not found at {self.gfa_path}")
        self.assertTrue(os.path.exists(self.gff_path), f"GFF file not found at {self.gff_path}")
        self.assertTrue(os.path.exists(self.fasta_path), f"FASTA file not found at {self.fasta_path}")
        
        # Initialize a GFA parser and a path analyzer
        self.gfa_parser = GFAParser()
        self.analyzer = PathAnalyzer()
        
        # Parse the GFA file
        self.gfa = self.gfa_parser.parse(self.gfa_path)
        
        # Use mocks for specific tests that need controlled data
        self.setup_mock_data()
    
    def setup_mock_data(self):
        """Set up additional mock data for specific test cases."""
        # Create a mock GFA object
        self.mock_gfa = MagicMock()
        
        # Create mock paths with different naming patterns
        self.mock_paths = {
            'sample1_h1': MagicMock(),
            'sample1_h2': MagicMock(),
            'sample2_h1': MagicMock(),
            'sample2_h2': MagicMock(),
            'sample3.1': MagicMock(),
            'sample3.2': MagicMock(),
            'hap1_sample4': MagicMock(),
            'hap2_sample4': MagicMock(),
            'ref_chr1': MagicMock(),
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
        # Test with real generated data
        paths = self.analyzer.load_gfa(self.gfa)
        
        # We should have at least the reference path and sample paths
        self.assertGreaterEqual(len(paths), 5)  # ref_chr1 + sample1_h1, sample1_h2, sample2_h1, sample2_h2
        self.assertIn('ref_chr1', paths)
        self.assertIn('sample1_h1', paths)
        self.assertIn('sample1_h2', paths)
        self.assertIn('sample2_h1', paths)
        self.assertIn('sample2_h2', paths)
        
        # Also test with controlled mock data
        self.analyzer.load_gfa(self.mock_gfa)
        
        # Check that all mock paths were extracted
        self.assertEqual(len(self.analyzer.paths), len(self.mock_paths))
        for path_id in self.mock_paths:
            self.assertIn(path_id, self.analyzer.paths)
    
    def test_group_paths_by_sample(self):
        """Test that paths are correctly grouped by sample name."""
        # Test with real generated data
        self.analyzer.load_gfa(self.gfa)
        sample_groups = self.analyzer.group_paths_by_sample()
        
        # The data generator should have created 3 sample groups: ref_chr1, sample1, sample2
        self.assertGreaterEqual(len(sample_groups), 3)
        self.assertIn('ref_chr1', sample_groups)
        self.assertIn('sample1', sample_groups)
        self.assertIn('sample2', sample_groups)
        
        # Check that the right paths are in each group
        self.assertEqual(sample_groups['ref_chr1'], ['ref_chr1'])
        self.assertEqual(set(sample_groups['sample1']), {'sample1_h1', 'sample1_h2'})
        self.assertEqual(set(sample_groups['sample2']), {'sample2_h1', 'sample2_h2'})
        
        # Also test with controlled mock data
        self.analyzer.load_gfa(self.mock_gfa)
        mock_sample_groups = self.analyzer.group_paths_by_sample()
        
        # Check that we have the expected number of groups
        expected_samples = {'sample1', 'sample2', 'sample3', 'sample4', 'ref_chr1', 'complex_name_with_no_pattern'}
        self.assertEqual(len(mock_sample_groups), len(expected_samples))
        
        # Check that each group contains the correct paths
        self.assertEqual(set(mock_sample_groups['sample1']), {'sample1_h1', 'sample1_h2'})
        self.assertEqual(set(mock_sample_groups['sample2']), {'sample2_h1', 'sample2_h2'})
        self.assertEqual(set(mock_sample_groups['sample3']), {'sample3.1', 'sample3.2'})
        self.assertEqual(set(mock_sample_groups['sample4']), {'hap1_sample4', 'hap2_sample4'})
        
        # Paths without a recognizable pattern should be in their own group
        self.assertEqual(mock_sample_groups['ref_chr1'], ['ref_chr1'])
        self.assertEqual(mock_sample_groups['complex_name_with_no_pattern'], ['complex_name_with_no_pattern'])
    
    def test_identify_haplotypes(self):
        """Test that haplotype relationships are correctly identified."""
        # Test with real generated data
        self.analyzer.load_gfa(self.gfa)
        self.analyzer.group_paths_by_sample()
        haplotype_groups = self.analyzer.identify_haplotypes()
        
        # Check real data haplotypes
        self.assertIn('sample1', haplotype_groups)
        self.assertIn('sample2', haplotype_groups)
        
        sample1_haps = dict(haplotype_groups['sample1'])
        self.assertEqual(sample1_haps['sample1_h1'], '1')
        self.assertEqual(sample1_haps['sample1_h2'], '2')
        
        sample2_haps = dict(haplotype_groups['sample2'])
        self.assertEqual(sample2_haps['sample2_h1'], '1')
        self.assertEqual(sample2_haps['sample2_h2'], '2')
        
        # Reference should have default haplotype "1"
        ref_haps = dict(haplotype_groups['ref_chr1'])
        self.assertEqual(ref_haps['ref_chr1'], '1')
        
        # Also test with controlled mock data
        self.analyzer.load_gfa(self.mock_gfa)
        self.analyzer.group_paths_by_sample()
        mock_haplotype_groups = self.analyzer.identify_haplotypes()
        
        # Check that haplotypes are correctly identified for each mock sample
        sample1_haplotypes = dict(mock_haplotype_groups['sample1'])
        self.assertEqual(sample1_haplotypes['sample1_h1'], '1')
        self.assertEqual(sample1_haplotypes['sample1_h2'], '2')
        
        sample2_haplotypes = dict(mock_haplotype_groups['sample2'])
        self.assertEqual(sample2_haplotypes['sample2_h1'], '1')
        self.assertEqual(sample2_haplotypes['sample2_h2'], '2')
        
        sample3_haplotypes = dict(mock_haplotype_groups['sample3'])
        self.assertEqual(sample3_haplotypes['sample3.1'], '1')
        self.assertEqual(sample3_haplotypes['sample3.2'], '2')
        
        sample4_haplotypes = dict(mock_haplotype_groups['sample4'])
        self.assertEqual(sample4_haplotypes['hap1_sample4'], '1')
        self.assertEqual(sample4_haplotypes['hap2_sample4'], '2')
        
        # Paths without a recognizable haplotype pattern should default to haplotype "1"
        reference_haplotypes = dict(mock_haplotype_groups['ref_chr1'])
        self.assertEqual(reference_haplotypes['ref_chr1'], '1')
        
        complex_haplotypes = dict(mock_haplotype_groups['complex_name_with_no_pattern'])
        self.assertEqual(complex_haplotypes['complex_name_with_no_pattern'], '1')
    
    def test_select_paths_by_sample(self):
        """Test selecting paths by sample name."""
        # Test with real generated data
        self.analyzer.load_gfa(self.gfa)
        self.analyzer.group_paths_by_sample()
        self.analyzer.identify_haplotypes()
        
        # Select sample1 paths
        selected_paths = self.analyzer.select_paths(sample_names=['sample1'])
        
        # Check selection from real data
        self.assertEqual(len(selected_paths), 2)
        self.assertIn('sample1_h1', selected_paths)
        self.assertIn('sample1_h2', selected_paths)
        
        # Also test with controlled mock data
        self.analyzer.load_gfa(self.mock_gfa)
        self.analyzer.group_paths_by_sample()
        self.analyzer.identify_haplotypes()
        
        # Select paths for specific samples
        mock_selected_paths = self.analyzer.select_paths(sample_names=['sample1', 'sample3'])
        
        # Check that only paths from the specified samples are selected
        expected_paths = {'sample1_h1', 'sample1_h2', 'sample3.1', 'sample3.2'}
        self.assertEqual(set(mock_selected_paths), expected_paths)
    
    def test_select_paths_by_haplotype(self):
        """Test selecting paths by haplotype ID."""
        # Test with real generated data
        self.analyzer.load_gfa(self.gfa)
        self.analyzer.group_paths_by_sample()
        self.analyzer.identify_haplotypes()
        
        # Select haplotype 1 paths
        selected_paths = self.analyzer.select_paths(haplotype_ids=['1'])
        
        # Check selection from real data
        self.assertGreaterEqual(len(selected_paths), 3)  # At least ref_chr1, sample1_h1, sample2_h1
        self.assertIn('ref_chr1', selected_paths)
        self.assertIn('sample1_h1', selected_paths)
        self.assertIn('sample2_h1', selected_paths)
        
        # Also test with controlled mock data
        self.analyzer.load_gfa(self.mock_gfa)
        self.analyzer.group_paths_by_sample()
        self.analyzer.identify_haplotypes()
        
        # Select paths for specific haplotypes
        mock_selected_paths = self.analyzer.select_paths(haplotype_ids=['1'])
        
        # Check that only haplotype 1 paths are selected
        expected_paths = {'sample1_h1', 'sample2_h1', 'sample3.1', 'hap1_sample4', 'ref_chr1', 'complex_name_with_no_pattern'}
        self.assertEqual(set(mock_selected_paths), expected_paths)
    
    def test_select_paths_by_sample_and_haplotype(self):
        """Test selecting paths by both sample name and haplotype ID."""
        # Test with real generated data
        self.analyzer.load_gfa(self.gfa)
        self.analyzer.group_paths_by_sample()
        self.analyzer.identify_haplotypes()
        
        # Select paths for specific samples and haplotypes
        selected_paths = self.analyzer.select_paths(
            sample_names=['sample1', 'sample2'],
            haplotype_ids=['2']
        )
        
        # Check selection from real data
        self.assertEqual(len(selected_paths), 2)
        self.assertIn('sample1_h2', selected_paths)
        self.assertIn('sample2_h2', selected_paths)
        
        # Also test with controlled mock data
        self.analyzer.load_gfa(self.mock_gfa)
        self.analyzer.group_paths_by_sample()
        self.analyzer.identify_haplotypes()
        
        # Select paths for specific samples and haplotypes
        mock_selected_paths = self.analyzer.select_paths(
            sample_names=['sample1', 'sample2'],
            haplotype_ids=['2']
        )
        
        # Check that only the specified sample+haplotype combinations are selected
        expected_paths = {'sample1_h2', 'sample2_h2'}
        self.assertEqual(set(mock_selected_paths), expected_paths)
    
    def test_select_specific_paths(self):
        """Test selecting specific paths by ID."""
        # Test with real generated data
        self.analyzer.load_gfa(self.gfa)
        
        # Select specific paths by ID
        selected_paths = self.analyzer.select_paths(
            path_ids=['sample1_h1', 'ref_chr1', 'nonexistent_path']
        )
        
        # Check selection from real data
        self.assertEqual(len(selected_paths), 2)
        self.assertIn('sample1_h1', selected_paths)
        self.assertIn('ref_chr1', selected_paths)
        
        # Also test with controlled mock data
        self.analyzer.load_gfa(self.mock_gfa)
        
        # Select specific paths by ID
        mock_selected_paths = self.analyzer.select_paths(
            path_ids=['sample1_h1', 'ref_chr1', 'nonexistent_path']
        )
        
        # Check that only the existing specified paths are selected
        expected_paths = {'sample1_h1', 'ref_chr1'}
        self.assertEqual(set(mock_selected_paths), expected_paths)
    
    def test_get_path_segments(self):
        """Test retrieving segments for a specific path."""
        # Test with real generated data
        self.analyzer.load_gfa(self.gfa)
        
        # Get segments for a specific path
        segments = self.analyzer.get_path_segments('sample1_h1')
        
        # Check segments from real data
        self.assertGreaterEqual(len(segments), 1)
        # Each segment should start with 's' followed by a number as per data generator
        for segment in segments:
            self.assertTrue(segment.startswith('s') and segment[1:].split('_')[0].isdigit(), 
                           f"Segment ID {segment} doesn't match expected format")
        
        # Also test with controlled mock data
        self.analyzer.load_gfa(self.mock_gfa)
        
        # Get segments for a specific path
        mock_segments = self.analyzer.get_path_segments('sample1_h1')
        
        # Check that the correct segments are returned
        expected_segments = ['seg1_sample1_h1', 'seg2_sample1_h1', 'seg3_sample1_h1']
        self.assertEqual(mock_segments, expected_segments)
    
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
