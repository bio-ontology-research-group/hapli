import os
import unittest
import tempfile
from src.parsers.gfa_parser import GFAParser

class TestGFAParser(unittest.TestCase):
    """Tests for the GFA parser."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.parser = GFAParser()
        self.data_dir = "data"
        self.example_gfa = os.path.join(self.data_dir, "example.gfa")
        
        # Create temp files for malformed input tests
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # Empty file
        self.empty_file = os.path.join(self.temp_dir.name, "empty.gfa")
        with open(self.empty_file, 'w') as f:
            pass
            
        # Malformed file (missing required fields) - GFA2 format
        self.malformed_file = os.path.join(self.temp_dir.name, "malformed.gfa")
        with open(self.malformed_file, 'w') as f:
            f.write("H\tVN:Z:2.0\n")
            f.write("S\n")  # Missing required fields
            f.write("E\t*\ts1+\ts2+\t10\t10\t0\t0\t*\n")
            
    def tearDown(self):
        """Clean up test fixtures."""
        self.temp_dir.cleanup()
        
    def test_parse_valid_file(self):
        """Test parsing a valid GFA file."""
        self.assertTrue(os.path.exists(self.example_gfa), f"Example GFA file not found at {self.example_gfa}")
        
        gfa = self.parser.parse(self.example_gfa)
        self.assertIsNotNone(gfa)
        
        # Check that we have segments
        segments = self.parser.get_segments()
        self.assertGreater(len(segments), 0)
        
    def test_nonexistent_file(self):
        """Test parsing a non-existent file."""
        with self.assertRaises(FileNotFoundError):
            self.parser.parse("nonexistent_file.gfa")
            
    def test_empty_file(self):
        """Test parsing an empty file."""
        with self.assertRaises(ValueError):
            self.parser.parse(self.empty_file)
            
    def test_malformed_file(self):
        """Test parsing a malformed file."""
        # Don't log to stdout during tests
        import logging
        logging.getLogger('src.parsers.gfa_parser').setLevel(logging.CRITICAL)
        
        # Should not raise exception but return minimal GFA object
        gfa = self.parser.parse(self.malformed_file)
        self.assertIsNotNone(gfa)
        
        # Reset log level
        logging.getLogger('src.parsers.gfa_parser').setLevel(logging.WARNING)
        
    def test_get_segment_sequence(self):
        """Test retrieving segment sequences."""
        gfa = self.parser.parse(self.example_gfa)
        
        # Get the first segment ID
        segment_ids = list(self.parser.get_segments().keys())
        if segment_ids:
            segment_id = segment_ids[0]
            sequence = self.parser.get_segment_sequence(segment_id)
            # Either it has a sequence or it's explicitly None
            self.assertIsNotNone(sequence)
        
        # Test non-existent segment
        sequence = self.parser.get_segment_sequence("nonexistent_segment")
        self.assertIsNone(sequence)

if __name__ == '__main__':
    unittest.main()
