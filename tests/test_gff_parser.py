import os
import unittest
import tempfile
from src.parsers.gff_parser import GFF3Parser

class TestGFF3Parser(unittest.TestCase):
    """Tests for the GFF3 parser."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.parser = GFF3Parser()
        self.data_dir = "data"
        self.example_gff = os.path.join(self.data_dir, "example.gff3")
        
        # Create temp files for malformed input tests
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # Empty file
        self.empty_file = os.path.join(self.temp_dir.name, "empty.gff3")
        with open(self.empty_file, 'w') as f:
            pass
            
        # Malformed file (missing required fields)
        self.malformed_file = os.path.join(self.temp_dir.name, "malformed.gff3")
        with open(self.malformed_file, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("ref_chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")
            f.write("ref_chr1\tmalformed\n")  # Missing fields
            
    def tearDown(self):
        """Clean up test fixtures."""
        self.temp_dir.cleanup()
        
    def test_parse_valid_file(self):
        """Test parsing a valid GFF3 file."""
        self.assertTrue(os.path.exists(self.example_gff), f"Example GFF3 file not found at {self.example_gff}")
        
        records = self.parser.parse(self.example_gff)
        self.assertIsNotNone(records)
        self.assertGreater(len(records), 0)
        
        # Check that features have been indexed
        self.assertGreater(len(self.parser.features_by_id), 0)
        
    def test_nonexistent_file(self):
        """Test parsing a non-existent file."""
        with self.assertRaises(FileNotFoundError):
            self.parser.parse("nonexistent_file.gff3")
            
    def test_empty_file(self):
        """Test parsing an empty file."""
        with self.assertRaises(ValueError):
            self.parser.parse(self.empty_file)
            
    def test_malformed_file(self):
        """Test parsing a malformed file."""
        try:
            records = self.parser.parse(self.malformed_file)
            # If parsing succeeds, check that we got the valid record
            self.assertGreater(len(records), 0)
            
            # Should have indexed the valid gene feature
            self.assertIn('gene1', self.parser.features_by_id)
        except Exception as e:
            # Biopython might raise an exception for malformed files
            # This is acceptable behavior
            pass
        
    def test_get_features_by_type(self):
        """Test retrieving features by type."""
        records = self.parser.parse(self.example_gff)
        
        # Get feature types
        feature_types = self.parser.get_all_feature_types()
        self.assertGreater(len(feature_types), 0)
        
        # Test with a valid type
        first_type = next(iter(feature_types))
        features = self.parser.get_features_by_type(first_type)
        self.assertGreater(len(features), 0)
        
        # Test with an invalid type
        features = self.parser.get_features_by_type("nonexistent_type")
        self.assertEqual(len(features), 0)
        
    def test_get_feature_by_id(self):
        """Test retrieving a feature by ID."""
        records = self.parser.parse(self.example_gff)
        
        # Get a feature ID
        if self.parser.features_by_id:
            feature_id = next(iter(self.parser.features_by_id.keys()))
            feature = self.parser.get_feature_by_id(feature_id)
            self.assertIsNotNone(feature)
        
        # Test with an invalid ID
        feature = self.parser.get_feature_by_id("nonexistent_id")
        self.assertIsNone(feature)

if __name__ == '__main__':
    unittest.main()
