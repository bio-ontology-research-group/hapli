import os
import unittest
import tempfile
from src.parsers.fasta_parser import FastaParser

class TestFastaParser(unittest.TestCase):
    """Tests for the FASTA parser."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.parser = FastaParser()
        self.data_dir = "data"
        self.reference_fasta = os.path.join(self.data_dir, "reference.fasta")
        
        # Create temp files for malformed input tests
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # Empty file
        self.empty_file = os.path.join(self.temp_dir.name, "empty.fasta")
        with open(self.empty_file, 'w') as f:
            pass
            
        # Malformed file (not in FASTA format)
        self.malformed_file = os.path.join(self.temp_dir.name, "malformed.fasta")
        with open(self.malformed_file, 'w') as f:
            f.write("This is not a FASTA file\n")
            
    def tearDown(self):
        """Clean up test fixtures."""
        self.temp_dir.cleanup()
        
    def test_parse_valid_file(self):
        """Test parsing a valid FASTA file."""
        self.assertTrue(os.path.exists(self.reference_fasta), 
                        f"Reference FASTA file not found at {self.reference_fasta}")
        
        sequences = self.parser.parse(self.reference_fasta)
        self.assertIsNotNone(sequences)
        self.assertGreater(len(sequences), 0)
        
    def test_nonexistent_file(self):
        """Test parsing a non-existent file."""
        with self.assertRaises(FileNotFoundError):
            self.parser.parse("nonexistent_file.fasta")
            
    def test_empty_file(self):
        """Test parsing an empty file."""
        with self.assertRaises(ValueError):
            self.parser.parse(self.empty_file)
            
    def test_malformed_file(self):
        """Test parsing a malformed file."""
        with self.assertRaises(ValueError):
            self.parser.parse(self.malformed_file)
        
    def test_get_sequence(self):
        """Test retrieving a sequence by ID."""
        sequences = self.parser.parse(self.reference_fasta)
        
        # Get the first sequence ID
        seq_id = next(iter(sequences.keys()))
        sequence = self.parser.get_sequence(seq_id)
        self.assertIsNotNone(sequence)
        
        # Test with a non-existent ID
        sequence = self.parser.get_sequence("nonexistent_id")
        self.assertIsNone(sequence)
        
    def test_get_sequence_string(self):
        """Test retrieving a sequence as a string."""
        sequences = self.parser.parse(self.reference_fasta)
        
        # Get the first sequence ID
        seq_id = next(iter(sequences.keys()))
        sequence_str = self.parser.get_sequence_string(seq_id)
        self.assertIsNotNone(sequence_str)
        self.assertIsInstance(sequence_str, str)
        
        # Test with a non-existent ID
        sequence_str = self.parser.get_sequence_string("nonexistent_id")
        self.assertIsNone(sequence_str)
        
    def test_get_sequence_length(self):
        """Test retrieving the length of a sequence."""
        sequences = self.parser.parse(self.reference_fasta)
        
        # Get the first sequence ID
        seq_id = next(iter(sequences.keys()))
        length = self.parser.get_sequence_length(seq_id)
        self.assertIsNotNone(length)
        self.assertGreater(length, 0)
        
        # Test with a non-existent ID
        length = self.parser.get_sequence_length("nonexistent_id")
        self.assertIsNone(length)
        
    def test_get_all_sequences(self):
        """Test retrieving all sequences."""
        sequences = self.parser.parse(self.reference_fasta)
        all_sequences = self.parser.get_all_sequences()
        self.assertEqual(len(sequences), len(all_sequences))

if __name__ == '__main__':
    unittest.main()
