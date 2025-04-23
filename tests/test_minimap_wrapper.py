#!/usr/bin/env python3
"""
Unit tests for the minimap2 wrapper.
"""

import os
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.alignment.minimap_wrapper import MinimapAligner, AlignmentResult

# Test data paths
TEST_DATA_DIR = Path(__file__).parent / "data"
REFERENCE_FILE = TEST_DATA_DIR / "reference.fasta"
REGIONS_FILE = TEST_DATA_DIR / "regions.fasta"
SNP_VARIANTS_FILE = TEST_DATA_DIR / "snp_variants.fasta"
INSERTION_VARIANTS_FILE = TEST_DATA_DIR / "insertion_variants.fasta"
DELETION_VARIANTS_FILE = TEST_DATA_DIR / "deletion_variants.fasta"
COMPLEX_VARIANTS_FILE = TEST_DATA_DIR / "complex_variants.fasta"


class TestMinimapWrapper(unittest.TestCase):
    """Test cases for the MinimapAligner class."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test data if it doesn't exist."""
        # Check if test data exists
        if not all(os.path.exists(f) for f in [
            REFERENCE_FILE, REGIONS_FILE, SNP_VARIANTS_FILE,
            INSERTION_VARIANTS_FILE, DELETION_VARIANTS_FILE, COMPLEX_VARIANTS_FILE
        ]):
            # Run the prepare_test_sequences.py script
            print("Test data not found. Generating test data...")
            prepare_script = TEST_DATA_DIR / "prepare_test_sequences.py"
            if os.path.exists(prepare_script):
                os.system(f"python {prepare_script}")
            else:
                raise FileNotFoundError(f"Test data preparation script not found: {prepare_script}")
        
        # Load test data
        cls.reference = SeqIO.read(REFERENCE_FILE, "fasta")
        cls.regions = list(SeqIO.parse(REGIONS_FILE, "fasta"))
        cls.snp_variants = list(SeqIO.parse(SNP_VARIANTS_FILE, "fasta"))
        cls.insertion_variants = list(SeqIO.parse(INSERTION_VARIANTS_FILE, "fasta"))
        cls.deletion_variants = list(SeqIO.parse(DELETION_VARIANTS_FILE, "fasta"))
        cls.complex_variants = list(SeqIO.parse(COMPLEX_VARIANTS_FILE, "fasta"))
    
    def setUp(self):
        """Set up a fresh aligner for each test."""
        self.aligner = MinimapAligner(preset="map-ont")
    
    def test_load_reference_string(self):
        """Test loading a reference sequence as a string."""
        self.aligner.load_reference(str(self.reference.seq))
        self.assertIsNotNone(self.aligner.aligner)
    
    def test_load_reference_seqrecord(self):
        """Test loading a reference sequence as a SeqRecord."""
        self.aligner.load_reference(self.reference)
        self.assertIsNotNone(self.aligner.aligner)
    
    def test_load_reference_dict(self):
        """Test loading a reference sequence as a dictionary of SeqRecords."""
        self.aligner.load_reference({self.reference.id: self.reference})
        self.assertIsNotNone(self.aligner.aligner)
    
    def test_load_reference_file(self):
        """Test loading a reference sequence from a file."""
        self.aligner.load_reference_file(str(REFERENCE_FILE))
        self.assertIsNotNone(self.aligner.aligner)
    
    def test_load_reference_empty(self):
        """Test loading an empty reference sequence."""
        with self.assertRaises(ValueError):
            self.aligner.load_reference("")
    
    def test_load_reference_file_not_found(self):
        """Test loading a non-existent reference file."""
        with self.assertRaises(FileNotFoundError):
            self.aligner.load_reference_file("nonexistent_file.fasta")
    
    def test_align_sequence_basic(self):
        """Test basic sequence alignment."""
        self.aligner.load_reference(self.reference)
        
        # Align the first region (should align perfectly)
        alignments = self.aligner.align_sequence(str(self.regions[0].seq))
        
        self.assertGreater(len(alignments), 0)
        self.assertIsInstance(alignments[0], AlignmentResult)
        self.assertGreater(alignments[0].score, 0)
    
    def test_align_sequence_without_reference(self):
        """Test alignment without loading a reference."""
        with self.assertRaises(ValueError):
            self.aligner.align_sequence("ACGT")
    
    def test_align_empty_sequence(self):
        """Test alignment with an empty sequence."""
        self.aligner.load_reference(self.reference)
        with self.assertRaises(ValueError):
            self.aligner.align_sequence("")
    
    def test_align_snp_variants(self):
        """Test alignment of sequences with SNPs."""
        self.aligner.load_reference(self.reference)
        
        for i, variant in enumerate(self.snp_variants):
            # Get the original region
            original = self.regions[i]
            
            # Align the variant
            alignments = self.aligner.align_sequence(str(variant.seq))
            
            # Check that we got an alignment
            self.assertGreater(len(alignments), 0)
            
            # Check that the identity is high but not perfect
            self.assertLess(alignments[0].identity, 1.0)
            self.assertGreater(alignments[0].identity, 0.95)  # Expect >95% identity for SNP variants
    
    def test_align_insertion_variants(self):
        """Test alignment of sequences with insertions."""
        self.aligner.load_reference(self.reference)
        
        for i, variant in enumerate(self.insertion_variants):
            # Align the variant
            alignments = self.aligner.align_sequence(str(variant.seq))
            
            # Check that we got an alignment
            self.assertGreater(len(alignments), 0)
            
            # Check that the query length is greater than the target length
            self.assertGreater(
                alignments[0].query_end - alignments[0].query_start,
                alignments[0].target_end - alignments[0].target_start
            )
    
    def test_align_deletion_variants(self):
        """Test alignment of sequences with deletions."""
        self.aligner.load_reference(self.reference)
        
        for i, variant in enumerate(self.deletion_variants):
            # Align the variant
            alignments = self.aligner.align_sequence(str(variant.seq))
            
            # Check that we got an alignment
            self.assertGreater(len(alignments), 0)
            
            # Check that the query length is less than the target length
            self.assertLess(
                alignments[0].query_end - alignments[0].query_start,
                alignments[0].target_end - alignments[0].target_start
            )
    
    def test_align_complex_variants(self):
        """Test alignment of sequences with complex changes."""
        self.aligner.load_reference(self.reference)
        
        for i, variant in enumerate(self.complex_variants):
            # Align the variant
            alignments = self.aligner.align_sequence(str(variant.seq))
            
            # Check that we got an alignment
            self.assertGreater(len(alignments), 0)
            
            # Check that the identity is lower than for SNP variants
            self.assertLess(alignments[0].identity, 0.95)
    
    def test_alignment_parameters(self):
        """Test setting different alignment parameters."""
        # Test with different k-mer size
        aligner1 = MinimapAligner(preset="map-ont", k=15)
        aligner1.load_reference(self.reference)
        
        # Test with different window size
        aligner2 = MinimapAligner(preset="map-ont", w=10)
        aligner2.load_reference(self.reference)
        
        # Test with custom scoring
        aligner3 = MinimapAligner(preset="map-ont", scoring=(2, 4, 4, 2))  # match, mismatch, gap_open, gap_extend
        aligner3.load_reference(self.reference)
        
        # Align the same sequence with all aligners
        seq = str(self.regions[0].seq)
        
        alignments1 = aligner1.align_sequence(seq)
        alignments2 = aligner2.align_sequence(seq)
        alignments3 = aligner3.align_sequence(seq)
        
        # All should produce alignments
        self.assertGreater(len(alignments1), 0)
        self.assertGreater(len(alignments2), 0)
        self.assertGreater(len(alignments3), 0)
    
    def test_align_multiple_sequences(self):
        """Test aligning multiple sequences at once."""
        self.aligner.load_reference(self.reference)
        
        # Create a list of sequences
        sequences = [str(region.seq) for region in self.regions]
        
        # Align all sequences
        results = self.aligner.align_multiple_sequences(sequences)
        
        # Check that we got results for all sequences
        self.assertEqual(len(results), len(sequences))
        
        # Check that each sequence has at least one alignment
        for i, alignments in results.items():
            self.assertGreater(len(alignments), 0)
    
    def test_get_alignment_statistics(self):
        """Test getting alignment statistics."""
        self.aligner.load_reference(self.reference)
        
        # Align a sequence
        alignments = self.aligner.align_sequence(str(self.regions[0].seq))
        
        # Get statistics for the first alignment
        stats = self.aligner.get_alignment_statistics(alignments[0])
        
        # Check that we got the expected statistics
        self.assertIn('query_length', stats)
        self.assertIn('target_length', stats)
        self.assertIn('score', stats)
        self.assertIn('identity', stats)
        self.assertIn('strand', stats)
        self.assertIn('cigar', stats)
    
    def test_calculate_identity_from_cigar(self):
        """Test calculating identity from a CIGAR string."""
        # Test with a perfect match
        identity = self.aligner._calculate_identity_from_cigar("100M")
        self.assertEqual(identity, 1.0)
        
        # Test with mismatches/indels
        identity = self.aligner._calculate_identity_from_cigar("90M5I5D")
        self.assertEqual(identity, 0.9)  # 90 matches out of 100 total (90M + 5I + 5D)


if __name__ == "__main__":
    unittest.main()
