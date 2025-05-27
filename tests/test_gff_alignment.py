#!/usr/bin/env python3
"""
Test GFF3 topological sorting and sequence extraction functionality.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import gffutils

from hapli.gff_alignment import GFFAligner, setup_logging


class TestGFFAlignment(unittest.TestCase):
    """Test GFF alignment functionality."""
    
    def setUp(self):
        """Set up test data."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.reference_file = self.test_dir / "reference.fa"
        self.gff_file = self.test_dir / "features.gff3"
        
        # Create test reference genome
        self.create_test_reference()
        
        # Create test GFF3 file
        self.create_test_gff()
    
    def tearDown(self):
        """Clean up test data."""
        shutil.rmtree(self.test_dir)
    
    def create_test_reference(self):
        """Create a simple test reference genome."""
        # Create a simple chromosome with known sequence
        sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 50  # 2100 bp
        
        record = SeqRecord(
            Seq(sequence),
            id="chr1",
            description="Test chromosome 1"
        )
        
        with open(self.reference_file, 'w') as f:
            SeqIO.write([record], f, "fasta")
    
    def create_test_gff(self):
        """Create a test GFF3 file with hierarchical features."""
        gff_content = """##gff-version 3
##sequence-region chr1 1 2100
chr1	test	gene	100	1000	.	+	.	ID=gene1;Name=test_gene1
chr1	test	mRNA	100	1000	.	+	.	ID=mRNA1;Parent=gene1;Name=transcript1
chr1	test	exon	100	300	.	+	.	ID=exon1;Parent=mRNA1;exon_number=1
chr1	test	exon	500	700	.	+	.	ID=exon2;Parent=mRNA1;exon_number=2
chr1	test	exon	800	1000	.	+	.	ID=exon3;Parent=mRNA1;exon_number=3
chr1	test	CDS	150	300	.	+	0	ID=cds1;Parent=mRNA1
chr1	test	CDS	500	700	.	+	0	ID=cds2;Parent=mRNA1
chr1	test	CDS	800	950	.	+	0	ID=cds3;Parent=mRNA1
chr1	test	gene	1200	1800	.	-	.	ID=gene2;Name=test_gene2
chr1	test	mRNA	1200	1800	.	-	.	ID=mRNA2;Parent=gene2;Name=transcript2
chr1	test	exon	1200	1400	.	-	.	ID=exon4;Parent=mRNA2;exon_number=1
chr1	test	exon	1600	1800	.	-	.	ID=exon5;Parent=mRNA2;exon_number=2
chr1	test	CDS	1250	1400	.	-	0	ID=cds4;Parent=mRNA2
chr1	test	CDS	1600	1750	.	-	0	ID=cds5;Parent=mRNA2
"""
        
        with open(self.gff_file, 'w') as f:
            f.write(gff_content)
    
    def test_load_gff_database(self):
        """Test loading GFF3 database."""
        aligner = GFFAligner(self.gff_file, self.reference_file)
        db = aligner.load_gff_database()
        
        self.assertIsNotNone(db)
        self.assertIsInstance(db, gffutils.FeatureDB)
        
        # Check that features are loaded
        genes = list(db.features_of_type('gene'))
        self.assertEqual(len(genes), 2)
    
    def test_load_reference_genome(self):
        """Test loading reference genome."""
        aligner = GFFAligner(self.gff_file, self.reference_file)
        genome = aligner.load_reference_genome()
        
        self.assertIsNotNone(genome)
        self.assertIn('chr1', genome)
        self.assertEqual(len(genome['chr1'].seq), 2100)
    
    def test_topological_sort_features(self):
        """Test topological sorting of features."""
        aligner = GFFAligner(self.gff_file, self.reference_file)
        aligner.load_gff_database()
        
        sorted_features = aligner.topological_sort_features()
        
        self.assertIsNotNone(sorted_features)
        self.assertGreater(len(sorted_features), 0)
        
        # Check that genes come before their children in the sort
        feature_types = [f.featuretype for f in sorted_features]
        
        # Genes should appear before mRNAs, which should appear before exons/CDS
        gene_indices = [i for i, ft in enumerate(feature_types) if ft == 'gene']
        mrna_indices = [i for i, ft in enumerate(feature_types) if ft == 'mRNA']
        exon_indices = [i for i, ft in enumerate(feature_types) if ft == 'exon']
        
        # Check hierarchical ordering (genes before mRNAs before exons)
        if gene_indices and mrna_indices:
            self.assertLess(min(gene_indices), min(mrna_indices))
        if mrna_indices and exon_indices:
            self.assertLess(min(mrna_indices), min(exon_indices))
    
    def test_extract_feature_sequence(self):
        """Test extracting sequence for a feature."""
        aligner = GFFAligner(self.gff_file, self.reference_file)
        aligner.load_gff_database()
        aligner.load_reference_genome()
        
        # Get a gene feature
        gene = list(aligner.db.features_of_type('gene'))[0]
        sequence = aligner.extract_feature_sequence(gene)
        
        self.assertIsNotNone(sequence)
        self.assertIsInstance(sequence, str)
        self.assertGreater(len(sequence), 0)
        
        # Check sequence length matches feature length
        expected_length = gene.end - gene.start + 1
        self.assertEqual(len(sequence), expected_length)
    
    def test_get_feature_sequences(self):
        """Test getting sequences for multiple features."""
        aligner = GFFAligner(self.gff_file, self.reference_file)
        aligner.load_gff_database()
        aligner.load_reference_genome()
        aligner.topological_sort_features()
        
        feature_sequences = aligner.get_feature_sequences()
        
        self.assertIsNotNone(feature_sequences)
        self.assertGreater(len(feature_sequences), 0)
        
        # Check that each item is a (feature, sequence) tuple
        for feature, sequence in feature_sequences:
            self.assertIsNotNone(feature)
            self.assertIsInstance(sequence, str)
            self.assertGreater(len(sequence), 0)
    
    def test_process_gff_alignment(self):
        """Test complete workflow."""
        aligner = GFFAligner(self.gff_file, self.reference_file)
        
        # Run complete workflow
        feature_sequences = aligner.process_gff_alignment(max_features=5)
        
        self.assertIsNotNone(feature_sequences)
        self.assertGreater(len(feature_sequences), 0)
        self.assertLessEqual(len(feature_sequences), 5)
        
        # Verify sequences are extracted correctly
        for feature, sequence in feature_sequences:
            expected_length = feature.end - feature.start + 1
            self.assertEqual(len(sequence), expected_length)
    
    def test_print_feature_sequences(self):
        """Test printing feature sequences (basic functionality)."""
        aligner = GFFAligner(self.gff_file, self.reference_file)
        aligner.process_gff_alignment()
        
        # This should not raise an exception
        try:
            aligner.print_feature_sequences(max_features=3, max_seq_length=50)
        except Exception as e:
            self.fail(f"print_feature_sequences raised an exception: {e}")


def run_simple_test():
    """Run a simple test that prints sequences from test data."""
    print("Running simple GFF alignment test...")
    
    # Set up logging
    setup_logging(verbose=True)
    
    # Create temporary test data
    test_dir = Path(tempfile.mkdtemp())
    reference_file = test_dir / "reference.fa"
    gff_file = test_dir / "features.gff3"
    
    try:
        # Create test reference
        sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 100  # 4200 bp
        record = SeqRecord(Seq(sequence), id="chr1", description="Test chromosome")
        
        with open(reference_file, 'w') as f:
            SeqIO.write([record], f, "fasta")
        
        # Create test GFF3
        gff_content = """##gff-version 3
##sequence-region chr1 1 4200
chr1	test	gene	500	2000	.	+	.	ID=large_gene;Name=large_test_gene
chr1	test	mRNA	500	2000	.	+	.	ID=large_mRNA;Parent=large_gene;Name=large_transcript
chr1	test	exon	500	800	.	+	.	ID=large_exon1;Parent=large_mRNA;exon_number=1
chr1	test	exon	1200	1500	.	+	.	ID=large_exon2;Parent=large_mRNA;exon_number=2
chr1	test	exon	1800	2000	.	+	.	ID=large_exon3;Parent=large_mRNA;exon_number=3
chr1	test	gene	2500	3000	.	-	.	ID=small_gene;Name=small_test_gene
chr1	test	mRNA	2500	3000	.	-	.	ID=small_mRNA;Parent=small_gene;Name=small_transcript
chr1	test	exon	2500	3000	.	-	.	ID=small_exon1;Parent=small_mRNA;exon_number=1
"""
        
        with open(gff_file, 'w') as f:
            f.write(gff_content)
        
        # Run the alignment test
        aligner = GFFAligner(gff_file, reference_file)
        feature_sequences = aligner.process_gff_alignment()
        
        # Print results
        aligner.print_feature_sequences(max_features=10, max_seq_length=80)
        
        print(f"\nTest completed successfully!")
        print(f"Processed {len(feature_sequences)} features")
        
    finally:
        # Clean up
        shutil.rmtree(test_dir)


if __name__ == '__main__':
    # Run simple test
    run_simple_test()
    
    # Run unit tests
    unittest.main(argv=[''], exit=False, verbosity=2)
