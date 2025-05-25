"""
Unit tests for HGVS formatter and variant models.
"""

import unittest
from scripts.test_data.variant_models import (
    HGVSFormatter, SNV, Insertion, Deletion, Inversion, 
    Duplication, Translocation, ComplexVariant
)


class TestHGVSFormatter(unittest.TestCase):
    """Test HGVS formatter functionality."""
    
    def setUp(self):
        self.formatter = HGVSFormatter()
    
    def test_chromosome_accession_mapping(self):
        """Test chromosome to RefSeq accession mapping."""
        # Test with chr prefix
        self.assertEqual(self.formatter.get_accession('chr1'), 'NC_000001.11')
        self.assertEqual(self.formatter.get_accession('chrX'), 'NC_000023.11')
        self.assertEqual(self.formatter.get_accession('chrM'), 'NC_012920.1')
        
        # Test without chr prefix
        self.assertEqual(self.formatter.get_accession('1'), 'NC_000001.11')
        self.assertEqual(self.formatter.get_accession('X'), 'NC_000023.11')
        self.assertEqual(self.formatter.get_accession('M'), 'NC_012920.1')
    
    def test_genomic_substitution(self):
        """Test SNV/substitution formatting."""
        result = self.formatter.format_genomic_substitution('chr1', 12345, 'A', 'G')
        expected = 'NC_000001.11:g.12345A>G'
        self.assertEqual(result, expected)
    
    def test_genomic_deletion(self):
        """Test deletion formatting."""
        # Single base deletion
        result = self.formatter.format_genomic_deletion('chr1', 12345, 12345)
        expected = 'NC_000001.11:g.12345del'
        self.assertEqual(result, expected)
        
        # Multi-base deletion
        result = self.formatter.format_genomic_deletion('chr1', 12345, 12350)
        expected = 'NC_000001.11:g.12345_12350del'
        self.assertEqual(result, expected)
    
    def test_genomic_insertion(self):
        """Test insertion formatting."""
        result = self.formatter.format_genomic_insertion('chr1', 12345, 12346, 'ATCG')
        expected = 'NC_000001.11:g.12345_12346insATCG'
        self.assertEqual(result, expected)
    
    def test_genomic_duplication(self):
        """Test duplication formatting."""
        # Single base duplication
        result = self.formatter.format_genomic_duplication('chr1', 12345, 12345)
        expected = 'NC_000001.11:g.12345dup'
        self.assertEqual(result, expected)
        
        # Multi-base duplication
        result = self.formatter.format_genomic_duplication('chr1', 12345, 12350)
        expected = 'NC_000001.11:g.12345_12350dup'
        self.assertEqual(result, expected)
    
    def test_genomic_inversion(self):
        """Test inversion formatting."""
        result = self.formatter.format_genomic_inversion('chr1', 12345, 12350)
        expected = 'NC_000001.11:g.12345_12350inv'
        self.assertEqual(result, expected)
    
    def test_genomic_delins(self):
        """Test deletion-insertion formatting."""
        result = self.formatter.format_genomic_delins('chr1', 12345, 12350, 'ATCG')
        expected = 'NC_000001.11:g.12345_12350delinsATCG'
        self.assertEqual(result, expected)
    
    def test_genomic_copy_number(self):
        """Test copy number variant formatting."""
        result = self.formatter.format_genomic_copy_number('chr1', 12345, 12350, 5)
        expected = 'NC_000001.11:g.12345_12350[5]'
        self.assertEqual(result, expected)
    
    def test_hgvs_validation(self):
        """Test HGVS syntax validation."""
        # Valid HGVS strings
        valid_cases = [
            'NC_000001.11:g.12345A>G',
            'NC_000001.11:g.12345del',
            'NC_000001.11:g.12345_12350del',
            'NC_000001.11:g.12345_12346insATCG',
            'NC_000001.11:g.12345dup',
            'NC_000001.11:g.12345_12350dup',
            'NC_000001.11:g.12345_12350inv',
            'NC_000001.11:g.12345delinsATCG',
            'NC_000001.11:g.12345_12350[5]'
        ]
        
        for case in valid_cases:
            with self.subTest(case=case):
                self.assertTrue(self.formatter.validate_hgvs_syntax(case))
        
        # Invalid HGVS strings
        invalid_cases = [
            'invalid_format',
            'NC_000001.11:g.12345',
            'chr1:g.12345A>G',
            'NC_000001.11:12345A>G'
        ]
        
        for case in invalid_cases:
            with self.subTest(case=case):
                self.assertFalse(self.formatter.validate_hgvs_syntax(case))


class TestVariantHGVS(unittest.TestCase):
    """Test HGVS notation generation for variant classes."""
    
    def setUp(self):
        self.formatter = HGVSFormatter()
    
    def test_snv_hgvs(self):
        """Test SNV HGVS notation."""
        snv = SNV(
            chromosome='chr1',
            position=12345,
            ref_allele='A',
            alt_allele='G'
        )
        
        result = snv.to_hgvs_notation(self.formatter)
        expected = 'NC_000001.11:g.12345A>G'
        self.assertEqual(result, expected)
    
    def test_insertion_hgvs(self):
        """Test insertion HGVS notation."""
        insertion = Insertion(
            chromosome='chr1',
            position=12345,
            ref_allele='',
            alt_allele='ATCG',
            insert_sequence='ATCG'
        )
        
        result = insertion.to_hgvs_notation(self.formatter)
        expected = 'NC_000001.11:g.12345_12346insATCG'
        self.assertEqual(result, expected)
    
    def test_deletion_hgvs(self):
        """Test deletion HGVS notation."""
        deletion = Deletion(
            chromosome='chr1',
            position=12345,
            ref_allele='ATCG',
            alt_allele='',
            deleted_length=4
        )
        
        result = deletion.to_hgvs_notation(self.formatter)
        expected = 'NC_000001.11:g.12345_12348del'
        self.assertEqual(result, expected)
    
    def test_duplication_hgvs(self):
        """Test duplication HGVS notation."""
        # Simple duplication (copy number = 2)
        dup = Duplication(
            chromosome='chr1',
            position=12345,
            ref_allele='ATCG',
            alt_allele='ATCGATCG',
            end_position=12348,
            copy_number=2
        )
        
        result = dup.to_hgvs_notation(self.formatter)
        expected = 'NC_000001.11:g.12345_12348dup'
        self.assertEqual(result, expected)
        
        # Higher copy number
        dup_high = Duplication(
            chromosome='chr1',
            position=12345,
            ref_allele='ATCG',
            alt_allele='ATCGATCGATCGATCG',
            end_position=12348,
            copy_number=4
        )
        
        result = dup_high.to_hgvs_notation(self.formatter)
        expected = 'NC_000001.11:g.12345_12348[4]'
        self.assertEqual(result, expected)
    
    def test_inversion_hgvs(self):
        """Test inversion HGVS notation."""
        inversion = Inversion(
            chromosome='chr1',
            position=12345,
            ref_allele='ATCG',
            alt_allele='CGAT',
            end_position=12348
        )
        
        result = inversion.to_hgvs_notation(self.formatter)
        expected = 'NC_000001.11:g.12345_12348inv'
        self.assertEqual(result, expected)
    
    def test_complex_variant_hgvs(self):
        """Test complex variant (delins) HGVS notation."""
        complex_var = ComplexVariant(
            chromosome='chr1',
            position=12345,
            ref_allele='ATCG',
            alt_allele='GGCC',
            deleted_length=4,
            insert_sequence='GGCC'
        )
        
        result = complex_var.to_hgvs_notation(self.formatter)
        expected = 'NC_000001.11:g.12345_12348delinsGGCC'
        self.assertEqual(result, expected)
    
    def test_translocation_hgvs(self):
        """Test translocation HGVS notation."""
        translocation = Translocation(
            chromosome='chr1',
            position=12345,
            ref_allele='A',
            alt_allele='A',
            chr2='chr2',
            pos2=67890
        )
        
        result = translocation.to_hgvs_notation(self.formatter)
        self.assertIn('Complex translocation', result)
        self.assertIn('NC_000001.11:g.12345', result)
        self.assertIn('NC_000002.12:g.67890', result)


if __name__ == '__main__':
    unittest.main()
