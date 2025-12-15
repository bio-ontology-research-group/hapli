import pytest
import pysam
from hapli.variation.haplotype import HaplotypeGenerator

def test_generate_haplotypes(dummy_fasta, dummy_vcf, tmp_path):
    # Dummy FASTA: chr1, 100bp (50A, 50T)
    # Dummy VCF: 
    # pos 20: A->C (GT 0|1) -> Het
    # pos 70: T->G (GT 1|1) -> HomAlt
    
    gen = HaplotypeGenerator(dummy_fasta, dummy_vcf)
    out_file = tmp_path / "hap.fa"
    
    # Generate for region chr1:1-100
    gen.generate_haplotypes_for_region("chr1", 0, 100, "sample1", out_file)
    
    assert out_file.exists()
    
    # Check sequences
    fasta = pysam.FastaFile(str(out_file))
    
    # Expected Hap1 (0 from 0|1 at 20, 1 from 1|1 at 70)
    # Pos 20 (1-based) is index 19. Ref is T. GT=0 -> T.
    # Pos 70 (1-based) is index 69. Ref is A. GT=1 -> G.
    h1 = fasta.fetch(fasta.references[0])
    assert h1[19] == 'T' # Ref
    assert h1[69] == 'G' # Alt
    
    # Expected Hap2 (1 from 0|1 at 20, 1 from 1|1 at 70)
    # Pos 20 -> C
    # Pos 70 -> G
    h2 = fasta.fetch(fasta.references[1])
    assert h2[19] == 'C' # Alt
    assert h2[69] == 'G' # Alt
    
    # Check length
    assert len(h1) == 100
