import pytest
import shutil
import json
from pathlib import Path
from hapli.core.io import GFFProcessor, SequenceExtractor
from hapli.alignment.hierarchical import HierarchicalAligner

@pytest.mark.skipif(not shutil.which("minimap2"), reason="minimap2 not found")
def test_hierarchical_alignment(dummy_gff, dummy_fasta, tmp_path):
    # Setup
    # Create a haplotype FASTA that matches the reference exactly
    # so we expect perfect alignment.
    hap_path = tmp_path / "haplotypes.fa"
    shutil.copy(dummy_fasta, hap_path)
    
    # Initialize objects
    gff_proc = GFFProcessor(dummy_gff, target_gene="TestGene")
    seq_ext = SequenceExtractor(dummy_fasta)
    aligner = HierarchicalAligner(gff_proc, seq_ext, threads=1)
    
    gene = gff_proc.features_by_id["gene1"]
    
    # Run alignment
    results = aligner.align_gene(gene, hap_path)
    
    # Check results
    assert "chr1" in results
    res = results["chr1"]
    
    assert res.feature_id == "gene1"
    assert res.identity == 1.0
    assert len(res.children) == 1
    
    mrna = res.children[0]
    assert mrna.feature_id == "mrna1"
    assert len(mrna.children) == 2
    
    # Check Exons
    exons = sorted(mrna.children, key=lambda x: x.target_start)
    assert exons[0].feature_id == "exon1"
    assert exons[1].feature_id == "exon2"
    assert exons[0].identity == 1.0
    assert exons[1].identity == 1.0

@pytest.mark.skipif(not shutil.which("minimap2"), reason="minimap2 not found")
def test_alignment_with_mutation(dummy_gff, dummy_fasta, tmp_path):
    # Create a haplotype with a mutation in exon1
    # Exon1 is 10-40 (indices 9-39).
    # Mutate index 20 (within exon1)
    
    with open(dummy_fasta, 'r') as f:
        lines = f.readlines()
    seq = lines[1].strip()
    
    # Mutate
    mut_seq = list(seq)
    mut_seq[20] = 'C' if mut_seq[20] != 'C' else 'G'
    mut_seq = "".join(mut_seq)
    
    hap_path = tmp_path / "mut_hap.fa"
    with open(hap_path, "w") as f:
        f.write(f">hap1\n{mut_seq}\n")
        
    # Run alignment
    gff_proc = GFFProcessor(dummy_gff, target_gene="TestGene")
    seq_ext = SequenceExtractor(dummy_fasta)
    aligner = HierarchicalAligner(gff_proc, seq_ext, threads=1)
    
    gene = gff_proc.features_by_id["gene1"]
    results = aligner.align_gene(gene, hap_path)
    
    res = results["hap1"]
    # Gene identity should be < 1.0
    # Length is 81. 1 mismatch. Identity approx 80/81 = 0.987
    assert res.identity < 1.0
    assert res.identity > 0.95
    
    mrna = res.children[0]
    # Check exon1
    exon1 = [c for c in mrna.children if c.feature_id == "exon1"][0]
    assert exon1.identity < 1.0
    
    # Exon2 should be perfect
    exon2 = [c for c in mrna.children if c.feature_id == "exon2"][0]
    assert exon2.identity == 1.0
