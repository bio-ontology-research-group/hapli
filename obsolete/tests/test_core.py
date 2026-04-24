import pytest
from hapli.core.models import AlignmentResult
from hapli.core.io import GFFProcessor, SequenceExtractor

def test_alignment_result_to_dict():
    res = AlignmentResult(
        feature_id="f1", feature_type="gene", mapq=60,
        target_start=100, target_end=200, identity=0.99, cigar="100M"
    )
    child = AlignmentResult(
        feature_id="f2", feature_type="exon", mapq=60,
        target_start=110, target_end=150, identity=1.0, cigar="40M"
    )
    res.children.append(child)
    
    d = res.to_dict()
    assert d['feature_id'] == "f1"
    assert len(d['children']) == 1
    assert d['children'][0]['feature_id'] == "f2"

def test_gff_processor_load(dummy_gff):
    proc = GFFProcessor(dummy_gff, target_gene="TestGene")
    assert "gene1" in proc.features_by_id
    assert "mrna1" in proc.features_by_id
    assert "exon1" in proc.features_by_id
    
    gene = proc.features_by_id["gene1"]
    children = list(proc.get_children(gene))
    assert len(children) == 1
    assert children[0].id == "mrna1"

def test_gff_processor_target_optimization(dummy_gff):
    # Test that it correctly identifies the target gene and loads its hierarchy
    proc = GFFProcessor(dummy_gff, target_gene="TestGene")
    assert proc.features_by_id['gene1'].featuretype == 'gene'

def test_sequence_extractor(dummy_fasta, dummy_gff):
    proc = GFFProcessor(dummy_gff, target_gene="TestGene")
    ext = SequenceExtractor(dummy_fasta)
    
    # Gene is at 10-90 (1-based)
    # Sequence is "ACGT" repeated...
    
    gene = proc.features_by_id["gene1"]
    seq = ext.get_sequence(gene)
    
    assert len(seq) == 81
    # Check start/end based on "ACGT" pattern
    # index 9 (10th base): 
    # 0123 4567 89
    # ACGT ACGT AC
    # So index 9 is 'C'.
    assert seq.startswith("C")

def test_reverse_complement():
    seq = "ATCG"
    rc = SequenceExtractor._reverse_complement(seq)
    assert rc == "CGAT"
