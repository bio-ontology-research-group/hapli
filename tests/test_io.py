"""
Unit tests for `hapli.core.io` — GFF processor + sequence extractor.

These are exercised transitively by test_paper_cases.py, but direct tests
pin the corner cases that matter:
 - single-gene targeted scan stops after the gene region (doesn't parse the
   whole GFF)
 - Parent=<gene> / Parent=<mRNA> ancestry resolution
 - reverse-complement helper on SequenceExtractor
"""
from __future__ import annotations

from pathlib import Path

import pytest

from hapli.core.io import GFFProcessor, Feature, SequenceExtractor


def _write_gff(path: Path, lines: list[str]) -> Path:
    path.write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
    return path


def _find_gene(proc: GFFProcessor):
    """Find the gene feature in features_by_id (mirrors pipeline._locate_gene)."""
    for f in proc.features_by_id.values():
        if f.featuretype == "gene":
            return f
    return None


def test_gffprocessor_loads_targeted_gene_and_children(tmp_path: Path):
    gff = _write_gff(tmp_path / "anno.gff3", [
        "chr1\tsynth\tgene\t1001\t2000\t.\t+\t.\tID=G1;Name=G1;biotype=protein_coding",
        "chr1\tsynth\tmRNA\t1001\t2000\t.\t+\t.\tID=T1;Parent=G1",
        "chr1\tsynth\texon\t1001\t1200\t.\t+\t.\tID=T1.e1;Parent=T1",
        "chr1\tsynth\texon\t1400\t2000\t.\t+\t.\tID=T1.e2;Parent=T1",
        "chr1\tsynth\tCDS\t1001\t1200\t.\t+\t0\tID=T1.c1;Parent=T1",
        "chr1\tsynth\tCDS\t1400\t2000\t.\t+\t0\tID=T1.c2;Parent=T1",
    ])
    proc = GFFProcessor(gff, target_gene="G1")
    g = _find_gene(proc)
    assert g is not None and g.id == "G1"
    assert g.start == 1001 and g.end == 2000
    assert g.strand == "+"

    children = list(proc.get_children("G1"))
    assert len(children) == 1
    mrna = children[0]
    assert mrna.featuretype == "mRNA" and mrna.id == "T1"

    exons = [c for c in proc.get_children("T1") if c.featuretype == "exon"]
    cds = [c for c in proc.get_children("T1") if c.featuretype == "CDS"]
    assert len(exons) == 2
    assert len(cds) == 2
    assert {e.id for e in exons} == {"T1.e1", "T1.e2"}


def test_gffprocessor_finds_gene_by_name_attribute(tmp_path: Path):
    """User passes a gene Name (common case: BRCA1) not an ID (ENSG…).
    Processor must resolve Name → full feature set."""
    gff = _write_gff(tmp_path / "anno.gff3", [
        "chr1\tsynth\tgene\t1001\t2000\t.\t+\t.\tID=ENSG0001;Name=BRCA1",
        "chr1\tsynth\tmRNA\t1001\t2000\t.\t+\t.\tID=T1;Parent=ENSG0001",
    ])
    proc = GFFProcessor(gff, target_gene="BRCA1")
    g = _find_gene(proc)
    assert g is not None and g.id == "ENSG0001"


def test_gffprocessor_finds_gene_by_gencode_gene_name_attribute(tmp_path: Path):
    """Regression: GENCODE GFF3 uses `gene_name=BRCA1` (lowercase), not `Name=BRCA1`.
    The IBEX run on real GENCODE v45 surfaced this — old code only matched on
    Name= and missed every gene in GENCODE.
    """
    gff = _write_gff(tmp_path / "anno.gff3", [
        "chr17\tHAVANA\tgene\t1001\t2000\t.\t+\t.\tID=ENSG00000012048.25;gene_id=ENSG00000012048.25;gene_type=protein_coding;gene_name=BRCA1",
        "chr17\tHAVANA\tmRNA\t1001\t2000\t.\t+\t.\tID=T1;Parent=ENSG00000012048.25",
    ])
    proc = GFFProcessor(gff, target_gene="BRCA1")
    g = _find_gene(proc)
    assert g is not None and g.id == "ENSG00000012048.25"


def test_gffprocessor_missing_gene_leaves_features_empty(tmp_path: Path):
    gff = _write_gff(tmp_path / "anno.gff3", [
        "chr1\tsynth\tgene\t1001\t2000\t.\t+\t.\tID=G1;Name=G1",
    ])
    proc = GFFProcessor(gff, target_gene="NONEXISTENT")
    assert _find_gene(proc) is None


def test_sequence_extractor_reverse_complement():
    # Static helper; no FASTA required
    rc = SequenceExtractor._reverse_complement
    assert rc("ACGT") == "ACGT"            # palindrome
    assert rc("ATCG") == "CGAT"
    assert rc("AAAANNNNTTTT") == "AAAANNNNTTTT"
    # Lowercase supported
    assert rc("acgt") == "acgt"


def test_feature_dataclass_has_gffutils_compatible_attrs():
    """Feature mimics gffutils.Feature so downstream code treats them identically."""
    f = Feature(
        seqid="chr1", source="synth", featuretype="mRNA",
        start=100, end=200, score=".", strand="+", frame=".",
        attributes={"ID": ["T1"], "Parent": ["G1"]}, id="T1",
    )
    # The attrs pipeline.py reads:
    for attr in ("seqid", "start", "end", "strand", "featuretype", "id", "attributes"):
        assert hasattr(f, attr), f"Feature missing {attr}"
    assert f.id == "T1"
    assert f.attributes["Parent"] == ["G1"]
