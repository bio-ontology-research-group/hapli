"""Tests for hapli.core.protein — splicing, translation, reverse-strand.

Includes the **deferred Phase 0 reverse-strand regression test** that guards
the Phase 0 fix to the splicing logic (plus-strand extract + concat + RC the
whole thing, rather than RC each exon individually — the latter is wrong
for multi-exon genes on the minus strand).
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from hapli.core.protein import (
    GffRecord,
    extract_proteins_for_gene,
    load_gff,
    reverse_complement,
    splice_cds,
    translate_cds,
    translate_transcript,
    write_protein_fasta,
)


def _write_fasta(tmp_path: Path, name: str, records: dict[str, str]) -> Path:
    fa = tmp_path / f"{name}.fa"
    with fa.open("w") as f:
        for seqid, seq in records.items():
            f.write(f">{seqid}\n{seq}\n")
    pysam.faidx(str(fa))
    return fa


def _write_gff(tmp_path: Path, name: str, lines: list[str]) -> Path:
    gff = tmp_path / f"{name}.gff3"
    gff.write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
    return gff


def test_reverse_complement_preserves_case_and_n():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("acgt") == "acgt"
    assert reverse_complement("ATGN") == "NCAT"
    assert reverse_complement("nAT") == "ATn"


def test_translate_cds_trims_trailing_partial():
    prot, trailing = translate_cds("ATGGCT" + "A", phase=0)
    assert prot == "MA"
    assert trailing == 1


def test_translate_cds_honours_phase():
    # If phase=1, skip the first base
    prot, trailing = translate_cds("GATGGCTTAA", phase=1)   # skip 'G' → "ATGGCTTAA"
    assert prot == "MA*"
    assert trailing == 0


# ─────────────────────────────────────────────────────────────────────────────
# End-to-end: plus-strand single-exon gene
# ─────────────────────────────────────────────────────────────────────────────
def test_plus_strand_single_exon_translates_cleanly(tmp_path):
    fa = _write_fasta(tmp_path, "ref", {"chr1": "AAAA" + "ATGGCTGCTTAA" + "AAAA"})
    #                                    1234   5678901234
    # ATGGCTGCTTAA spans positions 5..16 (12 bp)
    gff = _write_gff(tmp_path, "ref", [
        "chr1\tsynth\tgene\t5\t16\t.\t+\t.\tID=G1;Name=G1",
        "chr1\tsynth\tmRNA\t5\t16\t.\t+\t.\tID=T1;Parent=G1",
        "chr1\tsynth\texon\t5\t16\t.\t+\t.\tID=T1.exon1;Parent=T1",
        "chr1\tsynth\tCDS\t5\t16\t.\t+\t0\tID=T1.cds1;Parent=T1",
    ])
    proteins = extract_proteins_for_gene(fa, gff, "G1")
    p = proteins["T1"]
    assert p.protein == "MAA*"
    assert p.start_codon_ok
    assert p.stop_codon_ok
    assert p.premature_stop_at is None
    assert p.strand == "+"


# ─────────────────────────────────────────────────────────────────────────────
# REGRESSION: reverse-strand multi-exon gene
#
# This is the Phase 0 deferred regression. The old `SequenceExtractor.get_sequence`
# reverse-complemented *each exon individually*; if we then concatenated the
# exons by genomic-start order and DID NOT reverse the concatenation, we would
# have each exon's bases in reversed order but the exons themselves in 5'→3'
# (forward) transcript order. That mangles the spliced protein.
#
# The correct pattern (enshrined in splice_cds): extract each CDS on the
# forward strand of the haplotype, concatenate in genomic order, then
# reverse-complement the whole concatenated block.
# ─────────────────────────────────────────────────────────────────────────────
def test_reverse_strand_multi_exon_splices_correctly(tmp_path):
    # On the FORWARD strand, place two CDS blocks:
    #   block A at 5..10   = "AGCGCC" (6 bp)
    #   block B at 15..20  = "CATTTA" (6 bp)
    # The minus-strand spliced mRNA should be the RC of block A + block B
    # concatenated in genomic order, i.e.:
    #   concat = "AGCGCC" + "CATTTA" = "AGCGCCCATTTA"
    #   RC(concat) = "TAAATGGGCGCT"
    # That RC'd sequence, read 5'→3', encodes: TAA ATG GGC GCT
    #                                          stop   M   G   A
    # But "TAA" as the first codon would be a stop. We want a real protein; shift
    # the exons so the transcript starts at ATG. Let's design the + strand so that:
    #   spliced + = "AGCGCC" + "CATTTA" = "AGCGCCCATTTA"
    #   RC = "TAAATGGGCGCT"
    # and use phase=3 on the first CDS? simpler: pick a different layout.
    #
    # Layout v2 — directly place "TAAAGCGCCCAT" so that its RC is "ATGGGCGCTTTA".
    # That decodes to "MGA*". Split it across two exons on the + strand.
    #
    #   Forward ref:  AAA TAA AGC GCC CAT AAA  (positions 1..18)
    #                 1234567890123456789
    #                 block A positions 4..9   = "TAAAGC"   (6 bp)
    #                 block B positions 10..15 = "GCCCAT"   (6 bp)
    #   concat (gen order) = "TAAAGC" + "GCCCAT" = "TAAAGCGCCCAT"
    #   RC = "ATGGGCGCTTTA" → translates to "MGA*" (start, G, A, stop)
    ref_seq = "AAA" + "TAAAGC" + "GCCCAT" + "AAA"           # 18 bp
    fa = _write_fasta(tmp_path, "ref", {"chr1": ref_seq})
    gff = _write_gff(tmp_path, "ref", [
        "chr1\tsynth\tgene\t4\t15\t.\t-\t.\tID=G1;Name=G1",
        "chr1\tsynth\tmRNA\t4\t15\t.\t-\t.\tID=T1;Parent=G1",
        "chr1\tsynth\texon\t4\t9\t.\t-\t.\tID=T1.exon1;Parent=T1",
        "chr1\tsynth\tCDS\t4\t9\t.\t-\t0\tID=T1.cds1;Parent=T1",
        "chr1\tsynth\texon\t10\t15\t.\t-\t.\tID=T1.exon2;Parent=T1",
        "chr1\tsynth\tCDS\t10\t15\t.\t-\t0\tID=T1.cds2;Parent=T1",
    ])
    proteins = extract_proteins_for_gene(fa, gff, "G1")
    p = proteins["T1"]
    # CORRECT spliced CDS (reverse-strand): RC("TAAAGC" + "GCCCAT") = "ATGGGCGCTTTA"
    # Translation: ATG GGC GCT TTA → "MGAL"
    assert p.protein == "MGAL", f"reverse-strand splice wrong: {p.protein!r}"
    assert p.start_codon_ok
    # No terminal stop in this construction — that's fine, the test is about splicing.
    assert p.strand == "-"


def test_buggy_per_exon_rc_pattern_produces_wrong_protein(tmp_path):
    """Demonstrate what the OLD (wrong) code would have produced, so that if
    anyone reintroduces the per-exon-RC pattern the test fails loudly.

    Does NOT use splice_cds — does the buggy pattern inline.
    """
    ref_seq = "AAA" + "TAAAGC" + "GCCCAT" + "AAA"
    fa = _write_fasta(tmp_path, "ref", {"chr1": ref_seq})

    # Emulate the old bug: RC each exon individually, then concatenate in
    # genomic order (without RC'ing the concat).
    with pysam.FastaFile(str(fa)) as pf:
        exon_a = pf.fetch("chr1", 3, 9)     # "TAAAGC" (+strand)
        exon_b = pf.fetch("chr1", 9, 15)    # "GCCCAT"
    buggy = reverse_complement(exon_a) + reverse_complement(exon_b)
    # That would be: RC("TAAAGC") + RC("GCCCAT") = "GCTTTA" + "ATGGGC" = "GCTTTAATGGGC"
    # Translated: GCT TTA ATG GGC → "ALMG" — completely different and wrong.
    from Bio.Seq import Seq
    buggy_protein = str(Seq(buggy).translate())
    assert buggy_protein != "MGAL", (
        "Sanity: the buggy per-exon-RC pattern must produce a DIFFERENT protein "
        "than the correct splice_cds output, otherwise the regression test is vacuous."
    )
    assert buggy_protein == "ALMG"


# ─────────────────────────────────────────────────────────────────────────────
# Multi-exon plus strand
# ─────────────────────────────────────────────────────────────────────────────
def test_plus_strand_multi_exon_splices_correctly(tmp_path):
    # Place coding bases so that ATG (exon1) + GCT (exon1) + TAA (exon2) = MA*
    ref_seq = "AAA" + "ATGGCT" + "NN" + "TAA" + "AAA"   # intron represented by NN
    #         1234567890123456789
    #         exon1 positions 4..9  = "ATGGCT"
    #         intron 10..11 = "NN"
    #         exon2 positions 12..14 = "TAA"
    fa = _write_fasta(tmp_path, "ref", {"chr1": ref_seq})
    gff = _write_gff(tmp_path, "ref", [
        "chr1\tsynth\tgene\t4\t14\t.\t+\t.\tID=G1;Name=G1",
        "chr1\tsynth\tmRNA\t4\t14\t.\t+\t.\tID=T1;Parent=G1",
        "chr1\tsynth\texon\t4\t9\t.\t+\t.\tID=T1.exon1;Parent=T1",
        "chr1\tsynth\tCDS\t4\t9\t.\t+\t0\tID=T1.cds1;Parent=T1",
        "chr1\tsynth\texon\t12\t14\t.\t+\t.\tID=T1.exon2;Parent=T1",
        "chr1\tsynth\tCDS\t12\t14\t.\t+\t0\tID=T1.cds2;Parent=T1",
    ])
    proteins = extract_proteins_for_gene(fa, gff, "G1")
    p = proteins["T1"]
    assert p.protein == "MA*"


def test_transcript_without_cds_yields_empty_protein(tmp_path):
    fa = _write_fasta(tmp_path, "ref", {"chr1": "A" * 100})
    gff = _write_gff(tmp_path, "ref", [
        "chr1\tsynth\tgene\t10\t50\t.\t+\t.\tID=G1;Name=G1",
        "chr1\tsynth\tmRNA\t10\t50\t.\t+\t.\tID=T1;Parent=G1",
        "chr1\tsynth\texon\t10\t50\t.\t+\t.\tID=T1.exon1;Parent=T1",
        # no CDS record
    ])
    proteins = extract_proteins_for_gene(fa, gff, "G1")
    p = proteins["T1"]
    assert p.protein == ""
    assert "no CDS records" in p.warnings


def test_write_protein_fasta_handles_multiple_records(tmp_path):
    fa = _write_fasta(tmp_path, "ref", {"chr1": "AAA" + "ATGGCTGCTTAA" + "AAA"})
    gff = _write_gff(tmp_path, "ref", [
        "chr1\tsynth\tgene\t4\t15\t.\t+\t.\tID=G1;Name=G1",
        "chr1\tsynth\tmRNA\t4\t15\t.\t+\t.\tID=T1;Parent=G1",
        "chr1\tsynth\tCDS\t4\t15\t.\t+\t0\tID=T1.cds1;Parent=T1",
        "chr1\tsynth\tmRNA\t4\t15\t.\t+\t.\tID=T2;Parent=G1",
        "chr1\tsynth\tCDS\t4\t15\t.\t+\t0\tID=T2.cds1;Parent=T2",
    ])
    proteins = extract_proteins_for_gene(fa, gff, "G1")
    assert set(proteins) == {"T1", "T2"}
    out = write_protein_fasta(proteins.values(), tmp_path / "prot.fa", header_prefix="hap1.")
    text = out.read_text()
    assert ">hap1.T1" in text
    assert ">hap1.T2" in text


def test_gencode_style_gene_lookup_via_gene_name_attr(tmp_path):
    """Regression: GENCODE GFFs use `gene_name=BRCA1` (lowercase) for the gene
    symbol; reference / Liftoff GFFs use `Name=BRCA1`. transcripts_of_gene
    must match either. The IBEX HG00097/BRCA1 run surfaced this — old code
    only matched on Name=, returning 0 transcripts on every GENCODE gene
    and producing the spurious "No reference protein records" message for
    every (sample, gene) pair.
    """
    fa = _write_fasta(tmp_path, "ref", {"chr17": "AAA" + "ATGGCTGCTTAA" + "AAA"})
    gff = _write_gff(tmp_path, "ref", [
        # GENCODE-style: gene_name= attribute, transcript featuretype
        "chr17\tHAVANA\tgene\t4\t15\t.\t+\t.\tID=ENSG00000012048.25;gene_id=ENSG00000012048.25;gene_type=protein_coding;gene_name=BRCA1",
        "chr17\tHAVANA\ttranscript\t4\t15\t.\t+\t.\tID=ENST00000357654.9;Parent=ENSG00000012048.25;transcript_type=protein_coding",
        "chr17\tHAVANA\tCDS\t4\t15\t.\t+\t0\tID=CDS:ENST00000357654.9:1;Parent=ENST00000357654.9",
    ])
    # User passes the symbol "BRCA1"; lookup must resolve via gene_name.
    proteins = extract_proteins_for_gene(fa, gff, "BRCA1")
    assert "ENST00000357654.9" in proteins, (
        f"GENCODE gene_name lookup failed; found transcripts: {list(proteins)}"
    )
    assert proteins["ENST00000357654.9"].protein == "MAA*"


def test_liftoff_copy_suffix_gene_id_matches(tmp_path):
    """Liftoff renames gene IDs like G1 → G1_0. The extractor must still resolve
    them when the caller asks for the original name.
    """
    fa = _write_fasta(tmp_path, "ref", {"chr1": "AAA" + "ATGGCTGCTTAA" + "AAA"})
    gff = _write_gff(tmp_path, "ref", [
        "chr1\tLiftoff\tgene\t4\t15\t.\t+\t.\tID=G1_0;Name=G1",
        "chr1\tLiftoff\tmRNA\t4\t15\t.\t+\t.\tID=T1;Parent=G1_0",
        "chr1\tLiftoff\tCDS\t4\t15\t.\t+\t0\tID=T1.cds1;Parent=T1",
    ])
    # Caller passes the original gene name "G1" — extractor should find it even
    # though the GFF uses "G1_0".
    proteins = extract_proteins_for_gene(fa, gff, "G1")
    assert "T1" in proteins
    assert proteins["T1"].protein == "MAA*"
