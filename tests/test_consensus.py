"""
Regression tests for the Phase 0 consensus-based HaplotypeBuilder.

Protects the two bugs that were in the hand-rolled `_apply_variants`:

1. Symbolic SV alleles (`<DEL>`) were inserted as literal strings into the
   haplotype FASTA instead of being materialised.
2. pysam's `sample_rec.alleles` (genotype-ordered) was indexed with
   `allele_indices` (VCF-ordered), so phased SNVs with GT=1|0 silently
   kept REF on hap1.

Also pins the BND preflight rejection contract (Mode A cannot represent
breakends — must error, not produce garbage).
"""

from __future__ import annotations

import pytest

from hapli.external.consensus import (
    Region,
    UnsupportedVariantError,
    consensus_region,
)


REF_SEQ = (
    # 100 bp, mostly A; unique bases at positions 51 (T) and 71 (C)
    "A" * 50 + "T" + "A" * 19 + "C" + "A" * 29
)


def test_region_parses():
    r = Region.parse("chr1:10-20")
    assert r.chrom == "chr1" and r.start == 10 and r.end == 20
    assert str(r) == "chr1:10-20"


def test_region_rejects_garbage():
    with pytest.raises(ValueError):
        Region.parse("not-a-region")


def test_phased_snv_applies_to_correct_haplotype(tmp_fasta, tmp_vcf):
    """GT=1|0 means hap1 carries ALT, hap2 keeps REF.

    The old hand-rolled code was inverted on this and kept REF on hap1.
    """
    fa = tmp_fasta("ref", {"chr1": REF_SEQ})
    vcf = tmp_vcf(
        "v",
        header_lines=[
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "##contig=<ID=chr1,length=100>",
        ],
        record_lines=[
            # 1-based position 51: REF=T, ALT=G, phase 1|0
            "chr1\t51\t.\tT\tG\t.\tPASS\t.\tGT\t1|0",
        ],
    )
    haps = consensus_region(fa, vcf, sample="S1", region="chr1:1-100")
    # pos 51 → offset 50 in 0-based string
    assert haps["hap1"][50] == "G", "hap1 must carry the ALT on GT=1|0"
    assert haps["hap2"][50] == "T", "hap2 must keep REF on GT=1|0"
    # Rest of sequence unchanged
    assert haps["hap1"][:50] == REF_SEQ[:50]
    assert haps["hap1"][51:] == REF_SEQ[51:]


def test_symbolic_del_materialises_as_deletion_not_literal_string(tmp_fasta, tmp_vcf):
    """A `<DEL>` record with END= must remove bases, not insert the text '<DEL>'.

    This is the Phase 0 headline bug.
    """
    fa = tmp_fasta("ref", {"chr1": REF_SEQ})
    vcf = tmp_vcf(
        "v",
        header_lines=[
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
            '##ALT=<ID=DEL,Description="Deletion">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "##contig=<ID=chr1,length=100>",
        ],
        record_lines=[
            # Symbolic DEL from 10..20 inclusive (11 bp removed), homozygous
            "chr1\t10\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\tGT\t1|1",
        ],
    )
    haps = consensus_region(fa, vcf, sample="S1", region="chr1:1-100")
    for name in ("hap1", "hap2"):
        assert "<DEL>" not in haps[name], f"{name} must not contain the literal string '<DEL>'"
    # VCF symbolic-DEL convention: the anchor base at POS is retained; the
    # deletion covers POS+1..END, so 10 bp are removed (11..20 inclusive).
    assert len(haps["hap1"]) == 100 - 10
    assert len(haps["hap2"]) == 100 - 10


def test_bnd_is_rejected_with_clear_error(tmp_fasta, tmp_vcf):
    """Mode A cannot linearise a translocation. Must error, not silently corrupt."""
    fa = tmp_fasta("ref", {"chr1": REF_SEQ})
    vcf = tmp_vcf(
        "v",
        header_lines=[
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "##contig=<ID=chr1,length=100>",
            "##contig=<ID=chr2,length=100>",
        ],
        record_lines=[
            # Breakend: N joined to chr2:50
            "chr1\t51\t.\tN\tN[chr2:50[\t.\tPASS\tSVTYPE=BND\tGT\t1|0",
        ],
    )
    with pytest.raises(UnsupportedVariantError, match="breakend"):
        consensus_region(fa, vcf, sample="S1", region="chr1:1-100")


def test_symbolic_ins_without_sequence_is_rejected(tmp_fasta, tmp_vcf):
    """A symbolic <INS> with no resolved sequence has nothing for bcftools to insert."""
    fa = tmp_fasta("ref", {"chr1": REF_SEQ})
    vcf = tmp_vcf(
        "v",
        header_lines=[
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
            '##ALT=<ID=INS,Description="Insertion">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "##contig=<ID=chr1,length=100>",
        ],
        record_lines=[
            "chr1\t51\t.\tA\t<INS>\t.\tPASS\tSVTYPE=INS\tGT\t1|0",
        ],
    )
    with pytest.raises(UnsupportedVariantError, match="<INS>"):
        consensus_region(fa, vcf, sample="S1", region="chr1:1-100")


def test_missing_sample_raises(tmp_fasta, tmp_vcf):
    fa = tmp_fasta("ref", {"chr1": REF_SEQ})
    vcf = tmp_vcf(
        "v",
        header_lines=[
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "##contig=<ID=chr1,length=100>",
        ],
        record_lines=[
            "chr1\t51\t.\tT\tG\t.\tPASS\t.\tGT\t1|0",
        ],
    )
    with pytest.raises(ValueError, match="not in VCF"):
        consensus_region(fa, vcf, sample="NOT_PRESENT", region="chr1:1-100")


def test_compound_phased_variants_stay_on_their_haplotypes(tmp_fasta, tmp_vcf):
    """Two variants in the same region with opposite phases must not cross-contaminate.

    Variant A at pos 51: 1|0 (hap1 only)
    Variant B at pos 71: 0|1 (hap2 only)

    A faulty applier that applies both variants to both haplotypes — or drops the
    allele-index mapping on one — will produce wrong sequences.
    """
    fa = tmp_fasta("ref", {"chr1": REF_SEQ})
    vcf = tmp_vcf(
        "v",
        header_lines=[
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "##contig=<ID=chr1,length=100>",
        ],
        record_lines=[
            "chr1\t51\t.\tT\tG\t.\tPASS\t.\tGT\t1|0",
            "chr1\t71\t.\tC\tA\t.\tPASS\t.\tGT\t0|1",
        ],
    )
    haps = consensus_region(fa, vcf, sample="S1", region="chr1:1-100")
    # hap1: only the pos-51 change
    assert haps["hap1"][50] == "G"
    assert haps["hap1"][70] == "C"
    # hap2: only the pos-71 change
    assert haps["hap2"][50] == "T"
    assert haps["hap2"][70] == "A"


def test_symbolic_inversion_is_resolved_and_applied(tmp_fasta, tmp_vcf):
    """bcftools consensus (through ≥1.21) refuses symbolic <INV>. The hapli
    wrapper rewrites INV records into explicit REF/ALT sequence alleles
    (ALT = reverse-complement of the inverted range) before calling bcftools.

    We place a unique 10-bp tag `AACCGGTTAC` at positions 41..50. A symbolic
    <INV> with END=50 on hap1 should reverse-complement positions 41..50
    on hap1 → `GTAACCGGTT`, while hap2 keeps the original tag.
    """
    TAG = "AACCGGTTAC"
    RC_TAG = "GTAACCGGTT"  # reverse-complement of TAG
    # Build a 100 bp reference with the tag embedded at 41..50 (0-based 40..49)
    seq = "A" * 40 + TAG + "A" * 50
    assert len(seq) == 100

    fa = tmp_fasta("invref", {"chr1": seq})
    vcf = tmp_vcf(
        "inv",
        header_lines=[
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End pos">',
            '##ALT=<ID=INV,Description="Inversion">',
            "##contig=<ID=chr1,length=100>",
        ],
        record_lines=[
            # Invert 41..50 inclusive on hap1.
            "chr1\t41\t.\tA\t<INV>\t.\tPASS\tSVTYPE=INV;END=50\tGT\t1|0",
        ],
    )
    haps = consensus_region(fa, vcf, sample="S1", region="chr1:1-100")
    assert haps["hap1"][40:50] == RC_TAG, f"hap1 tag should be reversed-complemented; got {haps['hap1'][40:50]}"
    assert haps["hap2"][40:50] == TAG, "hap2 must keep the original (reference) tag"
    # Flanking sequence unaffected
    assert haps["hap1"][:40] == "A" * 40
    assert haps["hap2"][:40] == "A" * 40


def test_symbolic_tandem_duplication_is_resolved_and_applied(tmp_fasta, tmp_vcf):
    """Symbolic <DUP:TANDEM> must be rewritten into an explicit sequence
    duplication. A unique 10-bp tag duplicated on hap1 should appear twice
    back-to-back on hap1's consensus; hap2 keeps one copy.
    """
    TAG = "AACCGGTTAC"
    seq = "A" * 40 + TAG + "A" * 50
    assert len(seq) == 100

    fa = tmp_fasta("dupref", {"chr1": seq})
    vcf = tmp_vcf(
        "dup",
        header_lines=[
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End pos">',
            '##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">',
            "##contig=<ID=chr1,length=100>",
        ],
        record_lines=[
            # Tandem-duplicate positions 41..50 on hap1.
            "chr1\t41\t.\tA\t<DUP:TANDEM>\t.\tPASS\tSVTYPE=DUP;END=50\tGT\t1|0",
        ],
    )
    haps = consensus_region(fa, vcf, sample="S1", region="chr1:1-110")
    # hap1: 40 A's + TAG + TAG + 50 A's = 110 bases
    assert haps["hap1"][:40] == "A" * 40
    assert haps["hap1"][40:50] == TAG
    assert haps["hap1"][50:60] == TAG, (
        f"expected the tag again as the tandem duplicate; got {haps['hap1'][50:60]!r}"
    )
    # hap2: original 100 bp (but the region string asked for 1-110 — bcftools
    # trims to whatever actually got consensus). Accept length 100 or 110.
    assert TAG in haps["hap2"]
    # Count: hap1 has TAG twice, hap2 has TAG once.
    assert haps["hap1"].count(TAG) == 2
    assert haps["hap2"].count(TAG) == 1
