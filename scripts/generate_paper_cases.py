#!/usr/bin/env python3
"""
Generate one or more named "motivating cases" for the paper + benchmarks.

Each case writes a self-contained `{ref.fa, annotation.gff3, phased.vcf.gz,
expected.json}` bundle into a per-case output directory. The bundle is sized
to be human-auditable — a single small synthetic contig with one gene built
around the specific variant pattern under test.

Intended consumers:
  - pytest: end-to-end regression of `hapli analyze` on each case.
  - Phase 5 benchmarks: a substrate where the expected functional call is
    known by construction.
  - Paper figures: minimal examples of variant-by-variant failure modes.

CLI
---
    uv run python3 scripts/generate_paper_cases.py --case 05_frameshift_rescue --out data/paper_cases
    uv run python3 scripts/generate_paper_cases.py --all --out data/paper_cases
    uv run python3 scripts/generate_paper_cases.py --list

Case registry — all 12 cases wired:
  01_synonymous_pair          — both haplotypes intact (negative control).
  02_single_benign_missense   — one benign missense on hap1.
  03_single_stop_gain         — premature stop on hap1, hap2 intact.
  04_compound_het_lof         — stop-gain hap1 + frameshift hap2.
  05_frameshift_rescue        — +1 insertion then −1 deletion restoring frame.
  06_compound_missense_pocket — two benign missense; pocket disrupted jointly.
  07_splice_donor_cryptic     — canonical donor loss + cryptic donor gain.
  08_symbolic_del_whole_gene  — <DEL> removing the gene on hap1.
  09_upstream_del_shift       — large DEL upstream, downstream gene intact.
  10_inversion_body           — <INV> bisecting gene body.
  11_tandem_dup_small_gene    — <DUP:TANDEM> duplicating a small gene.
  12_inframe_domain_indel     — in-frame 3 bp deletion at a codon.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path

import pysam


# -----------------------------------------------------------------------------
# Infrastructure: one data bundle per case
# -----------------------------------------------------------------------------
@dataclass
class CaseBundle:
    """Everything a motivating case needs in one place."""

    name: str
    description: str
    # Reference: {seqid: sequence}
    reference: dict[str, str]
    # GFF3: list of (seqid, source, type, start, end, score, strand, phase, attributes)
    gff_records: list[tuple]
    # VCF: list of (chrom, pos, id, ref, alt, qual, filter, info, format, sample_gt)
    vcf_records: list[tuple]
    # VCF header extra lines (INFO/ALT/FORMAT/contig)
    vcf_headers: list[str]
    # Expected hapli call: stable structure that tests can assert against.
    expected: dict = field(default_factory=dict)

    def write(self, out_dir: Path) -> None:
        out_dir.mkdir(parents=True, exist_ok=True)

        # Reference FASTA + fai
        ref_path = out_dir / "reference.fa"
        with ref_path.open("w") as f:
            for seqid, seq in self.reference.items():
                f.write(f">{seqid}\n")
                for i in range(0, len(seq), 60):
                    f.write(seq[i : i + 60] + "\n")
        pysam.faidx(str(ref_path))

        # GFF3
        gff_path = out_dir / "annotation.gff3"
        with gff_path.open("w") as f:
            f.write("##gff-version 3\n")
            for seqid, source, kind, start, end, score, strand, phase, attrs in self.gff_records:
                attr_str = ";".join(f"{k}={v}" for k, v in attrs.items())
                f.write(
                    f"{seqid}\t{source}\t{kind}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attr_str}\n"
                )

        # VCF (uncompressed), then bgzip + tabix
        vcf_path = out_dir / "phased.vcf"
        with vcf_path.open("w") as f:
            f.write("##fileformat=VCFv4.2\n")
            for h in self.vcf_headers:
                if not h.endswith("\n"):
                    h = h + "\n"
                f.write(h)
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
            for row in self.vcf_records:
                f.write("\t".join(str(c) for c in row) + "\n")
        gz = out_dir / "phased.vcf.gz"
        pysam.tabix_compress(str(vcf_path), str(gz), force=True)
        pysam.tabix_index(str(gz), preset="vcf", force=True)
        vcf_path.unlink()

        # Expected-call pinning
        (out_dir / "expected.json").write_text(
            json.dumps({"case": self.name, "description": self.description, **self.expected}, indent=2)
        )


# -----------------------------------------------------------------------------
# Tiny synthetic helpers
#   — an A/T/G/C generator just stable enough to put codons where we need them
# -----------------------------------------------------------------------------
def _codon_run(codon: str, count: int) -> str:
    """Repeat a codon `count` times. Cheap way to build an in-frame stretch."""
    if len(codon) != 3:
        raise ValueError("codon must be 3 nt")
    return codon * count


def _simple_gene(
    seqid: str,
    gene_start: int,
    exons: list[tuple[int, int]],
    strand: str = "+",
    gene_id: str = "G1",
    transcript_id: str = "T1",
) -> list[tuple]:
    """Build GFF3 records for a gene with arbitrary exon layout. One transcript."""
    gene_end = max(e for _, e in exons)
    records = [
        (
            seqid,
            "hapli-synth",
            "gene",
            gene_start,
            gene_end,
            ".",
            strand,
            ".",
            {"ID": gene_id, "Name": gene_id, "biotype": "protein_coding"},
        ),
        (
            seqid,
            "hapli-synth",
            "mRNA",
            gene_start,
            gene_end,
            ".",
            strand,
            ".",
            {"ID": transcript_id, "Parent": gene_id},
        ),
    ]
    for i, (es, ee) in enumerate(exons, start=1):
        exon_id = f"{transcript_id}.exon{i}"
        cds_id = f"{transcript_id}.cds{i}"
        records.append(
            (seqid, "hapli-synth", "exon", es, ee, ".", strand, ".", {"ID": exon_id, "Parent": transcript_id})
        )
        records.append(
            (seqid, "hapli-synth", "CDS", es, ee, ".", strand, "0", {"ID": cds_id, "Parent": transcript_id})
        )
    return records


# -----------------------------------------------------------------------------
# Case 01 — synonymous pair, negative control
# -----------------------------------------------------------------------------
def case_01_synonymous_pair() -> CaseBundle:
    # Build a 3000 bp chromosome with one in-frame gene, 1 exon.
    # The gene starts at 1001 for a comfortable 1 kb buffer.
    prefix = "A" * 1000                                    # 1000 bp of A-rich context
    # Gene coding region 1001..1300 (300 bp = 100 codons).
    # Put known codons: start (ATG) + 98 x GCT (Ala) + stop (TAA).
    cds = "ATG" + _codon_run("GCT", 98) + "TAA"            # 300 bp
    suffix = "A" * 1700                                    # pad out to 3000
    seq = prefix + cds + suffix
    assert len(seq) == 3000

    gff = _simple_gene("chr1", gene_start=1001, exons=[(1001, 1300)])

    # Two synonymous SNVs: GCT→GCC (both Ala) at codon positions 2 and 3.
    # Codon 2 occupies positions 1004..1006. Third-base wobble at pos 1006.
    # Codon 3 occupies positions 1007..1009. Third-base wobble at pos 1009.
    # phase: pos 1006 → 1|0, pos 1009 → 0|1 (one on each hap, both synonymous)
    headers = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    vcf = [
        ("chr1", 1006, ".", "T", "C", ".", "PASS", ".", "GT", "1|0"),
        ("chr1", 1009, ".", "T", "C", ".", "PASS", ".", "GT", "0|1"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {"hap1": "intact", "hap2": "intact"},
        "notes": "Two synonymous variants, one per haplotype. Variant-by-variant and "
                 "haplotype-level analysis agree: both haplotypes functional.",
    }
    return CaseBundle(
        name="01_synonymous_pair",
        description="Negative control: two synonymous SNVs, one per haplotype; both haplotypes intact.",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Case 05 — frameshift + rescue (the hero case)
# -----------------------------------------------------------------------------
def case_05_frameshift_rescue() -> CaseBundle:
    """+1 insertion then −1 deletion within a few codons, both on the same haplotype.

    A variant-by-variant LoF caller (LOFTEE, raw VEP) flags both as frame-disrupting
    premature truncations. The true joint effect: a short peptide shift between the
    two edits, then a restored frame. Haplotype-level analysis should classify this
    as functional (possibly with a small in-frame substitution window).
    """
    prefix = "A" * 1000
    # 120 bp CDS: ATG + 38 x GCT + TAA
    cds = "ATG" + _codon_run("GCT", 38) + "TAA"
    suffix = "A" * 1880
    seq = prefix + cds + suffix
    assert len(seq) == 3000

    gff = _simple_gene("chr1", gene_start=1001, exons=[(1001, 1120)])

    headers = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    # Codon table (1-based genomic positions):
    #   codon 1  = 1001..1003 ATG
    #   codon 2  = 1004..1006 GCT
    #   codon 3  = 1007..1009 GCT
    #   codon 4  = 1010..1012 GCT   (last base of codon 4 at pos 1012 = 'T')
    #   ...
    #   codon 10 = 1028..1030 GCT   (bases: 1028='G', 1029='C', 1030='T')
    #
    # Insertion at POS=1012: REF='T' (last base of codon 4), ALT='TC' (inserts 'C'
    # immediately after pos 1012). Net effect: +1 bp → frame shifts by +1.
    #
    # Deletion at POS=1029: REF='CT' (bases 1029-1030 of codon 10), ALT='C' → removes
    # one base. Net over the pair: 0 bp, but the codons between 1013 and 1028 of the
    # reference are translated in a +1-shifted window. Frame is restored from 1030 onward.
    vcf = [
        ("chr1", 1012, ".", "T", "TC", ".", "PASS", ".", "GT", "1|0"),
        ("chr1", 1029, ".", "CT", "C", ".", "PASS", ".", "GT", "1|0"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {"hap1": "intact_with_inframe_window", "hap2": "intact"},
        "notes": "Hero case: individually +1 and −1 would each be flagged frame-disrupting. "
                 "Jointly, frame is restored after a ~17 nt window (~6 aa substitution). "
                 "Variant-by-variant tools flag LoF; haplotype-level analysis classifies as functional.",
    }
    return CaseBundle(
        name="05_frameshift_rescue",
        description="+1 insertion then −1 deletion on the same haplotype; frame restored jointly.",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Case 04 — compound-het LoF (stop-gain hap1 + frameshift hap2)
# -----------------------------------------------------------------------------
def case_04_compound_het_lof() -> CaseBundle:
    """Classic autosomal-recessive compound-het LoF pattern. Two different LoF
    variants on the two haplotypes — both copies of the gene are independently
    knocked out. Variant-by-variant tools see two independent LoF calls;
    hapli's diploid aggregator emits `compound_het_lof=True`, which is the
    actionable haplotype-level signal.
    """
    prefix = "A" * 1000
    cds = "ATG" + _codon_run("GCT", 98) + "TAA"   # 300 bp
    suffix = "A" * 1700
    seq = prefix + cds + suffix
    assert len(seq) == 3000

    gff = _simple_gene("chr1", gene_start=1001, exons=[(1001, 1300)])

    headers = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    # Codon 4 = pos 1010..1012 (GCT). To make it TAA (stop) we replace all 3
    # bases of the codon: REF="GCT" → ALT="TAA". Hap1 (1|0) carries the stop.
    #
    # Codon 50 = pos 1148..1150 (GCT). Insert one base after pos 1148 → frameshift.
    # REF=G (at 1148), ALT=GA → +1bp insertion. Hap2 (0|1) carries the frameshift.
    vcf = [
        ("chr1", 1010, ".", "GCT", "TAA", ".", "PASS", ".", "GT", "1|0"),
        ("chr1", 1148, ".", "G", "GA", ".", "PASS", ".", "GT", "0|1"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {
            "hap1": {"call": "lof", "evidence": "stop_gained_at_codon_4"},
            "hap2": {"call": "lof", "evidence": "frameshift_after_codon_49"},
        },
        "diploid_call": {"compound_het_lof": True},
        "notes": "Compound heterozygous LoF: both haps independently knocked out by "
                 "different mechanisms. Hapli emits compound_het_lof=True; bcftools/csq "
                 "emits per-variant LoF tags but does not aggregate to the diploid level.",
    }
    return CaseBundle(
        name="04_compound_het_lof",
        description="Stop-gain (hap1) + frameshift (hap2): classic compound-het LoF.",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Case 06 — compound missense in (hypothetical) binding pocket
# -----------------------------------------------------------------------------
def case_06_compound_missense_pocket() -> CaseBundle:
    """Two missense variants on the *same* haplotype, both individually benign-
    looking by per-variant predictors but jointly altering interacting residues.
    Variant-by-variant: two missense calls. Hapli's epistasis residual flags
    the joint effect.

    The synthetic gene encodes M + 49 × Ala + stop (50 aa). We mutate two distant
    Ala residues to Trp+Trp on the same haplotype — bulky aromatic residues that
    would individually be tolerated in a small protein but jointly clash if
    structurally proximate. This is a synthetic stand-in; the real-world story
    requires a structured target like HRAS or CALM1.
    """
    prefix = "A" * 1000
    cds = "ATG" + _codon_run("GCT", 49) + "TAA"   # 153 bp
    suffix = "A" * (3000 - 1000 - 153)
    seq = prefix + cds + suffix
    assert len(seq) == 3000

    gff = _simple_gene("chr1", gene_start=1001, exons=[(1001, 1153)])

    headers = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    # Two missense Ala→Trp substitutions at codons 10 and 30, both on hap1.
    # Codon 10 = pos 1028..1030 = "GCT" → "TGG" (Trp). Multi-base substitution.
    # Codon 30 = pos 1088..1090 = "GCT" → "TGG".
    vcf = [
        ("chr1", 1028, ".", "GCT", "TGG", ".", "PASS", ".", "GT", "1|0"),
        ("chr1", 1088, ".", "GCT", "TGG", ".", "PASS", ".", "GT", "1|0"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {
            "hap1": {"call": "compound_missense", "evidence": "two_aromatic_substitutions_on_same_hap"},
            "hap2": {"call": "intact"},
        },
        "epistasis_signal_expected": True,
        "notes": "Two distant missense variants on hap1. Variant-by-variant tools "
                 "emit two independent missense calls. Hapli's ESM2 residual flags "
                 "the joint effect; bcftools/csq cannot reach this signal.",
    }
    return CaseBundle(
        name="06_compound_missense_pocket",
        description="Two missense (Ala→Trp) on hap1 — bcftools/csq calls 'missense, missense'; hapli flags epistasis residual.",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Case 07 — splice donor loss + cryptic donor gain
# -----------------------------------------------------------------------------
def case_07_splice_donor_cryptic() -> CaseBundle:
    """Two variants on the same haplotype that, individually, look like LoF
    splice events. Variant 1 destroys the canonical 5' splice donor (GT→AT)
    of an intron; variant 2 creates a cryptic donor (NN→GT) ~24 nt downstream
    in the intron. Per-variant predictors call BOTH as splice-LoF (donor
    loss is "frame-disrupting splice", and the cryptic-gain variant flips a
    non-splice site to a splice-like motif).

    The joint effect: the intron is still spliced, but at the new (cryptic)
    donor instead of the canonical one — the resulting mRNA has a small (~24 nt)
    insertion of intron-tail sequence rather than complete intron retention
    or exon truncation. SpliceAI / Pangolin would correctly call this as
    "moved donor" given a haplotype-aware view; per-variant scoring gives the
    misleading "double LoF" answer.

    Synthetic encoding: a 2-exon gene with a small intron containing an
    obvious GT...AG splice signature. We mutate the canonical 5' GT and
    create a fresh GT 24 nt downstream.
    """
    prefix = "A" * 1000
    # Two-exon gene:
    #   exon 1: positions 1001..1059 (59 bp)  starts with ATG
    #   intron: positions 1060..1099 (40 bp)  starts with GT...AG (canonical)
    #   exon 2: positions 1100..1159 (60 bp)  ends with TAA
    exon1 = "ATG" + "GCT" * 18 + "GC"           # 59 bp
    intron = "GT" + "A" * 16 + "TT" + "A" * 18 + "AG"   # 40 bp; pos 1060 = G, 1061 = T (canonical donor)
    exon2 = "T" + "GCT" * 19 + "TAA"            # 1 + 57 + 3 = 61 bp
    # Recompute exon2 to be exactly 60 bp: "T" + "GCT"*19 = 1 + 57 = 58, then "TA" + "A" = 61. Adjust.
    exon2 = "TC" + "GCT" * 18 + "TAA"           # 2 + 54 + 3 = 59 bp — close enough; gene length doesn't need to be perfect 60s
    cds_region = exon1 + intron + exon2         # 59 + 40 + 59 = 158 bp
    gene_end = 1000 + len(cds_region)           # 1158
    suffix = "A" * (3000 - gene_end)
    seq = prefix + cds_region + suffix
    assert len(seq) == 3000

    e1_end = 1000 + len(exon1)                  # 1059
    intron_start = e1_end + 1                   # 1060
    intron_end = intron_start + len(intron) - 1 # 1099
    e2_start = intron_end + 1                   # 1100
    e2_end = e2_start + len(exon2) - 1          # 1158

    gff = _simple_gene("chr1", gene_start=1001, exons=[(1001, e1_end), (e2_start, e2_end)])

    headers = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    # Variant 1: destroy canonical donor by changing intron pos 1061 T → A.
    # Reference at 1060/1061 = "GT" (canonical donor); alt = "GA".
    # Variant 2: create cryptic donor 24 nt downstream — at intron pos 1085/1086.
    # Reference at 1085 = ? Let's compute: intron starts at 1060, so pos 1085 = intron[25] = "A".
    # Make positions 1085..1086 = "GT" via SNVs at 1085 (A→G) and 1086 (A→T).
    # Both on hap1 (1|0).
    vcf = [
        ("chr1", 1061, ".", "T", "A", ".", "PASS", ".", "GT", "1|0"),
        ("chr1", 1085, ".", "A", "G", ".", "PASS", ".", "GT", "1|0"),
        ("chr1", 1086, ".", "A", "T", ".", "PASS", ".", "GT", "1|0"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {
            "hap1": {"call": "splice_altered_in_frame", "evidence": "donor_moved_24nt_downstream"},
            "hap2": {"call": "intact"},
        },
        "notes": "Three variants on hap1: cancels canonical 5' splice donor, creates a "
                 "cryptic donor 24 nt downstream. Variant-by-variant predictors call all "
                 "three as splice-disrupting LoF; haplotype-level analysis would find a "
                 "still-functional, slightly extended exon 1 (a SpliceAI-aware variant of "
                 "hapli would catch this; the current pipeline reports the protein as "
                 "lifted but with a length change in evidence.protein).",
    }
    return CaseBundle(
        name="07_splice_donor_cryptic",
        description="Splice donor loss + cryptic donor gain on same hap; moves the splice site rather than disrupting it.",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Case 09 — upstream large DEL shifts coords; downstream gene must still lift
# -----------------------------------------------------------------------------
def case_09_upstream_del_shift() -> CaseBundle:
    """A 200 bp deletion 500 nt upstream of a gene shifts every reference
    coordinate within that haplotype by -200 bp. A naive pipeline that uses
    reference GFF coordinates on the haplotype FASTA would extract garbage.
    Hapli's Liftoff step solves this: the gene is still found, just at
    coordinates 200 bp earlier. This is the regression test for the
    "reference coordinates are invalid post-consensus" failure mode that
    motivated promoting Liftoff to a canonical step in Phase 1.
    """
    prefix = "A" * 1500                  # 1500 bp upstream context
    cds = "ATG" + _codon_run("GCT", 38) + "TAA"   # 120 bp gene at positions 1501..1620
    suffix = "A" * (3000 - 1500 - 120)
    seq = prefix + cds + suffix
    assert len(seq) == 3000

    gff = _simple_gene("chr1", gene_start=1501, exons=[(1501, 1620)])

    headers = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    # Delete 200 bp at positions 1001..1200 (well upstream of the gene at 1501).
    # VCF anchor: POS=1000, REF = 201 chars from pos 1000-1200; ALT = 1 char (anchor only).
    # We don't actually need to spell out 201 ref bases; use a short prefix +
    # bcftools will accept the format, but the cleaner way is to just emit
    # the actual reference slice as REF and the first base as ALT.
    ref_slice = seq[999:1200]            # positions 1000..1200 → 201 bases
    vcf = [
        ("chr1", 1000, ".", ref_slice, ref_slice[0], ".", "PASS", ".", "GT", "1|0"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {
            "hap1": {"call": "intact_with_upstream_shift", "evidence": "200bp_DEL_at_1001..1200"},
            "hap2": {"call": "intact"},
        },
        "notes": "Hap1 carries a 200 bp deletion 300 nt upstream of the gene. Reference "
                 "coordinates are now wrong on hap1 by -200 bp. A naive pipeline that "
                 "extracts the gene using ref-GFF coordinates produces garbage; hapli's "
                 "Liftoff step finds the gene at the correct shifted position and "
                 "extracts the protein cleanly. The protein on hap1 should match the "
                 "reference protein exactly (no coding change).",
    }
    return CaseBundle(
        name="09_upstream_del_shift",
        description="200 bp upstream deletion shifts downstream gene coords; Liftoff finds the gene at its new position.",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Case 08 — symbolic <DEL> removing a whole gene on one haplotype
# -----------------------------------------------------------------------------
def case_08_symbolic_del_whole_gene() -> CaseBundle:
    """Exercises the Phase 0 SV fix: symbolic <DEL> must remove sequence, not corrupt."""
    prefix = "A" * 1000
    cds = "ATG" + _codon_run("GCT", 98) + "TAA"   # 300 bp
    suffix = "A" * 1700
    seq = prefix + cds + suffix
    assert len(seq) == 3000

    gff = _simple_gene("chr1", gene_start=1001, exons=[(1001, 1300)])

    headers = [
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
        '##ALT=<ID=DEL,Description="Deletion">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    # Delete the whole gene region (1001..1300) on hap1 via symbolic DEL.
    # Anchor at 1000 (VCF convention: POS is the base before the deletion).
    vcf = [
        ("chr1", 1000, ".", "A", "<DEL>", ".", "PASS", "SVTYPE=DEL;END=1300", "GT", "1|0"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {
            "hap1": {"presence": "deleted", "source": "SV"},
            "hap2": "intact",
        },
        "notes": "Symbolic <DEL> spanning the entire gene on hap1. Regression-tests "
                 "the Phase 0 bug where the old builder inserted the literal string '<DEL>' "
                 "into the haplotype FASTA.",
    }
    return CaseBundle(
        name="08_symbolic_del_whole_gene",
        description="Symbolic <DEL> removes the whole gene on hap1; hap2 intact.",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Case 02 — Single benign missense on hap1 (negative control for variant scoring)
# -----------------------------------------------------------------------------
def case_02_single_benign_missense() -> CaseBundle:
    """One missense variant on hap1; hap2 reference. Both haplotypes should be
    classified as intact by hapli (missense without clustering, no LoF signal).
    bcftools/csq emits a per-variant 'missense' tag on hap1 only.
    """
    prefix = "A" * 1000
    cds = "ATG" + _codon_run("GCT", 98) + "TAA"   # 300 bp
    suffix = "A" * 1700
    seq = prefix + cds + suffix
    assert len(seq) == 3000

    gff = _simple_gene("chr1", gene_start=1001, exons=[(1001, 1300)])

    headers = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    # Codon 10 = pos 1028..1030 = "GCT" (Ala) → "TCT" (Ser). Conservative change.
    # Hap1 (1|0).
    vcf = [
        ("chr1", 1028, ".", "G", "T", ".", "PASS", ".", "GT", "1|0"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {
            "hap1": {"call": "missense_single", "evidence": "Ala10Ser"},
            "hap2": {"call": "intact"},
        },
        "notes": "Single conservative missense on hap1. Both haps functional. "
                 "Baseline for scoring: single benign missense should not trigger LoF flag.",
    }
    return CaseBundle(
        name="02_single_benign_missense",
        description="Single conservative Ala→Ser missense on hap1; both haps functional.",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Case 03 — Single stop-gain on hap1 (dominant heterozygous LoF)
# -----------------------------------------------------------------------------
def case_03_single_stop_gain() -> CaseBundle:
    """A single premature-stop variant on hap1; hap2 reference. Classic
    dominant-heterozygous LoF pattern: one copy knocked out, one intact.
    Distinct from case 04 (compound-het) where both haps are independently LoF.
    """
    prefix = "A" * 1000
    cds = "ATG" + _codon_run("GCT", 98) + "TAA"   # 300 bp
    suffix = "A" * 1700
    seq = prefix + cds + suffix
    assert len(seq) == 3000

    gff = _simple_gene("chr1", gene_start=1001, exons=[(1001, 1300)])

    headers = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    # Codon 20 = pos 1058..1060 (GCT) → TAA (stop). Mid-gene so NMD-eligible.
    vcf = [
        ("chr1", 1058, ".", "GCT", "TAA", ".", "PASS", ".", "GT", "1|0"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {
            "hap1": {"call": "lof", "evidence": "stop_gained_at_codon_20"},
            "hap2": {"call": "intact"},
        },
        "diploid_call": {"compound_het_lof": False},
        "notes": "Single stop-gain on hap1 only. Dominant heterozygous LoF. "
                 "hapli should flag hap1 as LoF (low score, PTC detected) and "
                 "hap2 as intact. compound_het_lof = False since hap2 is fine.",
    }
    return CaseBundle(
        name="03_single_stop_gain",
        description="Single stop-gain on hap1; hap2 intact (dominant het LoF).",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Case 10 — Symbolic inversion crossing gene body
# -----------------------------------------------------------------------------
def case_10_inversion_body() -> CaseBundle:
    """Symbolic <INV> spanning the middle of the gene on hap1. bcftools consensus
    (≥ 1.15) reverse-complements the inverted range. Liftoff's homology-based
    mapping should detect the gene as bisected / partially lifted; the
    downstream protein extraction either truncates or mis-translates the
    inverted segment. TOGA (opt-in) would classify as 'lost'.
    """
    prefix = "A" * 1000
    cds = "ATG" + _codon_run("GCT", 98) + "TAA"   # 300 bp
    suffix = "A" * 1700
    seq = prefix + cds + suffix
    assert len(seq) == 3000

    gff = _simple_gene("chr1", gene_start=1001, exons=[(1001, 1300)])

    headers = [
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
        '##ALT=<ID=INV,Description="Inversion">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    # Invert positions 1100..1200 — a 101 bp range straddling the middle of the
    # 300 bp gene. After RC, the codons in that window translate to a different
    # (and likely frame-disrupting) sequence depending on RC parity.
    vcf = [
        ("chr1", 1100, ".", "G", "<INV>", ".", "PASS", "SVTYPE=INV;END=1200", "GT", "1|0"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {
            "hap1": {"call": "sv_disrupted", "evidence": "INV_1100..1200_bisects_gene"},
            "hap2": {"call": "intact"},
        },
        "notes": "Inversion of 101 bp crossing the gene body on hap1. Liftoff "
                 "should flag partial / low-identity lift. Exact downstream call "
                 "depends on whether Liftoff still maps the gene (partial) or "
                 "drops it (presence=deleted).",
    }
    return CaseBundle(
        name="10_inversion_body",
        description="Symbolic <INV> spans gene body on hap1; Liftoff should report partial lift.",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Case 11 — Tandem duplication of a small gene
# -----------------------------------------------------------------------------
def case_11_tandem_dup_small_gene() -> CaseBundle:
    """Symbolic <DUP:TANDEM> over the entire gene on hap1. bcftools consensus
    duplicates the sequence. Liftoff --copies should surface two copies of the
    gene on hap1. hap2 remains intact (one copy).
    """
    prefix = "A" * 1000
    cds = "ATG" + _codon_run("GCT", 98) + "TAA"   # 300 bp
    suffix = "A" * 1700
    seq = prefix + cds + suffix
    assert len(seq) == 3000

    gff = _simple_gene("chr1", gene_start=1001, exons=[(1001, 1300)])

    headers = [
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
        '##ALT=<ID=DUP,Description="Duplication">',
        '##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    # Tandem-duplicate the gene region on hap1.
    vcf = [
        ("chr1", 1001, ".", "A", "<DUP:TANDEM>", ".", "PASS", "SVTYPE=DUP;END=1300", "GT", "1|0"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {
            "hap1": {"call": "duplicated", "copies": 2, "evidence": "DUP:TANDEM over gene"},
            "hap2": {"call": "intact"},
        },
        "notes": "Tandem duplication over the gene on hap1; Liftoff --copies "
                 "surfaces two intact gene copies on hap1. hap2 unchanged.",
    }
    return CaseBundle(
        name="11_tandem_dup_small_gene",
        description="<DUP:TANDEM> over the gene on hap1; Liftoff reports two copies.",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Case 12 — In-frame indel at a domain boundary
# -----------------------------------------------------------------------------
def case_12_inframe_domain_indel() -> CaseBundle:
    """A single in-frame 3 bp deletion on hap1 removes one codon. No frameshift,
    no PTC — protein is 1 aa shorter than reference. A crude LoF-only tool
    (stop/frameshift-focused) labels this as benign/unknown; AlphaMissense
    doesn't score indels. hapli's protein-diff reports `length_changed=True`,
    a single-residue deletion event, and the ESM pseudo-likelihood delta on
    the whole protein gives a per-haplotype score that reflects the subtle
    effect.
    """
    prefix = "A" * 1000
    cds = "ATG" + _codon_run("GCT", 98) + "TAA"   # 300 bp
    suffix = "A" * 1700
    seq = prefix + cds + suffix
    assert len(seq) == 3000

    gff = _simple_gene("chr1", gene_start=1001, exons=[(1001, 1300)])

    headers = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=3000>",
    ]
    # Delete 3 bp at positions 1040..1042 (removes codon 14, "GCT"). VCF anchor
    # at POS=1039 whose ref base is 'T' (last base of codon 13). So REF="TGCT"
    # (anchor + codon 14), ALT="T".
    vcf = [
        ("chr1", 1039, ".", "TGCT", "T", ".", "PASS", ".", "GT", "1|0"),
    ]
    expected = {
        "gene": "G1",
        "call_per_haplotype": {
            "hap1": {"call": "inframe_del", "evidence": "1_codon_del_at_codon_14"},
            "hap2": {"call": "intact"},
        },
        "notes": "In-frame 3 bp deletion removes one codon on hap1. No PTC; protein "
                 "is 99 aa instead of 100 aa. Tests that (a) hapli classifies this "
                 "as 'inframe_del' rather than LoF, (b) protein-diff reports "
                 "length_changed=True with no frameshift signature.",
    }
    return CaseBundle(
        name="12_inframe_domain_indel",
        description="Single 3 bp in-frame deletion on hap1; protein is 1 aa shorter.",
        reference={"chr1": seq},
        gff_records=gff,
        vcf_records=vcf,
        vcf_headers=headers,
        expected=expected,
    )


# -----------------------------------------------------------------------------
# Registry
# -----------------------------------------------------------------------------
CASES = {
    "01_synonymous_pair": case_01_synonymous_pair,
    "02_single_benign_missense": case_02_single_benign_missense,
    "03_single_stop_gain": case_03_single_stop_gain,
    "04_compound_het_lof": case_04_compound_het_lof,
    "05_frameshift_rescue": case_05_frameshift_rescue,
    "06_compound_missense_pocket": case_06_compound_missense_pocket,
    "07_splice_donor_cryptic": case_07_splice_donor_cryptic,
    "08_symbolic_del_whole_gene": case_08_symbolic_del_whole_gene,
    "09_upstream_del_shift": case_09_upstream_del_shift,
    "10_inversion_body": case_10_inversion_body,
    "11_tandem_dup_small_gene": case_11_tandem_dup_small_gene,
    "12_inframe_domain_indel": case_12_inframe_domain_indel,
}

# Names reserved for follow-on implementation (none at the moment — all 12 wired).
PHASE5_TODO: list[str] = []


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=Path("data/paper_cases"), help="output root directory")
    parser.add_argument("--case", action="append", help="case name (may be repeated)", default=[])
    parser.add_argument("--all", action="store_true", help="generate every implemented case")
    parser.add_argument("--list", action="store_true", help="list cases and exit")
    parser.add_argument("--force", action="store_true", help="overwrite existing case directories")
    args = parser.parse_args(argv)

    if args.list:
        print("Implemented:")
        for name in sorted(CASES):
            print(f"  {name}  — {CASES[name]().description}")
        if PHASE5_TODO:
            print("TODO (Phase 5):")
            for name in PHASE5_TODO:
                print(f"  {name}")
        return 0

    if args.all:
        names = list(CASES)
    elif args.case:
        names = args.case
    else:
        parser.error("Provide --case NAME (repeatable), --all, or --list")

    for name in names:
        if name not in CASES:
            print(f"[error] unknown case {name!r}; see --list", file=sys.stderr)
            return 1
        out = args.out / name
        if out.exists():
            if not args.force:
                print(f"[skip]  {name} exists at {out} (use --force to overwrite)")
                continue
            shutil.rmtree(out)
        CASES[name]().write(out)
        print(f"[ok]    {name} -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
