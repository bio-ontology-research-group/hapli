"""Tests for the ClinVar lookup wrapper."""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from hapli.external.clinvar import (
    ClinVarHit,
    ClinVarLookup,
    _normalise_clnsig,
)


# ─────────────────────────────────────────────────────────────────────────────
# Pure-python helpers
# ─────────────────────────────────────────────────────────────────────────────
def test_normalise_clnsig_collapses_priority():
    # Multi-significance string gets the most-pathogenic alias
    assert _normalise_clnsig("Pathogenic,Likely_pathogenic") == "pathogenic"
    assert _normalise_clnsig("Likely_benign,Benign") == "likely_benign"
    assert _normalise_clnsig("Uncertain_significance") == "vus"
    assert _normalise_clnsig("Benign") == "benign"
    assert _normalise_clnsig(None) is None
    assert _normalise_clnsig("") is None


# ─────────────────────────────────────────────────────────────────────────────
# VCF lookup with synthetic ClinVar VCF
# ─────────────────────────────────────────────────────────────────────────────
def _write_vcf(path: Path, records: list[tuple]) -> Path:
    """Write a tabix-indexed VCF with ClinVar-shaped INFO."""
    raw = path.with_suffix("")
    with raw.open("w") as f:
        f.write(
            "##fileformat=VCFv4.2\n"
            '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">\n'
            '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review status">\n'
            '##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene info">\n'
            '##INFO=<ID=CLNACC,Number=.,Type=String,Description="Accession">\n'
            "##contig=<ID=1,length=100>\n"
            "##contig=<ID=2,length=100>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        for chrom, pos, rid, ref, alt, info in records:
            f.write(f"{chrom}\t{pos}\t{rid}\t{ref}\t{alt}\t.\tPASS\t{info}\n")
    pysam.tabix_compress(str(raw), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return path


def test_clinvar_vcf_lookup_returns_normalised_hit(tmp_path):
    p = _write_vcf(tmp_path / "clinvar.vcf.gz", [
        ("1", 100, "rs1", "A", "G", "CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,multiple_submitters;GENEINFO=BRCA1:672"),
        ("1", 200, "rs2", "C", "T", "CLNSIG=Likely_benign;GENEINFO=BRCA1:672"),
        ("1", 300, "rs3", "G", "A", "CLNSIG=Pathogenic,Likely_pathogenic;GENEINFO=TP53:7157"),
    ])
    with ClinVarLookup.from_vcf(p) as cv:
        # query with chr-prefix → wrapper should strip
        h = cv.get("chr1", 100, "A", "G")
        assert h is not None
        assert h.clnsig == "pathogenic"
        assert h.raw_clnsig == "Pathogenic"
        assert h.gene_symbol == "BRCA1"

        h2 = cv.get("1", 200, "C", "T")
        assert h2 is not None
        assert h2.clnsig == "likely_benign"

        # Multi-significance → collapses to most pathogenic
        h3 = cv.get("1", 300, "G", "A")
        assert h3.clnsig == "pathogenic"


def test_clinvar_vcf_miss_returns_none(tmp_path):
    p = _write_vcf(tmp_path / "clinvar.vcf.gz", [
        ("1", 100, "rs1", "A", "G", "CLNSIG=Pathogenic"),
    ])
    with ClinVarLookup.from_vcf(p) as cv:
        # Wrong ALT
        assert cv.get("1", 100, "A", "T") is None
        # Wrong POS
        assert cv.get("1", 101, "A", "G") is None
        # Wrong contig
        assert cv.get("99", 100, "A", "G") is None


def test_clinvar_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        ClinVarLookup.from_vcf(tmp_path / "nope.vcf.gz")
