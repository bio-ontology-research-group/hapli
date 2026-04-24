"""Tests for hapli/external/csq.py wrapper."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from hapli.external.csq import (
    _decode_bitmask,
    _parse_bcsq_entry,
    group_by_haplotype,
    run_csq,
)

pytestmark = pytest.mark.skipif(
    shutil.which("bcftools") is None,
    reason="bcftools not on PATH",
)


def test_parse_bcsq_entry_splits_seven_fields():
    e = _parse_bcsq_entry("missense|G1|T1|protein_coding|+|4A>4R|100A>G")
    assert e["consequence"] == "missense"
    assert e["gene"] == "G1"
    assert e["transcript"] == "T1"
    assert e["strand"] == "+"
    assert e["amino_acid_change"] == "4A>4R"
    assert e["dna_change"] == "100A>G"


def test_parse_bcsq_entry_detects_backreference():
    e = _parse_bcsq_entry("@1012")
    assert e == {"backref": 1012}


def test_parse_bcsq_entry_pads_missing_tail_fields():
    e = _parse_bcsq_entry("synonymous|G1|T1|protein_coding|+")
    assert e["amino_acid_change"] == ""
    assert e["dna_change"] == ""


def test_decode_bitmask_interleaves_haplotypes():
    # 5 = 0b101 -> bit 0 (BCSQ #1 on hap1) + bit 2 (BCSQ #2 on hap1)
    out = _decode_bitmask(5, n_entries=3)
    assert out == {1: {1}, 2: {1}}

    # 3 = 0b11 -> bit 0 (BCSQ #1 on hap1) + bit 1 (BCSQ #1 on hap2)
    out = _decode_bitmask(3, n_entries=3)
    assert out == {1: {1, 2}}


@pytest.fixture(scope="module")
def case_05_bundle(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("csq_cases")
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).parent.parent / "scripts" / "generate_paper_cases.py"),
            "--case", "05_frameshift_rescue",
            "--out", str(out),
            "--force",
        ],
        check=True,
        capture_output=True,
    )
    return out / "05_frameshift_rescue"


def test_frameshift_rescue_joint_call_lands_only_on_hap1(case_05_bundle: Path, tmp_path: Path):
    """Both variants are 1|0 — both effects belong on hap1, none on hap2."""
    result = run_csq(
        vcf_path=case_05_bundle / "phased.vcf.gz",
        reference_fasta=case_05_bundle / "reference.fa",
        gff3=case_05_bundle / "annotation.gff3",
        sample="S1",
        out_vcf=tmp_path / "csq.vcf.gz",
    )
    assert len(result.consequences) > 0
    grouped = group_by_haplotype(result.consequences)
    assert grouped[1], "hap1 must receive consequences for two 1|0 variants"
    assert not grouped[2], "hap2 must receive no consequences (both variants are 1|0)"


def test_frameshift_rescue_gets_inframe_altering_tag(case_05_bundle: Path, tmp_path: Path):
    """bcftools csq's haplotype-aware logic recognises the +1/-1 pair as in-frame,
    not double-frameshift. We pass this recognition through verbatim and use it
    as the primary consequence label for the compound block.
    """
    result = run_csq(
        vcf_path=case_05_bundle / "phased.vcf.gz",
        reference_fasta=case_05_bundle / "reference.fa",
        gff3=case_05_bundle / "annotation.gff3",
        sample="S1",
        out_vcf=tmp_path / "csq.vcf.gz",
    )
    consequences = {c.consequence for c in result.consequences}
    assert "inframe_altering" in consequences, (
        "bcftools csq should annotate the +1/-1 pair as inframe_altering"
    )


def test_compound_id_groups_rescue_pair(case_05_bundle: Path, tmp_path: Path):
    """The compound_id is the joined dna_change string; both variants share it."""
    result = run_csq(
        vcf_path=case_05_bundle / "phased.vcf.gz",
        reference_fasta=case_05_bundle / "reference.fa",
        gff3=case_05_bundle / "annotation.gff3",
        sample="S1",
        out_vcf=tmp_path / "csq.vcf.gz",
    )
    compound_ids = {c.compound_id for c in result.consequences if c.compound_id}
    # Both variants should share exactly one compound_id
    assert len(compound_ids) == 1
    assert "+" in next(iter(compound_ids))   # joined notation e.g. "1012T>TC+1029CT>C"
