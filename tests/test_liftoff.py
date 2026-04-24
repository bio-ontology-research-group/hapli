"""Tests for hapli/external/liftoff.py wrapper.

Skips gracefully if the `liftoff` binary isn't on PATH, so CI without the
binary installed (or users who haven't run `uv tool install`) don't see
red. The end-to-end paper-cases tests in tests/test_paper_cases.py will
also catch regressions once Liftoff is integrated into the pipeline.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from hapli.external.consensus import consensus_region
from hapli.external.liftoff import (
    LiftoffNotAvailable,
    parse_lifted_gff,
    run_liftoff,
)

pytestmark = pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff binary not installed; run `uv tool install git+https://github.com/agshumate/Liftoff`",
)


@pytest.fixture(scope="module")
def case_08_bundle(tmp_path_factory) -> Path:
    """Generate case 08 (SV-deletes-whole-gene) once per module."""
    out = tmp_path_factory.mktemp("liftoff_cases")
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).parent.parent / "scripts" / "generate_paper_cases.py"),
            "--case", "08_symbolic_del_whole_gene",
            "--out", str(out),
            "--force",
        ],
        check=True,
        capture_output=True,
    )
    bundle = out / "08_symbolic_del_whole_gene"
    # Materialise haplotype FASTAs
    haps = consensus_region(
        bundle / "reference.fa",
        bundle / "phased.vcf.gz",
        sample="S1",
        region="chr1:1-3000",
    )
    for name, seq in haps.items():
        p = bundle / f"{name}.fa"
        with p.open("w") as f:
            f.write(f">chr1\n{seq}\n")
        pysam.faidx(str(p))
    return bundle


def test_deleted_gene_reported_as_deleted(case_08_bundle: Path, tmp_path: Path):
    """hap1 has the symbolic <DEL> that removes G1 — Liftoff must flag it absent."""
    result = run_liftoff(
        haplotype_fasta=case_08_bundle / "hap1.fa",
        reference_fasta=case_08_bundle / "reference.fa",
        reference_gff=case_08_bundle / "annotation.gff3",
        out_gff=tmp_path / "hap1.gff3",
    )
    assert "G1" in result.absent_gene_ids
    assert "G1" not in result.present_gene_ids
    assert result.presence["G1"].status == "deleted"


def test_intact_gene_reported_as_intact(case_08_bundle: Path, tmp_path: Path):
    """hap2 has no SV at the gene locus — Liftoff should lift G1 cleanly."""
    result = run_liftoff(
        haplotype_fasta=case_08_bundle / "hap2.fa",
        reference_fasta=case_08_bundle / "reference.fa",
        reference_gff=case_08_bundle / "annotation.gff3",
        out_gff=tmp_path / "hap2.gff3",
    )
    assert "G1" in result.present_gene_ids
    call = result.presence["G1"]
    assert call.status == "intact"
    assert call.coverage == 1.0
    assert call.sequence_identity == 1.0
    assert call.seqid == "chr1"


def test_missing_binary_raises_clear_error(monkeypatch, tmp_path: Path, case_08_bundle: Path):
    """If liftoff isn't on PATH the wrapper must surface a helpful error."""
    monkeypatch.setenv("PATH", "")
    with pytest.raises(LiftoffNotAvailable, match="not on PATH"):
        run_liftoff(
            haplotype_fasta=case_08_bundle / "hap1.fa",
            reference_fasta=case_08_bundle / "reference.fa",
            reference_gff=case_08_bundle / "annotation.gff3",
            out_gff=tmp_path / "x.gff3",
        )


def test_parse_lifted_gff_handles_partial_mapping(tmp_path: Path):
    """Directly exercise the attribute parser — no Liftoff run."""
    g = tmp_path / "partial.gff3"
    g.write_text(
        "##gff-version 3\n"
        "chr1\tLiftoff\tgene\t100\t500\t.\t+\t.\t"
        "ID=GENE1;Name=GENE1;coverage=0.6;sequence_ID=0.98;valid_ORFs=0;"
        "extra_copy_number=0;partial_mapping=True\n"
    )
    out = parse_lifted_gff(g)
    assert "GENE1" in out
    call = out["GENE1"]
    assert call.status == "partial"
    assert call.coverage == 0.6
    assert call.flags.get("partial_mapping") == "True"
