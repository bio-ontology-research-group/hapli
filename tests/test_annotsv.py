"""Tests for hapli/external/annotsv.py.

Skips the live-tool integration if `AnnotSV` is not on PATH (CI doesn't
have it). The TSV parser is exercised directly against frozen fixtures
so the column-mapping logic is regression-tested without needing the
binary.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from hapli.external.annotsv import (
    AnnotSVCall,
    AnnotSVNotAvailable,
    parse_annotsv_tsv,
    run_annotsv,
)


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> Path:
    with path.open("w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    return path


def test_parse_annotsv_tsv_full_mode_simple_del(tmp_path: Path):
    """A single SV in 'full' annotation mode — one row, one gene."""
    tsv = _write_tsv(
        tmp_path / "out.annotsv.tsv",
        header=[
            "AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_type",
            "Gene_name", "Location", "Annotation_mode", "ACMG_class",
            "Tx_start", "Tx_end",
        ],
        rows=[[
            "AnnotSV_1", "chr3", "10141635", "10153670", "DEL",
            "VHL", "txStart-txEnd", "full", "5",
            "10141635", "10153670",
        ]],
    )
    calls = list(parse_annotsv_tsv(tsv))
    assert len(calls) == 1
    c = calls[0]
    assert c.sv_chrom == "chr3"
    assert c.sv_start == 10141635
    assert c.sv_end == 10153670
    assert c.sv_type == "DEL"
    assert c.gene_name == "VHL"
    assert c.acmg_class == 5
    assert c.annotation_mode == "full"
    # Extra columns preserved verbatim
    assert c.extra["Tx_start"] == "10141635"
    assert c.extra["Tx_end"] == "10153670"


def test_parse_annotsv_tsv_split_mode_multi_gene(tmp_path: Path):
    """A larger SV that overlaps two genes — one 'full' row + two 'split' rows.
    AnnotSV emits the gene-level breakdown as separate rows; we yield one
    AnnotSVCall per row.
    """
    tsv = _write_tsv(
        tmp_path / "out.annotsv.tsv",
        header=["AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_type",
                "Gene_name", "Location", "Annotation_mode", "ACMG_class"],
        rows=[
            ["A_1", "chr1", "100000", "500000", "DEL", "GENEA;GENEB", "txStart-txEnd", "full", "4"],
            ["A_1", "chr1", "100000", "500000", "DEL", "GENEA",        "txStart-txEnd", "split", "4"],
            ["A_1", "chr1", "100000", "500000", "DEL", "GENEB",        "exon3-txEnd",   "split", "3"],
        ],
    )
    calls = list(parse_annotsv_tsv(tsv))
    assert len(calls) == 3
    assert calls[0].annotation_mode == "full"
    assert calls[1].annotation_mode == "split" and calls[1].gene_name == "GENEA"
    assert calls[2].annotation_mode == "split" and calls[2].gene_name == "GENEB"


def test_parse_annotsv_tsv_handles_NA_acmg(tmp_path: Path):
    """ACMG column may be 'NA' or empty — parser must coerce to None."""
    tsv = _write_tsv(
        tmp_path / "out.annotsv.tsv",
        header=["AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_type",
                "Gene_name", "Location", "Annotation_mode", "ACMG_class"],
        rows=[
            ["A_1", "chr2", "1", "2", "INS", "G1", "intronic", "split", "NA"],
            ["A_2", "chr2", "5", "6", "INV", "G2", "intronic", "split", ""],
        ],
    )
    calls = list(parse_annotsv_tsv(tsv))
    assert calls[0].acmg_class is None
    assert calls[1].acmg_class is None


def test_parse_annotsv_tsv_alias_for_acmg_column(tmp_path: Path):
    """v3.3 used `AnnotSV_ranking_score`, v3.4 uses `ACMG_class`. Tolerate either."""
    tsv = _write_tsv(
        tmp_path / "out.annotsv.tsv",
        header=["AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_type",
                "Gene_name", "Location", "Annotation_mode", "AnnotSV_ranking_score"],
        rows=[["A_1", "chr3", "1", "2", "DEL", "G", "exonic", "split", "3"]],
    )
    calls = list(parse_annotsv_tsv(tsv))
    assert calls[0].acmg_class == 3


def test_parse_annotsv_tsv_missing_essential_column_raises(tmp_path: Path):
    tsv = _write_tsv(
        tmp_path / "out.annotsv.tsv",
        header=["foo", "bar"],
        rows=[["1", "2"]],
    )
    with pytest.raises(ValueError, match="missing essential columns"):
        list(parse_annotsv_tsv(tsv))


def test_run_annotsv_missing_binary_raises(monkeypatch, tmp_path: Path):
    """If AnnotSV isn't on PATH, surface a helpful error."""
    monkeypatch.setenv("PATH", "")
    with pytest.raises(AnnotSVNotAvailable, match="not on PATH"):
        run_annotsv(
            vcf_path=tmp_path / "x.vcf.gz",
            out_tsv=tmp_path / "x.tsv",
        )


@pytest.mark.skipif(
    shutil.which("AnnotSV") is None,
    reason="AnnotSV binary not on PATH; conda install -c bioconda annotsv",
)
def test_run_annotsv_smoke_on_paper_case_08(tmp_path: Path):
    """End-to-end live-tool smoke if AnnotSV is installed.
    Generates paper case 08 (symbolic <DEL> over G1) and expects AnnotSV
    to recognise the DEL and assign G1 a non-empty location.
    """
    import subprocess
    import sys
    out = tmp_path / "case08"
    subprocess.run(
        [sys.executable,
         str(Path(__file__).parent.parent / "scripts" / "generate_paper_cases.py"),
         "--case", "08_symbolic_del_whole_gene",
         "--out", str(out),
         "--force"],
        check=True, capture_output=True,
    )
    bundle = out / "08_symbolic_del_whole_gene"
    result = run_annotsv(
        vcf_path=bundle / "phased.vcf.gz",
        out_tsv=tmp_path / "out.annotsv.tsv",
    )
    g1_rows = [c for c in result.per_sv_calls if c.gene_name == "G1"]
    assert g1_rows, "AnnotSV should produce at least one G1 row for the <DEL>"
    assert any(c.sv_type == "DEL" for c in g1_rows)
