"""Shared pytest fixtures for hapli tests."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pysam
import pytest


@pytest.fixture
def tmp_fasta(tmp_path: Path):
    """Factory for writing a small indexed FASTA inside tmp_path."""

    def _make(name: str, records: dict[str, str]) -> Path:
        fa = tmp_path / f"{name}.fa"
        with fa.open("w") as f:
            for seqid, seq in records.items():
                f.write(f">{seqid}\n")
                for i in range(0, len(seq), 60):
                    f.write(seq[i : i + 60] + "\n")
        subprocess.run(["samtools", "faidx", str(fa)], check=True, capture_output=True)
        return fa

    return _make


@pytest.fixture
def tmp_vcf(tmp_path: Path):
    """Factory for writing a tiny bgzipped + tabix-indexed VCF from a header + record lines.

    Uses pysam.tabix_compress / pysam.tabix_index so no standalone bgzip/tabix
    binary is required — samtools and pysam ship the implementation already.
    """

    def _make(name: str, header_lines: list[str], record_lines: list[str]) -> Path:
        vcf = tmp_path / f"{name}.vcf"
        with vcf.open("w") as f:
            f.write("##fileformat=VCFv4.2\n")
            for h in header_lines:
                if not h.endswith("\n"):
                    h = h + "\n"
                f.write(h)
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
            for r in record_lines:
                if not r.endswith("\n"):
                    r = r + "\n"
                f.write(r)
        gz = tmp_path / f"{name}.vcf.gz"
        pysam.tabix_compress(str(vcf), str(gz), force=True)
        pysam.tabix_index(str(gz), preset="vcf", force=True)
        return gz

    return _make


def _require_tool(name: str) -> None:
    if shutil.which(name) is None:
        pytest.skip(f"{name} not on PATH")


@pytest.fixture(autouse=True)
def _require_external_tools():
    for t in ("samtools", "bcftools"):
        _require_tool(t)
