"""Tests for the AlphaMissense lookup wrapper.

We ship a tiny synthetic TSV to exercise the tabix lookup without needing
the 11 GB production table. The real download procedure is documented in
`hapli.external.alphamissense.download_hint`.
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from hapli.external.alphamissense import (
    DEFAULT_COLUMNS,
    AlphaMissenseLookup,
)


def _write_fake_am(tmp_path: Path, chr_style: str = "chr", rows: list[tuple] | None = None) -> Path:
    """Write a tabix-indexed AlphaMissense-shaped TSV.

    `chr_style` = "chr" → CHROM column contains 'chr1'; "bare" → CHROM column contains '1'.
    """
    rows = rows or [
        # (CHROM, POS, REF, ALT, genome, uniprot, transcript, protein_variant, am_pathogenicity, am_class)
        ("1", 100, "A", "G", "hg38", "P12345", "ENST0001", "M34V", 0.12, "likely_benign"),
        ("1", 100, "A", "C", "hg38", "P12345", "ENST0001", "M34L", 0.88, "likely_pathogenic"),
        ("1", 200, "C", "T", "hg38", "P12345", "ENST0001", "R67W", 0.55, "ambiguous"),
        ("2", 500, "G", "A", "hg38", "P99999", "ENST0002", "E12K", 0.04, "likely_benign"),
    ]
    raw = tmp_path / "am.tsv"
    with raw.open("w") as f:
        f.write("#" + "\t".join(DEFAULT_COLUMNS) + "\n")
        for r in rows:
            chrom = r[0]
            if chr_style == "chr" and not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            elif chr_style == "bare" and chrom.startswith("chr"):
                chrom = chrom[3:]
            f.write("\t".join(str(c) if i != 0 else chrom for i, c in enumerate(r)) + "\n")
    gz = tmp_path / "am.tsv.gz"
    pysam.tabix_compress(str(raw), str(gz), force=True)
    pysam.tabix_index(
        str(gz),
        seq_col=0, start_col=1, end_col=1,
        zerobased=False, force=True, meta_char="#",
    )
    return gz


def test_basic_lookup_hits_and_misses(tmp_path):
    gz = _write_fake_am(tmp_path, chr_style="chr")
    with AlphaMissenseLookup(gz) as am:
        hit = am.get("chr1", 100, "A", "G")
        assert hit is not None
        assert hit.score == 0.12
        assert hit.am_class == "likely_benign"

        hit = am.get("chr1", 100, "A", "C")
        assert hit is not None
        assert hit.am_class == "likely_pathogenic"

        assert am.get("chr1", 100, "A", "T") is None   # no row for this ALT
        assert am.get("chr1", 101, "A", "G") is None   # wrong POS
        assert am.get("chr99", 1, "A", "G") is None    # absent contig


def test_chrom_prefix_auto_detection_strip(tmp_path):
    """Table uses bare '1'; query uses 'chr1'. Lookup should normalise."""
    gz = _write_fake_am(tmp_path, chr_style="bare")
    with AlphaMissenseLookup(gz) as am:
        hit = am.get("chr1", 100, "A", "G")
        assert hit is not None
        assert hit.score == 0.12


def test_chrom_prefix_auto_detection_add(tmp_path):
    """Table uses 'chr1'; query uses bare '1'. Lookup should normalise."""
    gz = _write_fake_am(tmp_path, chr_style="chr")
    with AlphaMissenseLookup(gz) as am:
        hit = am.get("1", 100, "A", "G")
        assert hit is not None
        assert hit.score == 0.12


def test_batch_lookup(tmp_path):
    gz = _write_fake_am(tmp_path, chr_style="chr")
    with AlphaMissenseLookup(gz) as am:
        result = am.get_many([
            ("chr1", 100, "A", "G"),
            ("chr1", 200, "C", "T"),
            ("chr1", 999, "A", "G"),   # miss
        ])
    assert result[("chr1", 100, "A", "G")].score == 0.12
    assert result[("chr1", 200, "C", "T")].am_class == "ambiguous"
    assert result[("chr1", 999, "A", "G")] is None


def test_missing_table_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        AlphaMissenseLookup(tmp_path / "nonexistent.tsv.gz")
