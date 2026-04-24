"""End-to-end test that AlphaMissense lookup is wired into the pipeline."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from hapli.external.alphamissense import DEFAULT_COLUMNS


def _build_synthetic_am_table(tmp_path: Path) -> Path:
    """Tabix-indexed AM table covering case 06's MNV positions on chr1."""
    raw = tmp_path / "am.tsv"
    with raw.open("w") as f:
        f.write("#" + "\t".join(DEFAULT_COLUMNS) + "\n")
        # Cover positions 1028-1030 and 1088-1090 with all SNVs (max_am = 0.9).
        for pos in (1028, 1029, 1030, 1088, 1089, 1090):
            for ref in "ACGT":
                for alt in "ACGT":
                    if ref == alt:
                        continue
                    f.write("\t".join([
                        "chr1", str(pos), ref, alt, "hg38",
                        "P_SYNTH", "T_SYNTH", f"X1{alt}",
                        "0.9", "likely_pathogenic",
                    ]) + "\n")
    gz = tmp_path / "am.tsv.gz"
    pysam.tabix_compress(str(raw), str(gz), force=True)
    pysam.tabix_index(str(gz),
                      seq_col=0, start_col=1, end_col=1,
                      zerobased=False, force=True, meta_char="#")
    return gz


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_alphamissense_aggregation_populates_missense_agg(tmp_path):
    """Run case 06 (compound missense) with --alphamissense-table and check
    that evidence.missense_agg has one record per (transcript, haplotype) with
    the synthetic per-MNV-max scores aggregated."""
    repo_root = Path(__file__).parent.parent
    am_table = _build_synthetic_am_table(tmp_path)

    # Generate case 06 fixture
    cases_dir = tmp_path / "paper_cases"
    subprocess.run(
        [sys.executable, str(repo_root / "scripts" / "generate_paper_cases.py"),
         "--case", "06_compound_missense_pocket",
         "--out", str(cases_dir), "--force"],
        check=True, capture_output=True,
    )
    case_dir = cases_dir / "06_compound_missense_pocket"
    out_dir = tmp_path / "results"

    subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1", "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
            "--alphamissense-table", str(am_table),
        ],
        cwd=repo_root, check=True, capture_output=True, timeout=120,
    )
    data = json.loads((out_dir / "S1_G1_analysis.json").read_text())
    agg = data["evidence"]["missense_agg"]
    # Both case-06 missense are on hap1 → expect one bucket
    assert len(agg) == 1, f"expected 1 bucket, got {agg}"
    bucket = agg[0]
    assert bucket["transcript"] == "T1"
    assert bucket["haplotype"] == 1
    assert bucket["n_missense"] == 2          # two MNVs aggregated
    assert bucket["max_am"] == 0.9
    assert bucket["mean_am"] == 0.9
    assert bucket["source"] == "alphamissense_v1"


def test_alphamissense_no_table_no_op(tmp_path):
    """Without --alphamissense-table, the pipeline must not populate missense_agg."""
    # We exercise the no-op via the unit test on _aggregate_alphamissense
    # without spinning up the full pipeline. Instantiate HapliPipeline with
    # alphamissense_table=None and confirm the call is a no-op.
    from hapli.workflow.pipeline import HapliPipeline
    from hapli.core.schema import GeneEvidence, ConsequenceCall

    pipeline = HapliPipeline(
        ref_path=tmp_path / "fake.fa",  # don't need to exist for this method
        vcf_path=None,
        gff_path=tmp_path / "fake.gff3",
        output_dir=tmp_path,
        alphamissense_table=None,
    )
    ev = GeneEvidence(gene="G1")
    ev.consequence.append(ConsequenceCall(
        chrom="chr1", pos=100, ref="A", alt="G",
        haplotype=1, transcript="T1", consequence="missense",
    ))
    pipeline._aggregate_alphamissense(ev)
    assert ev.missense_agg == []
