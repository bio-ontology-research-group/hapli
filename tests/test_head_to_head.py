"""Smoke test for the Phase 5-B head-to-head benchmark."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_head_to_head_runs_end_to_end(tmp_path: Path):
    """Run the head-to-head driver against the paper-cases bundles and verify
    each disagreement axis the paper claims is exercised at least once:

      * compound_het_lof_only_in_hapli  (case 04)
      * presence=deleted_only_in_hapli  (case 08)

    The epistasis axis requires --with-esm; we don't exercise it here to keep
    this test fast and CPU-only-friendly. The dedicated test
    `test_frameshift_rescue_emits_positive_epistasis_residual` covers it.
    """
    repo_root = Path(__file__).parent.parent

    # Generate fixtures
    cases_dir = tmp_path / "paper_cases"
    subprocess.run(
        [sys.executable, str(repo_root / "scripts" / "generate_paper_cases.py"),
         "--all", "--out", str(cases_dir), "--force"],
        check=True, capture_output=True,
    )

    out = tmp_path / "head_to_head.json"
    subprocess.run(
        ["uv", "run", "python3", "benchmarks/gnomad_flips/run_head_to_head.py",
         "--cases-dir", str(cases_dir),
         "--out", str(out)],
        cwd=repo_root, check=True, capture_output=True, timeout=600,
    )

    summary = json.loads(out.read_text())
    by_case = {r["case"]: r for r in summary["rows"]}

    # Case 04 must produce hapli's compound_het_lof signal
    chl_case = by_case["04_compound_het_lof"]
    assert chl_case["hapli"]["compound_het_lof"] is True
    assert "compound_het_lof_only_in_hapli" in chl_case["hapli_only_signals"]

    # Case 08 must produce hapli's presence-deleted signal
    sv_case = by_case["08_symbolic_del_whole_gene"]
    assert sv_case["hapli"]["presence_hap1"] == "deleted"
    assert any("deleted" in s for s in sv_case["hapli_only_signals"])

    # Tally must be sane
    assert summary["tally"]["hapli_compound_het_lof_calls"] >= 1
    assert summary["tally"]["csq_lof_haplotype_records"] >= 2  # case 04 has 2 LoF haps
