"""
Regression wrappers for the two population-scale benchmark smoke tests:

  benchmarks/1000g_acmg/smoke_test.py   — Mode A, 5 samples × 3 genes
  benchmarks/hprc/smoke_test.py         — Mode B, 3 samples × 2 genes

Each runs in <1 min and exercises the full Snakemake → analyze.py path
on a synthetic but biologically-shaped fixture; if either breaks, the
cluster-scale benchmark won't reach the cluster.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).parent.parent


@pytest.mark.skipif(
    shutil.which("liftoff") is None or shutil.which("snakemake") is None,
    reason="liftoff and/or snakemake not on PATH",
)
def test_1000g_acmg_smoke():
    r = subprocess.run(
        [sys.executable, str(REPO / "benchmarks" / "1000g_acmg" / "smoke_test.py")],
        capture_output=True, text=True, timeout=300,
    )
    assert r.returncode == 0, (
        f"1000g_acmg smoke failed (rc={r.returncode})\n"
        f"stdout:\n{r.stdout[-2000:]}\n"
        f"stderr:\n{r.stderr[-2000:]}"
    )
    assert "all assertions passed" in r.stderr


@pytest.mark.skipif(
    shutil.which("liftoff") is None or shutil.which("snakemake") is None,
    reason="liftoff and/or snakemake not on PATH",
)
def test_hprc_smoke():
    r = subprocess.run(
        [sys.executable, str(REPO / "benchmarks" / "hprc" / "smoke_test.py")],
        capture_output=True, text=True, timeout=300,
    )
    assert r.returncode == 0, (
        f"hprc smoke failed (rc={r.returncode})\n"
        f"stdout:\n{r.stdout[-2000:]}\n"
        f"stderr:\n{r.stderr[-2000:]}"
    )
    assert "all assertions passed" in r.stderr
