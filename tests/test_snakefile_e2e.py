"""
End-to-end smoke test for the Snakemake workflow.

Pins the two-stage pipeline: per-(sample, gene) hapli analyze, then the
aggregate rule producing per-sample and per-gene TSVs. This is the workflow
that scales the tool to population-level runs; without this test, a refactor
in the Snakefile or aggregator could silently break the workflow without
triggering any of the unit-level tests.

Uses the in-repo `data/paper_cases/04_compound_het_lof` bundle as the input —
it produces a clean compound_het_lof=True signal, so the per-gene TSV has a
predictable headline number to assert against.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest


@pytest.mark.skipif(
    shutil.which("liftoff") is None or shutil.which("snakemake") is None,
    reason="liftoff and/or snakemake not available",
)
def test_snakefile_end_to_end_on_paper_case_04(tmp_path: Path):
    repo_root = Path(__file__).parent.parent
    case_dir = repo_root / "data" / "paper_cases" / "04_compound_het_lof"
    if not case_dir.exists():
        # Generate on demand if the user hasn't run the case generator locally
        subprocess.run(
            [sys.executable, str(repo_root / "scripts" / "generate_paper_cases.py"),
             "--case", "04_compound_het_lof",
             "--out", str(repo_root / "data" / "paper_cases"),
             "--force"],
            check=True,
        )

    # Build a dedicated configfile in tmp_path that points at the case-04 bundle
    # and writes results into a run-specific subdirectory, so this test does
    # not stomp on anyone's results/workflow_output/.
    run_out = tmp_path / "run_out"
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "mode: A\n"
        f"output_dir: {run_out}\n"
        f"reference: {case_dir / 'reference.fa'}\n"
        f"gff: {case_dir / 'annotation.gff3'}\n"
        f"vcf: {case_dir / 'phased.vcf.gz'}\n"
        "samples:\n"
        "  - S1\n"
        "genes:\n"
        "  - G1\n"
        "with_esm: false\n"
        "threads_per_gene: 1\n"
        "mem_mb_per_gene: 2000\n"
    )

    # Stay in repo_root (so `uv run main.py` resolves) but write outputs under
    # tmp_path via the absolute-path output_dir configured above.
    result = subprocess.run(
        ["uv", "run", "snakemake",
         "-s", str(repo_root / "workflows" / "Snakefile"),
         "--configfile", str(cfg),
         "--cores", "2",
         "--forceall"],
        cwd=repo_root,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"snakemake failed:\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )

    # Per-sample TSV must exist with the one row we expect
    per_sample = run_out / "aggregate" / "per_sample.tsv"
    per_gene = run_out / "aggregate" / "per_gene.tsv"
    assert per_sample.exists(), f"missing {per_sample}"
    assert per_gene.exists(), f"missing {per_gene}"

    with per_sample.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 1, f"expected one sample row, got {len(rows)}"
    r = rows[0]
    assert r["sample"] == "S1"
    assert r["gene"] == "G1"
    # Case 04: compound-het LoF → both haps below threshold, flag=True
    assert r["compound_het_lof"] in ("1", "True", "true"), (
        f"case 04 must flag compound_het_lof; got row={r}"
    )

    # Per-gene TSV: one gene, one sample, compound_het_lof_freq == 1.0
    with per_gene.open() as f:
        grows = list(csv.DictReader(f, delimiter="\t"))
    assert len(grows) == 1
    g = grows[0]
    assert g["gene"] == "G1"
    assert float(g["compound_het_lof_freq"]) == 1.0


@pytest.mark.skipif(
    shutil.which("liftoff") is None or shutil.which("snakemake") is None,
    reason="liftoff and/or snakemake not available",
)
def test_snakefile_mode_b_end_to_end(tmp_path: Path):
    """Mode B (pre-assembled haplotype FASTAs) smoke test: materialise two
    haplotype FASTAs from case 04's VCF via bcftools consensus, then run the
    Snakefile with `mode: B` → aggregate. Verifies the HPRC-style code path
    end-to-end, not just the `hapli assess` subcommand in isolation.
    """
    import pysam
    from hapli.external.consensus import consensus_region

    repo_root = Path(__file__).parent.parent
    case_dir = repo_root / "data" / "paper_cases" / "04_compound_het_lof"
    if not case_dir.exists():
        subprocess.run(
            [sys.executable, str(repo_root / "scripts" / "generate_paper_cases.py"),
             "--case", "04_compound_het_lof",
             "--out", str(repo_root / "data" / "paper_cases"),
             "--force"],
            check=True,
        )

    # Materialise hap1/hap2 FASTAs (simulating HPRC assemblies)
    assembly_dir = tmp_path / "assembly"
    assembly_dir.mkdir()
    haps = consensus_region(
        reference_fasta=case_dir / "reference.fa",
        vcf_path=case_dir / "phased.vcf.gz",
        sample="S1",
        region="chr1:1-3000",
    )
    hap_paths = {}
    for name, seq in haps.items():
        p = assembly_dir / f"HG002.{name}.fa"
        with p.open("w") as f:
            f.write(f">chr1\n{seq}\n")
        pysam.faidx(str(p))
        hap_paths[name] = p

    run_out = tmp_path / "run_out"
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "mode: B\n"
        f"output_dir: {run_out}\n"
        f"reference: {case_dir / 'reference.fa'}\n"
        f"gff: {case_dir / 'annotation.gff3'}\n"
        "haps:\n"
        "  HG002:\n"
        f"    hap1: {hap_paths['hap1']}\n"
        f"    hap2: {hap_paths['hap2']}\n"
        "genes:\n"
        "  - G1\n"
        "with_esm: false\n"
        "threads_per_gene: 1\n"
        "mem_mb_per_gene: 2000\n"
    )

    result = subprocess.run(
        ["uv", "run", "snakemake",
         "-s", str(repo_root / "workflows" / "Snakefile"),
         "--configfile", str(cfg),
         "--cores", "2",
         "--forceall"],
        cwd=repo_root, capture_output=True, text=True, timeout=300,
    )
    assert result.returncode == 0, (
        f"snakemake Mode B failed:\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    # Aggregator outputs exist
    assert (run_out / "aggregate" / "per_sample.tsv").exists()
    assert (run_out / "aggregate" / "per_gene.tsv").exists()
    # Per-sample TSV has the HG002 row
    with (run_out / "aggregate" / "per_sample.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert any(r["sample"] == "HG002" and r["gene"] == "G1" for r in rows)
