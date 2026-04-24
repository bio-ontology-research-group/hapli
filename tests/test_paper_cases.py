"""
End-to-end smoke tests for the motivating cases in `data/paper_cases/`.

Each case bundle is generated on-demand by `scripts/generate_paper_cases.py`.
The tests then run `hapli analyze` against the bundle and assert (a) the
command exits 0, (b) the JSON is produced, (c) the schema-v1 contract holds.

Case 08 (symbolic <DEL> removes whole gene) is the Phase 0 SV-fix acceptance
test: pre-fix, this case crashed or produced garbage; post-fix, it runs clean
and the analysis JSON is emitted.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


CASES = [
    "01_synonymous_pair",
    "02_single_benign_missense",
    "03_single_stop_gain",
    "04_compound_het_lof",
    "05_frameshift_rescue",
    "06_compound_missense_pocket",
    "07_splice_donor_cryptic",
    "08_symbolic_del_whole_gene",
    "09_upstream_del_shift",
    "10_inversion_body",
    "11_tandem_dup_small_gene",
    "12_inframe_domain_indel",
]


@pytest.fixture(scope="module")
def paper_cases_dir(tmp_path_factory) -> Path:
    """Generate all paper cases once per test module into a tmp directory."""
    out = tmp_path_factory.mktemp("paper_cases")
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).parent.parent / "scripts" / "generate_paper_cases.py"),
            "--all",
            "--out",
            str(out),
            "--force",
        ],
        check=True,
        capture_output=True,
    )
    return out


def test_analyze_missing_gene_exits_cleanly(paper_cases_dir: Path, tmp_path: Path):
    """User-facing regression: `hapli analyze --gene NONEXISTENT` must exit 1
    with a clean "Error:" stderr message, not a Python traceback.
    """
    case_dir = paper_cases_dir / "01_synonymous_pair"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    result = subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "NONEXISTENT",   # no such gene in the fixture
            "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root, capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 1, f"expected exit 1; got {result.returncode}"
    # Clean single-line error (not a multi-line traceback)
    assert "NONEXISTENT" in (result.stdout + result.stderr)
    assert "Traceback" not in result.stderr, (
        f"traceback leaked to stderr:\n{result.stderr}"
    )


@pytest.mark.parametrize("case_name", CASES)
def test_analyze_runs_on_case(paper_cases_dir: Path, tmp_path: Path, case_name: str):
    case_dir = paper_cases_dir / case_name
    assert case_dir.exists(), f"case bundle missing: {case_dir}"

    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    result = subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1",
            "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, (
        f"hapli analyze failed on {case_name}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )

    json_path = out_dir / "S1_G1_analysis.json"
    assert json_path.exists(), f"expected {json_path} to be produced"
    data = json.loads(json_path.read_text())

    # Schema v1 contract
    assert data["gene"] == "G1"
    assert data["sample"] == "S1"
    assert "region" in data
    assert isinstance(data["transcripts"], list)
    for t in data["transcripts"]:
        assert "id" in t and "alignments" in t
        assert set(t["alignments"]) == {"hap1", "hap2"}
        for hap in ("hap1", "hap2"):
            rec = t["alignments"][hap]
            # Alignments can be "Unmapped" (for the SV-deletes-gene case on hap1)
            # or a dict with the four required fields.
            if isinstance(rec, dict):
                assert set(rec) >= {"identity", "nm", "cigar", "is_perfect"}


def test_symbolic_del_case_does_not_leak_literal_string(paper_cases_dir: Path, tmp_path: Path):
    """Specific acceptance test for the Phase 0 SV-fix: `<DEL>` must not appear
    as a literal string in the produced haplotype FASTA.
    """
    case_dir = paper_cases_dir / "08_symbolic_del_whole_gene"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1",
            "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root,
        check=True,
        timeout=60,
    )
    hap_fa = (out_dir / "S1_G1_haplotypes.fa").read_text()
    assert "<DEL>" not in hap_fa, "haplotype FASTA must not contain the literal string '<DEL>'"


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_frameshift_rescue_case_protein_is_functional(paper_cases_dir: Path, tmp_path: Path):
    """Phase 2 hero-case acceptance: the +1/-1 haplotype must produce a protein
    that is *functionally preserved* (high identity, no premature stop, the
    substitution cluster is flagged as a restored frameshift window), even
    though a variant-by-variant LoF caller would flag both variants as
    frame-disrupting.
    """
    case_dir = paper_cases_dir / "05_frameshift_rescue"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1",
            "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root,
        check=True,
        timeout=120,
    )
    data = json.loads((out_dir / "S1_G1_analysis.json").read_text())
    proteins = {(p["transcript"], p["haplotype"]): p for p in data["evidence"]["protein"]}
    hap1 = proteins.get(("T1", 1))
    hap2 = proteins.get(("T1", 2))

    assert hap1 is not None, "expected hap1 ProteinDiff for the rescued haplotype"
    assert hap2 is not None, "expected hap2 ProteinDiff for the reference haplotype"

    # hap2 has no variants — protein must match reference exactly.
    assert hap2["identity"] == 1.0
    assert hap2["premature_stop_at"] is None
    assert hap2["frameshift_region"] is None

    # hap1 has the +1/-1 pair. Key paper claims:
    #   1. No premature stop — the protein is full-length.
    #   2. A frameshift_region is detected and flagged as restored.
    #   3. Only a small window of substitutions — identity stays high (>0.7).
    assert hap1["premature_stop_at"] is None, (
        "hap1 must not have a premature stop — the rescue restores the frame"
    )
    assert hap1["frameshift_region"] is not None, (
        "hap1 must flag a frameshift region"
    )
    assert hap1["frameshift_region"]["restored"] is True, (
        "hap1 frameshift region must be recognised as restored"
    )
    assert hap1["identity"] >= 0.7, (
        f"hap1 identity unexpectedly low ({hap1['identity']}) — rescue window "
        f"should leave most of the protein intact"
    )
    # Per-haplotype protein FASTA artifact must be written.
    assert (out_dir / "S1_G1_hap1.pep.fa").exists()
    assert (out_dir / "S1_G1_hap2.pep.fa").exists()
    assert (out_dir / "S1_G1_ref.pep.fa").exists()


import importlib.util
_HAS_ESM = (
    importlib.util.find_spec("torch") is not None
    and importlib.util.find_spec("esm") is not None
)


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_mode_b_assess_produces_schema_v2_on_pre_assembled_haplotypes(paper_cases_dir: Path, tmp_path: Path):
    """Mode B acceptance: the `assess` subcommand must run end-to-end on pre-
    assembled haplotype FASTAs (no VCF), producing the same schema-v2 evidence
    layout as Mode A.

    We materialise "assembly" FASTAs from the case-05 VCF via bcftools consensus
    (simulating the output of hifiasm/Verkko) and feed them to `hapli assess`.
    """
    import pysam
    from hapli.external.consensus import consensus_region

    case_dir = paper_cases_dir / "05_frameshift_rescue"
    assembly_dir = tmp_path / "assembly"
    assembly_dir.mkdir()
    haps = consensus_region(
        reference_fasta=case_dir / "reference.fa",
        vcf_path=case_dir / "phased.vcf.gz",
        sample="S1",
        region="chr1:1-3000",
    )
    for name, seq in haps.items():
        p = assembly_dir / f"{name}.assembly.fa"
        with p.open("w") as f:
            f.write(f">chr1\n{seq}\n")
        pysam.faidx(str(p))

    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    result = subprocess.run(
        [
            "uv", "run", "main.py", "assess",
            "--gene", "G1",
            "--sample", "HG002",
            "--hap1", str(assembly_dir / "hap1.assembly.fa"),
            "--hap2", str(assembly_dir / "hap2.assembly.fa"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root, capture_output=True, text=True, timeout=180,
    )
    assert result.returncode == 0, (
        f"mode B assess failed\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    data = json.loads((out_dir / "HG002_G1_analysis.json").read_text())
    # Schema-v2 evidence must be populated even without a VCF
    presence = data["evidence"]["presence"]
    assert presence["hap1"]["status"] == "intact"
    assert presence["hap2"]["status"] == "intact"
    # Mode B deliberately skips bcftools csq; consequence list is empty
    assert data["evidence"]["consequence"] == []
    # But the protein diff must still detect the rescue window on hap1
    proteins = {(p["transcript"], p["haplotype"]): p for p in data["evidence"]["protein"]}
    hap1 = proteins.get(("T1", 1))
    assert hap1 is not None
    assert hap1["frameshift_region"] is not None
    assert hap1["frameshift_region"]["restored"] is True
    # Diploid reporter works uniformly on Mode B output
    dip = data["evidence"]["diploid"]
    assert dip["hap2_score"] == 1.0
    assert dip["compound_het_lof"] is False


@pytest.mark.skipif(
    shutil.which("liftoff") is None or not _HAS_ESM,
    reason="requires liftoff + fair-esm (the ml extra)",
)
def test_frameshift_rescue_emits_positive_epistasis_residual(paper_cases_dir: Path, tmp_path: Path):
    """Phase 3 research-contribution acceptance: on the +1/-1 rescue haplotype,
    running `analyze --with-esm` must produce an epistasis record whose residual
    is large and positive (joint far less deleterious than additive predicts),
    with flagged=True.
    """
    case_dir = paper_cases_dir / "05_frameshift_rescue"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1",
            "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
            "--with-esm",
        ],
        cwd=repo_root,
        check=True,
        timeout=300,
    )
    data = json.loads((out_dir / "S1_G1_analysis.json").read_text())
    # Epistasis block must be populated.
    epi = {(e["transcript"], e["haplotype"]): e for e in data["evidence"]["epistasis"]}
    hap1 = epi.get(("T1", 1))
    assert hap1 is not None, "expected epistasis record for T1 hap1"
    assert hap1["n_variants"] == 6
    assert hap1["s_joint"] > hap1["s_additive"], (
        "joint must be less deleterious than additive for a rescued window"
    )
    assert hap1["residual"] > 20.0, (
        f"rescue residual unexpectedly small: {hap1['residual']}"
    )
    assert hap1["flagged"] is True
    # hap2 has no substitutions — either absent or n_variants=0.
    hap2 = epi.get(("T1", 2))
    if hap2 is not None:
        assert hap2["n_variants"] == 0
        assert hap2["flagged"] is False


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_upstream_del_shift_protein_is_unchanged(paper_cases_dir: Path, tmp_path: Path):
    """Phase 1 architectural-validation acceptance: a 200 bp deletion 300 nt
    upstream of the gene shifts every reference coordinate on hap1, but hapli's
    Liftoff step finds the gene at its new position and extracts the protein
    identical to the reference. This is the test that confirms the
    "reference-coordinates-invalid-post-consensus" failure mode is handled.
    """
    case_dir = paper_cases_dir / "09_upstream_del_shift"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1", "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root, check=True, timeout=120,
    )
    data = json.loads((out_dir / "S1_G1_analysis.json").read_text())
    presence = data["evidence"]["presence"]
    assert presence["hap1"]["status"] == "intact", (
        f"hap1 should still be intact despite upstream DEL; got {presence['hap1']}"
    )
    assert presence["hap2"]["status"] == "intact"

    proteins = {(p["transcript"], p["haplotype"]): p for p in data["evidence"]["protein"]}
    hap1_protein = proteins.get(("T1", 1))
    assert hap1_protein is not None
    assert hap1_protein["identity"] == 1.0, (
        f"hap1 protein should match reference exactly; got {hap1_protein['identity']}"
    )
    assert hap1_protein["premature_stop_at"] is None
    # Both haps should score 1.0 since neither carries a coding change
    dip = data["evidence"]["diploid"]
    assert dip["hap1_score"] == 1.0 and dip["hap2_score"] == 1.0


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_compound_het_lof_case_emits_diploid_compound_het_lof_true(paper_cases_dir: Path, tmp_path: Path):
    """Phase 5-B acceptance: case 04 (stop-gain hap1 + frameshift hap2) must
    produce diploid.compound_het_lof=True. This is the haplotype-level signal
    that bcftools/csq cannot emit by design.
    """
    case_dir = paper_cases_dir / "04_compound_het_lof"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1",
            "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root, check=True, timeout=120,
    )
    data = json.loads((out_dir / "S1_G1_analysis.json").read_text())
    dip = data["evidence"]["diploid"]
    assert dip["compound_het_lof"] is True, (
        f"case 04 must produce compound_het_lof=True; got dip={dip}"
    )
    # Both haplotypes must score below the LoF threshold
    assert dip["hap1_score"] is not None and dip["hap1_score"] <= 0.5
    assert dip["hap2_score"] is not None and dip["hap2_score"] <= 0.5
    # bcftools csq must independently flag both haps as LoF-class
    consequences = data["evidence"]["consequence"]
    by_hap_lof = {1: False, 2: False}
    for c in consequences:
        if c["consequence"] in ("stop_gained", "frameshift", "splice_acceptor",
                                "splice_donor", "stop_lost", "start_lost"):
            by_hap_lof[c["haplotype"]] = True
    assert by_hap_lof[1], "csq should call hap1 LoF"
    assert by_hap_lof[2], "csq should call hap2 LoF"


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_symbolic_del_case_presence_says_deleted(paper_cases_dir: Path, tmp_path: Path):
    """Phase 1 acceptance: Liftoff must identify G1 as absent on hap1 (where the
    symbolic <DEL> removed it) and intact on hap2. This tests the end-to-end
    Liftoff integration, not just the wrapper.
    """
    case_dir = paper_cases_dir / "08_symbolic_del_whole_gene"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1",
            "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root,
        check=True,
        timeout=120,
    )
    data = json.loads((out_dir / "S1_G1_analysis.json").read_text())
    presence = data["evidence"]["presence"]
    assert presence["hap1"]["status"] == "deleted", (
        f"hap1 should be flagged deleted by Liftoff; got {presence['hap1']}"
    )
    assert presence["hap2"]["status"] == "intact", (
        f"hap2 should be flagged intact by Liftoff; got {presence['hap2']}"
    )


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_single_stop_gain_flags_hap1_lof_only(paper_cases_dir: Path, tmp_path: Path):
    """Case 03 (dominant heterozygous LoF): hap1 carries the stop-gain and must
    score low with a premature-stop detected; hap2 must be intact. The diploid
    aggregator must NOT set compound_het_lof — hap2 is fine.
    """
    case_dir = paper_cases_dir / "03_single_stop_gain"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1", "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root, check=True, timeout=120,
    )
    data = json.loads((out_dir / "S1_G1_analysis.json").read_text())
    proteins = {(p["transcript"], p["haplotype"]): p for p in data["evidence"]["protein"]}
    hap1 = proteins.get(("T1", 1))
    hap2 = proteins.get(("T1", 2))
    assert hap1 is not None and hap2 is not None

    # hap1: premature stop detected, identity drops.
    assert hap1["premature_stop_at"] is not None, "hap1 should have PTC detected"
    assert hap1["identity"] < 1.0
    # hap2: intact.
    assert hap2["premature_stop_at"] is None
    assert hap2["identity"] == 1.0

    dip = data["evidence"]["diploid"]
    assert dip["compound_het_lof"] is False, (
        "single-hap LoF must not trigger compound_het_lof (hap2 is intact)"
    )
    assert dip["hap2_score"] == 1.0


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_single_benign_missense_keeps_both_haps_intact(paper_cases_dir: Path, tmp_path: Path):
    """Case 02 baseline: a single missense must NOT produce an LoF call. Both
    haplotype scores should be high; no PTC, no frameshift, no compound-het flag.
    """
    case_dir = paper_cases_dir / "02_single_benign_missense"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1", "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root, check=True, timeout=120,
    )
    data = json.loads((out_dir / "S1_G1_analysis.json").read_text())
    proteins = {(p["transcript"], p["haplotype"]): p for p in data["evidence"]["protein"]}
    hap1 = proteins.get(("T1", 1))
    hap2 = proteins.get(("T1", 2))
    assert hap1 is not None and hap2 is not None

    # Neither hap has a PTC or frameshift
    assert hap1["premature_stop_at"] is None
    assert hap2["premature_stop_at"] is None
    assert hap1["frameshift_region"] is None
    assert hap2["frameshift_region"] is None
    # hap1 identity slightly below 1.0 (single missense), hap2 == 1.0
    assert hap2["identity"] == 1.0
    assert 0.9 < hap1["identity"] < 1.0, (
        f"single missense should give identity in (0.9, 1.0); got {hap1['identity']}"
    )

    dip = data["evidence"]["diploid"]
    assert dip["compound_het_lof"] is False
    assert dip["hap2_score"] == 1.0
    # hap1 score is identity-based when no PTC → should be near 1.0
    assert dip["hap1_score"] is not None and dip["hap1_score"] > 0.9


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_inframe_del_is_not_lof(paper_cases_dir: Path, tmp_path: Path):
    """Case 12: a 3 bp in-frame deletion must NOT be classified as LoF. Protein
    length on hap1 is 1 aa shorter than reference, no PTC, no frameshift.
    """
    case_dir = paper_cases_dir / "12_inframe_domain_indel"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1", "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root, check=True, timeout=120,
    )
    data = json.loads((out_dir / "S1_G1_analysis.json").read_text())
    proteins = {(p["transcript"], p["haplotype"]): p for p in data["evidence"]["protein"]}
    hap1 = proteins.get(("T1", 1))
    hap2 = proteins.get(("T1", 2))
    assert hap1 is not None and hap2 is not None

    # Key claim: no LoF signals on the in-frame del.
    assert hap1["premature_stop_at"] is None, "in-frame del must not trigger PTC"
    assert hap1["frameshift_region"] is None, "in-frame del must not trigger frameshift flag"
    # Length is changed by exactly 1 aa.
    ref_len = hap2["ref_length"] if "ref_length" in hap2 else None
    hap_len = hap1.get("hap_length")
    if ref_len is not None and hap_len is not None:
        assert ref_len - hap_len == 1, f"expected 1 aa shorter; got ref={ref_len} hap={hap_len}"
    # hap2 untouched.
    assert hap2["identity"] == 1.0


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_tandem_dup_does_not_crash(paper_cases_dir: Path, tmp_path: Path):
    """Case 11: tandem duplication of the gene on hap1. We assert only that
    the pipeline exits cleanly and produces valid JSON — the exact number of
    reported copies depends on bcftools and Liftoff versions (upstream tooling
    variability), and the pipeline's own schema currently collapses copies to
    one lift per haplotype. The acceptance here is "no crash, consistent schema".
    """
    case_dir = paper_cases_dir / "11_tandem_dup_small_gene"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    result = subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1", "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root, capture_output=True, text=True, timeout=180,
    )
    assert result.returncode == 0, (
        f"tandem-dup case crashed:\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    data = json.loads((out_dir / "S1_G1_analysis.json").read_text())
    # Schema sanity
    assert "evidence" in data and "presence" in data["evidence"]


@pytest.mark.skipif(
    shutil.which("liftoff") is None,
    reason="liftoff not installed",
)
def test_inversion_does_not_crash(paper_cases_dir: Path, tmp_path: Path):
    """Case 10: <INV> crossing the gene body on hap1. Acceptance as with case 11 —
    pipeline runs without crashing and emits valid JSON. Detailed 'bisected'
    semantics depend on Liftoff's partial-lift scoring; covered as a behavior
    regression rather than a ground-truth assertion.
    """
    case_dir = paper_cases_dir / "10_inversion_body"
    out_dir = tmp_path / "results"
    repo_root = Path(__file__).parent.parent
    result = subprocess.run(
        [
            "uv", "run", "main.py", "analyze",
            "--gene", "G1", "--sample", "S1",
            "--vcf", str(case_dir / "phased.vcf.gz"),
            "--reference", str(case_dir / "reference.fa"),
            "--gff", str(case_dir / "annotation.gff3"),
            "--output-dir", str(out_dir),
        ],
        cwd=repo_root, capture_output=True, text=True, timeout=180,
    )
    assert result.returncode == 0, (
        f"inversion case crashed:\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    data = json.loads((out_dir / "S1_G1_analysis.json").read_text())
    assert "evidence" in data and "presence" in data["evidence"]
