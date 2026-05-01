"""Tests for hapli/external/indigena.py.

Two layers:

1. **Pure-logic tests** (always run): validate _compare_vectorized math
   against hand-calculated synthetic embeddings + the HPO ID
   normalization helper.

2. **Live-model smoke tests** (skipped if bundle absent): load the
   actual INDIGENA bundle from `resources/indigena/v1/`, score a
   synthetic patient phenotype set against a small candidate gene
   panel, and assert the scoring shape + a sanity-check ranking.

Bundle is gitignored (~340 MB), distributed via Zenodo. Generate
locally via the indigena training pipeline (see indigena.py docstring).
"""

from __future__ import annotations

from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).parent.parent
BUNDLE_DIR = REPO_ROOT / "resources" / "indigena" / "v1"
HAS_BUNDLE = (BUNDLE_DIR / "entity_embeddings.pt").exists()


# ─────────────────────────────────────────────────────────────────────────
# Layer 1: pure-logic tests (no bundle, no torch dependency check)
# ─────────────────────────────────────────────────────────────────────────


def test_normalize_hpo_id_to_uri():
    from hapli.external.indigena import _normalize_hpo
    assert _normalize_hpo("HP:0001250") == "http://purl.obolibrary.org/obo/HP_0001250"


def test_normalize_hpo_passes_through_uri():
    from hapli.external.indigena import _normalize_hpo
    uri = "http://purl.obolibrary.org/obo/HP_0001263"
    assert _normalize_hpo(uri) == uri


def test_normalize_hpo_passes_through_non_hpo():
    """Gene/disease URIs pass through unchanged."""
    from hapli.external.indigena import _normalize_hpo
    assert _normalize_hpo("MGI:1234567") == "MGI:1234567"


def test_compare_vectorized_bma_perfect_match():
    """Two identical genes with identical phenotype set -> identical scores ~1.0
    when the gene phenotypes are the same as the disease phenotypes (sigmoid of
    a positive dot product saturates near 1)."""
    pytest.importorskip("torch")
    import torch as th
    from hapli.external.indigena import _compare_vectorized

    # Construct embeddings where pheno_a · pheno_a = +5 (sigmoid ~ 0.993)
    # and pheno_a · pheno_b = -5 (sigmoid ~ 0.007). Gene 0 has [pheno_a],
    # gene 1 has [pheno_b]. Disease has [pheno_a].
    pheno_a = th.tensor([1.0, 1.0, 1.0, 1.0, 1.0])
    pheno_b = -pheno_a
    # Pad to (num_genes=2, max_phenos=1, emb_dim=5)
    gene_matrix = th.stack([pheno_a, pheno_b]).unsqueeze(1)
    disease_vec = pheno_a.unsqueeze(0)
    counts = th.tensor([1.0, 1.0])

    scores = _compare_vectorized(gene_matrix, disease_vec, counts, aggregator="BMA")
    assert scores[0] > 0.99, f"identical-phenotype gene should score near 1, got {scores[0]}"
    assert scores[1] < 0.10, f"opposite-phenotype gene should score low, got {scores[1]}"


def test_compare_vectorized_bmm_returns_max():
    """BMM is max(gene-centric, disease-centric); BMA is the mean.
    On a single-phenotype query the two should differ when the gene-vs-disease
    centric scores are asymmetric — BMM can only be ≥ BMA."""
    pytest.importorskip("torch")
    import torch as th
    from hapli.external.indigena import _compare_vectorized

    # Asymmetric setup: gene has 4 phenotypes, only one matches the disease
    # → gene-centric average over 4 is dragged down; disease-centric stays high.
    match = th.tensor([2.0, 2.0, 2.0])
    miss = th.tensor([0.1, 0.1, 0.1])
    gene_matrix = th.stack([match, miss, miss, miss]).unsqueeze(0)  # (1, 4, 3)
    disease_vec = match.unsqueeze(0)
    counts = th.tensor([4.0])

    bma = _compare_vectorized(gene_matrix, disease_vec, counts, aggregator="BMA")
    bmm = _compare_vectorized(gene_matrix, disease_vec, counts, aggregator="BMM")
    assert bmm[0] >= bma[0], f"BMM ({bmm[0]}) should be >= BMA ({bma[0]})"


def test_compare_vectorized_unknown_aggregator_raises():
    pytest.importorskip("torch")
    import torch as th
    from hapli.external.indigena import _compare_vectorized

    g = th.zeros(1, 1, 3)
    d = th.zeros(1, 3)
    c = th.tensor([1.0])
    with pytest.raises(ValueError, match="Unknown aggregator"):
        _compare_vectorized(g, d, c, aggregator="GIC")  # noqa


def test_load_model_missing_bundle_raises(tmp_path):
    """If bundle dir lacks any required file, surface a helpful error."""
    pytest.importorskip("torch")
    from hapli.external.indigena import IndigenaNotAvailable, load_model

    empty = tmp_path / "indigena_empty"
    empty.mkdir()
    with pytest.raises(IndigenaNotAvailable, match="missing files"):
        load_model(empty)


# ─────────────────────────────────────────────────────────────────────────
# Layer 2: live-model smoke tests (skipped if bundle missing)
# ─────────────────────────────────────────────────────────────────────────


pytestmark_live = pytest.mark.skipif(
    not HAS_BUNDLE,
    reason=f"INDIGENA bundle missing at {BUNDLE_DIR}; run training + extraction"
)


@pytestmark_live
def test_score_genes_returns_full_ranking():
    """Score a small synthetic phenotype set against the full mouse gene panel.
    Sanity: every gene gets a rank, scores are bounded [0, 1], top-scoring
    gene is plausible (we don't pin a specific gene since training is
    fold-stochastic)."""
    from hapli.external.indigena import score_genes
    import json

    # Use a real disease's phenotypes as the query: pull the first OMIM disease
    # in the bundle's disease2pheno and use ITS phenotypes. The associated
    # gene(s) for that disease should rank competitively.
    with (BUNDLE_DIR / "disease2pheno.json").open() as f:
        disease2pheno = json.load(f)
    sample_disease, sample_phenos = next(iter(disease2pheno.items()))

    result = score_genes(
        phenotypes=sample_phenos[:20],   # cap at 20 to keep test fast
        model_dir=BUNDLE_DIR,
        aggregator="BMA",
    )
    assert len(result.scores) > 1000, f"expected to score >1000 mouse genes, got {len(result.scores)}"
    assert all(0.0 <= s.score <= 1.0 for s in result.scores)
    # Ranks must be 1..N strictly increasing
    for i, s in enumerate(result.scores, start=1):
        assert s.rank == i
    # Top score > median score (model is non-trivial)
    median = result.scores[len(result.scores) // 2].score
    assert result.scores[0].score > median


@pytestmark_live
def test_score_genes_candidate_panel_restricts():
    """Scoring with a candidate_genes panel must only score those genes."""
    from hapli.external.indigena import score_genes
    import json

    with (BUNDLE_DIR / "gene2pheno.json").open() as f:
        gene2pheno = json.load(f)
    panel = list(gene2pheno.keys())[:50]

    with (BUNDLE_DIR / "disease2pheno.json").open() as f:
        disease2pheno = json.load(f)
    _, phenos = next(iter(disease2pheno.items()))

    result = score_genes(
        phenotypes=phenos[:5],
        model_dir=BUNDLE_DIR,
        candidate_genes=panel,
    )
    assert len(result.scores) == 50
    assert {s.gene_id for s in result.scores} == set(panel)


@pytestmark_live
def test_score_genes_unknown_phenotype_raises():
    from hapli.external.indigena import score_genes

    with pytest.raises(ValueError, match="None of the .* query phenotypes"):
        score_genes(
            phenotypes=["HP:9999999", "HP:9999998"],   # implausible IDs
            model_dir=BUNDLE_DIR,
        )


@pytestmark_live
def test_score_genes_metadata_round_trips():
    """Metadata block exposes the training hyperparams + drop-counts."""
    from hapli.external.indigena import score_genes
    import json

    with (BUNDLE_DIR / "disease2pheno.json").open() as f:
        disease2pheno = json.load(f)
    _, phenos = next(iter(disease2pheno.items()))

    result = score_genes(
        phenotypes=phenos[:10],
        model_dir=BUNDLE_DIR,
    )
    assert "n_genes_scored" in result.metadata
    assert "n_query_pheno_in_ontology" in result.metadata
    # Training-config block should round-trip whatever extract_for_hapli.py wrote
    assert "method" in result.metadata.get("model", {})
