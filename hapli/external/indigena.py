"""
INDIGENA wrapper — inductive disease–gene scoring via UPheno-projected
knowledge-graph embeddings.

INDIGENA (Zhapa-Camacho & Hoehndorf, *Bioinformatics* under review) trains
a PyKEEN knowledge-graph embedding (TransD on Graph 4 inductive in the
production config) over UPheno + MGI gene–phenotype + HPO disease–phenotype.
At inference, the embedding plus a Best-Match-Average (BMA) aggregator
scores any (HPO patient phenotype set, candidate gene) pair — including
unseen diseases (sets of phenotypes that were not present at training time).

This wrapper is the inference-only path. The KGE training is a one-time
cost that produces an artifact bundle:

    resources/indigena/v1/
        entity_embeddings.pt      tensor [num_entities, dim]
        entity_to_id.json         {entity_uri: int_idx}
        gene2pheno.json           {mgi_gene_id: [pheno_uri, ...]}
        disease2pheno.json        {disease_id: [pheno_uri, ...]}
        metadata.json             training hyperparams + reproduction metrics

The bundle is gitignored (~340 MB) and distributed separately (Zenodo for
the paper). Generated via INDIGENA's `kge_transd.py` + the
`extract_for_hapli.py` companion script.

References
----------
INDIGENA repository: https://github.com/bio-ontology-research-group/indigena
The `compare_vectorized` core (BMA / BMM aggregation) is vendored from
INDIGENA's `evaluation.py` (commit on file when training was run).
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Literal, Optional


class IndigenaNotAvailable(RuntimeError):
    """Raised when the INDIGENA artifact bundle is missing or torch isn't installed."""


@dataclass
class IndigenaScore:
    """One (gene, score) row in the INDIGENA ranking."""
    gene_id: str            # MGI gene ID (caller maps to HGNC via HCOP)
    score: float            # similarity in [0, 1]
    rank: int               # within this query (1 = top)


@dataclass
class IndigenaResult:
    """Full ranking result from one INDIGENA query."""
    query_phenotypes: list[str]
    scores: list[IndigenaScore]                 # ranked desc by score, all eligible genes
    aggregator: str                             # "BMA" or "BMM"
    metadata: dict = field(default_factory=dict)

    def top(self, k: int = 10) -> list[IndigenaScore]:
        return self.scores[:k]

    def rank_of(self, gene_id: str) -> Optional[int]:
        """Return 1-based rank of `gene_id` in the result, or None if absent."""
        for s in self.scores:
            if s.gene_id == gene_id:
                return s.rank
        return None


@dataclass
class _IndigenaModel:
    """Loaded INDIGENA model bundle. Lazily constructed; pass around as a
    handle so repeat queries share the same loaded embeddings."""
    entity_embeddings: "Any"             # torch.Tensor [num_entities, dim]
    entity_to_id: dict[str, int]
    gene2pheno: dict[str, list[str]]
    disease2pheno: dict[str, list[str]]
    metadata: dict


@lru_cache(maxsize=4)
def load_model(model_dir: Path) -> _IndigenaModel:
    """Load the INDIGENA bundle from disk. Cached by `model_dir`.

    Raises IndigenaNotAvailable if any required file is missing or torch
    is not importable (training-side dependency).
    """
    try:
        import torch as th
    except ImportError as exc:
        raise IndigenaNotAvailable(
            "torch is not installed. INDIGENA inference needs torch — "
            "install it (CPU is fine) or skip --with-indigena."
        ) from exc
    model_dir = Path(model_dir)
    required = ["entity_embeddings.pt", "entity_to_id.json",
                "gene2pheno.json", "disease2pheno.json"]
    missing = [f for f in required if not (model_dir / f).exists()]
    if missing:
        raise IndigenaNotAvailable(
            f"INDIGENA bundle at {model_dir} is missing files: {missing}. "
            f"Generate via `python extract_for_hapli.py` after training "
            f"`kge_transd.py --fold 0 --mode inductive --graph4 ...` from "
            f"the indigena repo."
        )
    embeddings = th.load(model_dir / "entity_embeddings.pt", weights_only=True)
    with (model_dir / "entity_to_id.json").open() as f:
        entity_to_id = json.load(f)
    with (model_dir / "gene2pheno.json").open() as f:
        gene2pheno = json.load(f)
    with (model_dir / "disease2pheno.json").open() as f:
        disease2pheno = json.load(f)
    metadata = {}
    if (model_dir / "metadata.json").exists():
        with (model_dir / "metadata.json").open() as f:
            metadata = json.load(f)
    return _IndigenaModel(
        entity_embeddings=embeddings,
        entity_to_id=entity_to_id,
        gene2pheno=gene2pheno,
        disease2pheno=disease2pheno,
        metadata=metadata,
    )


def _compare_vectorized(
    all_genes_pheno_vectors,                # [num_genes, max_phenos, emb_dim]
    disease_phenos_vectors,                 # [num_disease_phenos, emb_dim]
    gene_pheno_counts,                      # [num_genes] — real (non-padded) phenotype counts
    aggregator: Literal["BMA", "BMM"] = "BMA",
):
    """BMA / BMM aggregation over per-phenotype embeddings.

    Vendored verbatim from INDIGENA's `evaluation.py::compare_vectorized`
    with cosmetic edits (function name, parameter naming) and a default
    aggregator. The mathematical content is unchanged.

    Original attribution: Zhapa-Camacho, F. and Hoehndorf, R., INDIGENA,
    Bioinformatics (under review).
    https://github.com/bio-ontology-research-group/indigena (evaluation.py)

    Returns a tensor of shape [num_genes] with similarity scores in [0, 1].
    """
    import torch as th

    num_genes, max_phenos, emb_dim = all_genes_pheno_vectors.shape

    # Pairwise sim_matrix = sigmoid(<g_pheno_emb, d_pheno_emb>)
    sim_matrix = th.matmul(
        all_genes_pheno_vectors.view(-1, emb_dim),
        disease_phenos_vectors.T,
    )
    # Padding rows are exactly zero — exclude them from the sigmoid by
    # pushing them to -inf so sigmoid drops to 0.
    sim_matrix[sim_matrix == 0] = -th.inf
    sim_matrix = th.sigmoid(sim_matrix).view(num_genes, max_phenos, -1)

    # Gene-centric: for each gene phenotype, take max over disease phenotypes,
    # then average over the gene's real phenotypes.
    gene_max_sim, _ = sim_matrix.max(dim=-1)
    gene_centric_sum = gene_max_sim.sum(dim=-1)
    safe_counts = th.max(gene_pheno_counts, th.tensor(1.0))
    gene_centric_scores = gene_centric_sum / safe_counts

    # Disease-centric: symmetric — for each disease phenotype, take max over
    # gene phenotypes, then average.
    disease_max_sim, _ = sim_matrix.max(dim=1)
    disease_centric_scores = disease_max_sim.mean(dim=-1)

    if aggregator == "BMA":
        return (gene_centric_scores + disease_centric_scores) / 2
    if aggregator == "BMM":
        return th.max(gene_centric_scores, disease_centric_scores)
    raise ValueError(f"Unknown aggregator {aggregator!r} (use BMA or BMM)")


def score_genes(
    phenotypes: list[str],
    *,
    model_dir: Path | str,
    candidate_genes: Optional[list[str]] = None,
    aggregator: Literal["BMA", "BMM"] = "BMA",
    logger: Optional[logging.Logger] = None,
) -> IndigenaResult:
    """Rank candidate genes against a patient HPO phenotype set.

    Parameters
    ----------
    phenotypes : list of HPO term URIs (e.g. ``http://purl.obolibrary.org/obo/HP_0001250``).
                 Bare ``HP:0001250`` IDs are auto-converted.
    model_dir  : path to the INDIGENA artifact bundle (resources/indigena/v1/).
    candidate_genes : restrict scoring to these MGI gene IDs (e.g. an ortholog
                      panel for human gene priorities). If None, all genes in
                      gene2pheno are scored.
    aggregator : "BMA" (paper default) or "BMM".

    Returns IndigenaResult with scores ranked descending.

    Raises
    ------
    IndigenaNotAvailable : if torch missing or the bundle is incomplete.
    ValueError           : if no query phenotype is in the trained ontology.
    """
    log = logger or logging.getLogger(__name__)
    import torch as th

    model = load_model(Path(model_dir))
    e2id = model.entity_to_id
    embeddings = model.entity_embeddings

    # Normalize HPO IDs ("HP:0001250" → full URI used at training time)
    normalized = [_normalize_hpo(p) for p in phenotypes]
    pheno_ids = [e2id[p] for p in normalized if p in e2id]
    if not pheno_ids:
        raise ValueError(
            f"None of the {len(phenotypes)} query phenotypes are in the trained "
            f"ontology. First few queries: {phenotypes[:3]}; first few normalized: {normalized[:3]}"
        )
    if len(pheno_ids) < len(phenotypes):
        log.warning(
            "Dropped %d/%d query phenotypes not in the trained ontology",
            len(phenotypes) - len(pheno_ids), len(phenotypes),
        )

    # Build the gene phenotype matrix
    eligible_genes = list(model.gene2pheno.keys()) if candidate_genes is None else \
        [g for g in candidate_genes if g in model.gene2pheno]
    if not eligible_genes:
        raise ValueError(
            f"No candidate genes are scorable. Of {len(candidate_genes) if candidate_genes else 0} "
            f"requested genes, none are in the model's gene2pheno (size={len(model.gene2pheno)}). "
            f"Map HGNC -> MGI via HCOP first if querying human genes."
        )
    max_phenos = max(len(model.gene2pheno[g]) for g in eligible_genes)
    emb_dim = embeddings.shape[1]
    gene_matrix = th.zeros(len(eligible_genes), max_phenos, emb_dim)
    counts = th.zeros(len(eligible_genes), dtype=th.float32)
    dropped = 0
    for i, g in enumerate(eligible_genes):
        gp_ids = [e2id[p] for p in model.gene2pheno[g] if p in e2id]
        if not gp_ids:
            dropped += 1
            continue
        gene_matrix[i, :len(gp_ids)] = embeddings[th.tensor(gp_ids)]
        counts[i] = float(len(gp_ids))
    if dropped:
        log.info("Skipped %d/%d genes whose phenotypes are absent from the trained ontology",
                 dropped, len(eligible_genes))

    disease_vectors = embeddings[th.tensor(pheno_ids)]
    raw_scores = _compare_vectorized(gene_matrix, disease_vectors, counts, aggregator=aggregator)

    # Build ranked list (descending; genes with count==0 get score==0 → tail)
    score_pairs = sorted(
        zip(eligible_genes, raw_scores.tolist()),
        key=lambda p: p[1],
        reverse=True,
    )
    scores = [
        IndigenaScore(gene_id=g, score=float(s), rank=rk)
        for rk, (g, s) in enumerate(score_pairs, start=1)
    ]
    return IndigenaResult(
        query_phenotypes=phenotypes,
        scores=scores,
        aggregator=aggregator,
        metadata={
            "model": model.metadata.get("training", {}),
            "n_query_pheno_in_ontology": len(pheno_ids),
            "n_query_pheno_dropped": len(phenotypes) - len(pheno_ids),
            "n_genes_scored": len(eligible_genes) - dropped,
            "n_genes_dropped": dropped,
        },
    )


def _normalize_hpo(term: str) -> str:
    """HP:NNNNNNN → http://purl.obolibrary.org/obo/HP_NNNNNNN. Pass-through if
    already a URI or non-HPO entity ID (gene URIs / OMIM disease URIs etc.)."""
    if term.startswith("HP:"):
        return f"http://purl.obolibrary.org/obo/HP_{term.split(':', 1)[1]}"
    return term
