"""
Diploid report assembly — the Phase 4 aggregator.

Collapses the per-gene Evidence bundle into a *two-number* summary that
deliberately refuses to commit to a single "diploid call". Downstream
clinical pipelines decide how to combine the two numbers with constraint
priors (pLI, haploinsufficiency, etc.) according to their own disease
model — a choice that belongs to the clinical pipeline, not to hapli.

Two numbers:

  hap1_score, hap2_score ∈ [0, 1]  — "function preservation" on each haplotype.
    1.0 : gene present, ORF intact, protein matches reference.
    0.0 : gene absent, or protein severely truncated at the N-terminus.
    intermediate values : length-preserving substitutions + late PTCs scale
                          linearly with alignment identity / stop-codon position.

The two-number output is complemented by:

  min_score, max_score   — convenience; trivial to recompute but saves consumers work.
  compound_het_lof       — both haplotypes independently called LoF (score below
                          `LOF_THRESHOLD`). Surfaces the autosomal-recessive
                          compound-het pattern even when the actual variants differ.
  constraints            — gnomAD / ClinGen priors for this gene (may be empty).

Decision rule helpers are exported so downstream code can apply its own
thresholds consistently:

  - score_haplotype(presence, protein_diff) → Optional[float]
  - is_lof(score, threshold)               → bool
"""

from __future__ import annotations

import logging
from typing import Optional

from ..core.schema import (
    ConstraintPriors,
    DiploidReport,
    GeneEvidence,
    PresenceCall,
    ProteinDiff,
)


LOF_THRESHOLD: float = 0.5
"""Haplotypes with score ≤ this are treated as LoF for compound-het detection."""


# ─────────────────────────────────────────────────────────────────────────────
# Scoring
# ─────────────────────────────────────────────────────────────────────────────
def score_haplotype(
    presence: PresenceCall | None,
    protein_diffs: list[ProteinDiff] | None = None,
    *,
    nterm_ptc_fraction: float = 0.10,
) -> Optional[float]:
    """Return a [0, 1] "function preservation" score for a single haplotype.

    The formula is intentionally simple and transparent so clinical consumers
    can audit it:

      * presence=deleted                     → 0.0                 (gene gone)
      * presence=not_run / no diffs          → None                (undetermined)
      * any transcript has N-terminal PTC    → 0.0                 (LoF)
      * any transcript has mid-protein PTC   → identity * (ptc_pos / ref_len)
      * otherwise                            → min transcript identity across isoforms

    A haplotype with multiple transcripts takes the MIN identity — the gene
    is only as functional as its worst transcript (conservative).
    """
    if presence is None or presence.status == "not_run":
        return None
    if presence.status == "deleted":
        return 0.0
    if not protein_diffs:
        # Lifted with some QC signal but no protein diff: fall back to presence.
        if presence.status == "intact":
            return 1.0
        if presence.status == "duplicated":
            return 1.0
        if presence.status in ("low_identity", "partial", "uncertain"):
            # Conservative middle — presence is known, function is not.
            return presence.sequence_identity or 0.5
        return None

    per_transcript_scores: list[float] = []
    for diff in protein_diffs:
        ref_len = diff.ref_length or 1
        pts = diff.premature_stop_at
        if pts is not None:
            frac = pts / ref_len
            if frac <= nterm_ptc_fraction:
                # Very early PTC → effectively no protein → LoF.
                per_transcript_scores.append(0.0)
            else:
                # Mid/late PTC: scale identity by how much of the protein is preserved.
                per_transcript_scores.append(diff.identity * frac)
        else:
            per_transcript_scores.append(diff.identity)

    protein_score = min(per_transcript_scores)

    # Cap by presence quality: Liftoff may report a degraded presence call
    # (low_identity, partial, uncertain) even when the protein extracted from
    # the lift happens to match reference — e.g. an inversion that Liftoff
    # mapped only over its unaffected flanks. In those cases the protein-only
    # score misses the structural signal. Take min(protein_score, presence_id)
    # so the structural disruption dampens the final call.
    if presence.status in ("low_identity", "partial"):
        cap = presence.sequence_identity
        if cap is not None:
            protein_score = min(protein_score, cap)
    return max(0.0, protein_score)


def is_lof(score: Optional[float], threshold: float = LOF_THRESHOLD) -> bool:
    """True iff the haplotype is called loss-of-function at the given threshold."""
    return score is not None and score <= threshold


# ─────────────────────────────────────────────────────────────────────────────
# Assembly
# ─────────────────────────────────────────────────────────────────────────────
def build_diploid_report(
    evidence: GeneEvidence,
    constraints: ConstraintPriors | None = None,
    *,
    lof_threshold: float = LOF_THRESHOLD,
    logger: logging.Logger | None = None,
) -> DiploidReport:
    """Collapse the per-gene Evidence into a DiploidReport.

    Does NOT commit to a single-number diploid call. Emits:
      - hap1_score, hap2_score (may be None if the stage wasn't run)
      - min_score, max_score (skipping Nones)
      - compound_het_lof: True iff both haplotypes independently call LoF
      - constraints: pLI / MisZ / o/e / ClinGen if provided

    Parameters
    ----------
    evidence     : the populated GeneEvidence bundle (post-Phase 2+).
    constraints  : optional ConstraintPriors from the gnomAD + ClinGen lookup.
    lof_threshold: score threshold for compound-het LoF detection.
    """
    logger = logger or logging.getLogger(__name__)

    diffs_by_hap: dict[int, list[ProteinDiff]] = {1: [], 2: []}
    for d in evidence.protein:
        diffs_by_hap.setdefault(d.haplotype, []).append(d)

    hap1_score = score_haplotype(evidence.presence.get("hap1"), diffs_by_hap.get(1))
    hap2_score = score_haplotype(evidence.presence.get("hap2"), diffs_by_hap.get(2))

    known = [s for s in (hap1_score, hap2_score) if s is not None]
    min_score = min(known) if known else None
    max_score = max(known) if known else None

    compound_het = (
        hap1_score is not None
        and hap2_score is not None
        and is_lof(hap1_score, lof_threshold)
        and is_lof(hap2_score, lof_threshold)
    )

    report = DiploidReport(
        hap1_score=hap1_score,
        hap2_score=hap2_score,
        min_score=min_score,
        max_score=max_score,
        compound_het_lof=compound_het,
        constraints=constraints or ConstraintPriors(),
    )
    logger.info(
        "%s diploid: hap1=%s hap2=%s min=%s max=%s compound_het_lof=%s",
        evidence.gene or "<gene>",
        f"{hap1_score:.3f}" if hap1_score is not None else "na",
        f"{hap2_score:.3f}" if hap2_score is not None else "na",
        f"{min_score:.3f}" if min_score is not None else "na",
        f"{max_score:.3f}" if max_score is not None else "na",
        compound_het,
    )
    return report


__all__ = [
    "LOF_THRESHOLD",
    "build_diploid_report",
    "is_lof",
    "score_haplotype",
]
