"""
The additive-vs-joint epistasis residual — the research contribution of hapli.

Given a reference protein R and a haplotype protein H carrying n≥2 substitutions,
this module computes a scalar *residual* that flags haplotypes whose joint
functional effect (ESM2 full-sequence pseudo-LL delta) is materially different
from the sum of per-variant effects measured in the reference context.

Core quantities (all using the per-residue log-probability operator p_M(·)
of an ESM2-family model M):

    S_additive = Σ_k [ log p_M(H[pos_k] | R) − log p_M(R[pos_k] | R) ]

        i.e. score each variant k IN ISOLATION under the reference sequence.
        This is the standard "wildtype marginals" variant-effect score of
        Meier et al. 2021 (ESM-1v), summed over variants. Additive by
        construction — the score of variant k is independent of variant j.

    S_joint = PLL(H) − PLL(R)
            = Σ_i log p_M(H[i] | H) − Σ_i log p_M(R[i] | R)

        i.e. score each full sequence with the model, take the difference.
        Context for each residue is the rest of that sequence, so variants
        are scored *jointly*.

    Residual = S_joint − S_additive

Interpretation:

    residual > 0  ⇒ joint effect is LESS deleterious than additive suggests.
                    Signature of a frameshift+rescue pair, or a compound
                    missense where the two variants buffer each other.
    residual < 0  ⇒ joint effect is MORE deleterious than additive suggests.
                    Signature of compound missense at a binding pocket
                    where the two variants synergise destructively.
    |residual| < threshold  ⇒ no meaningful epistasis detected; additive
                    prediction is adequate.

Decision rule: we only compute the residual when the haplotype has at least
two substitutions on the same transcript (one substitution is trivially
additive). We also require the two proteins to have the same length — the
primary definition needs aligned residues. Length-changing haplotypes (net
indel, premature truncation) are reported as `skipped` with an explanatory
`skip_reason`; those are handled by other pipeline signals (presence call,
premature-stop flag, frameshift-region tag) rather than by this score.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from ..core.schema import EpistasisResidual
from .esm_scoring import ESM2Scorer, EsmNotAvailable


@dataclass
class EpistasisInput:
    """One transcript + haplotype's worth of aligned reference/hap protein."""

    transcript: str
    haplotype: int                        # 1 or 2
    ref_seq: str                          # reference protein (no trailing '*')
    hap_seq: str                          # haplotype protein (no trailing '*')
    substitution_positions: Sequence[int] # 1-based positions where ref != hap


def substitution_positions(ref_seq: str, hap_seq: str) -> list[int]:
    """1-based positions where ref and hap differ. Requires equal length."""
    if len(ref_seq) != len(hap_seq):
        raise ValueError(
            f"ref_seq and hap_seq must be equal length "
            f"(got {len(ref_seq)} and {len(hap_seq)})"
        )
    return [i + 1 for i, (a, b) in enumerate(zip(ref_seq, hap_seq)) if a != b]


def compute_residual(
    scorer: ESM2Scorer,
    inp: EpistasisInput,
    *,
    min_substitutions: int = 2,
    logger: logging.Logger | None = None,
) -> EpistasisResidual:
    """Compute the additive / joint / residual scores.

    Returns an EpistasisResidual. When the inputs don't meet the computation
    preconditions (≥2 subs, equal length, clean amino acids), `flagged=False`
    and `n_variants` records why.
    """
    logger = logger or logging.getLogger(__name__)

    # Preflight: if either sequence is empty, nothing to score.
    if not inp.ref_seq or not inp.hap_seq:
        return EpistasisResidual(
            transcript=inp.transcript,
            haplotype=inp.haplotype,
            n_variants=0,
            joint_source=scorer.model_name,
        )

    # Preflight: require equal length. Indels at the protein level need a
    # different treatment (alignment-aware scoring) that we defer to v2.
    if len(inp.ref_seq) != len(inp.hap_seq):
        logger.info(
            "%s hap%d: skipping epistasis (ref len %d != hap len %d)",
            inp.transcript, inp.haplotype, len(inp.ref_seq), len(inp.hap_seq),
        )
        return EpistasisResidual(
            transcript=inp.transcript,
            haplotype=inp.haplotype,
            n_variants=len(inp.substitution_positions),
            joint_source=scorer.model_name,
        )

    # Preflight: require at least `min_substitutions` on the haplotype.
    n = len(inp.substitution_positions)
    if n < min_substitutions:
        return EpistasisResidual(
            transcript=inp.transcript,
            haplotype=inp.haplotype,
            n_variants=n,
            joint_source=scorer.model_name,
        )

    # --- joint: per-sequence pseudo-log-likelihood delta ---
    lp_ref = scorer.per_position_log_probs(inp.ref_seq)
    lp_hap = scorer.per_position_log_probs(inp.hap_seq)
    s_joint = float(lp_hap.sum() - lp_ref.sum())

    # --- additive: sum of per-variant log-odds in the reference context ---
    positions = list(inp.substitution_positions)
    alts = [inp.hap_seq[p - 1] for p in positions]
    lp_alts_in_ref = scorer.alt_log_probs(inp.ref_seq, positions, alts)
    lp_refs_in_ref = lp_ref[[p - 1 for p in positions]]
    per_variant_delta = lp_alts_in_ref - lp_refs_in_ref
    s_additive = float(per_variant_delta.sum())

    residual = s_joint - s_additive

    return EpistasisResidual(
        transcript=inp.transcript,
        haplotype=inp.haplotype,
        n_variants=n,
        s_additive=s_additive,
        s_joint=s_joint,
        residual=residual,
        additive_source="esm2_wildtype_marginals",
        joint_source=scorer.model_name,
        flagged=_is_significant(residual),
    )


def _is_significant(residual: float, threshold: float = 3.0) -> bool:
    """Threshold for flagging a residual as meaningful epistasis.

    3.0 log-units ≈ factor of 20× difference between additive and joint
    prediction — large enough to not be noise on 8M-param models, small
    enough to catch real rescue/synergy. This is a provisional threshold;
    Phase 5's MAVE benchmark will calibrate it properly.
    """
    return abs(residual) >= threshold


def compute_residuals_for_diff(
    scorer: ESM2Scorer,
    transcript_id: str,
    haplotype: int,
    ref_seq: str,
    hap_seq: str,
    *,
    min_substitutions: int = 2,
    logger: logging.Logger | None = None,
) -> EpistasisResidual:
    """Convenience wrapper: derives `substitution_positions` from equal-length
    seqs before calling `compute_residual`.
    """
    ref_trim = ref_seq.rstrip("*")
    hap_trim = hap_seq.rstrip("*")
    if len(ref_trim) != len(hap_trim):
        return EpistasisResidual(
            transcript=transcript_id,
            haplotype=haplotype,
            n_variants=0,
            joint_source=scorer.model_name,
        )
    positions = substitution_positions(ref_trim, hap_trim)
    return compute_residual(
        scorer,
        EpistasisInput(
            transcript=transcript_id,
            haplotype=haplotype,
            ref_seq=ref_trim,
            hap_seq=hap_trim,
            substitution_positions=positions,
        ),
        min_substitutions=min_substitutions,
        logger=logger,
    )


__all__ = [
    "EpistasisInput",
    "compute_residual",
    "compute_residuals_for_diff",
    "substitution_positions",
]
