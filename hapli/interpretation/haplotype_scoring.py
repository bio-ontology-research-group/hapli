"""
Unified haplotype-level scoring — single / double / N-way missense, nonsense,
frameshift, in-frame indel, delins, and SV-shaped variants.

Given a reference protein and an ordered list of HGVS-p variants defining one
haplotype, produce a structured `HaplotypeScore`:

  * category            — what kind of haplotype this is (missense_single,
                          missense_multi, indel, lof, synonymous, unmodified, mixed)
  * hap_length          — length of the constructed haplotype protein
  * identity            — fraction of ref residues preserved (with alignment)
  * premature_stop_at   — 1-based aa position of the earliest in-frame stop
  * s_additive          — Σ log p(alt_k | ref_context) − log p(ref_k | ref_context)
                          across *missense* variants only (the per-variant baseline)
  * s_joint             — PLL(hap) − PLL(ref), equal-length substitution-only haplotypes
  * residual            — s_joint − s_additive where both are defined, else None

Non-missense variants (nonsense, frameshift, indel, delins) are NOT given an
ESM2 log-odds — there's no amino-acid token to score — but they still land in
the output with `category=lof|indel` so downstream code can filter/aggregate
uniformly. This is the Phase-5 "works for any variant shape" contract.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

from ..core.hgvs import HaplotypeProtein, ProteinVariant, VariantKind, apply_variants


@dataclass
class HaplotypeScore:
    """One row of output per haplotype."""

    # inputs
    variants_raw: str                                       # original HGVS string
    variants: list[ProteinVariant]
    # haplotype construction
    hap_seq: str
    category: str                                           # overall shape flag
    hap_length: int
    length_changed: bool
    premature_stop_at: Optional[int]
    # scoring
    identity: float
    n_missense: int
    n_nonsense: int
    n_frameshift: int
    n_indel: int
    n_synonymous: int
    n_skipped: int
    s_additive: Optional[float] = None
    s_joint: Optional[float] = None
    residual: Optional[float] = None
    flagged: bool = False


def _count_identity(ref_seq: str, hap_seq: str) -> float:
    """Fraction of reference residues retained in the haplotype, counted by
    position when lengths match; falls back to longest-common-prefix/suffix
    matching for length-changing cases (a stable, cheap proxy — the rigorous
    protein-level diff lives in hapli.core.protein_diff)."""
    if not ref_seq:
        return 1.0 if not hap_seq else 0.0
    if len(ref_seq) == len(hap_seq):
        matches = sum(1 for a, b in zip(ref_seq, hap_seq) if a == b)
        return matches / len(ref_seq)
    # Length differs — approximate as (longest common prefix + suffix) / ref_len.
    lcp = 0
    for a, b in zip(ref_seq, hap_seq):
        if a != b:
            break
        lcp += 1
    lcs = 0
    for a, b in zip(reversed(ref_seq[lcp:]), reversed(hap_seq[lcp:])):
        if a != b:
            break
        lcs += 1
    matches = lcp + lcs
    return min(1.0, matches / len(ref_seq))


def score_haplotype(
    ref_seq: str,
    variants_raw: str,
    variants: Iterable[ProteinVariant],
    *,
    scorer=None,                                            # type: ESM2Scorer | None
    min_residual_variants: int = 2,
    residual_threshold: float = 3.0,
    logger: Optional[logging.Logger] = None,
) -> HaplotypeScore:
    """Construct the haplotype protein from `variants` applied to `ref_seq`,
    then compute the full suite of shape-appropriate signals.

    The ESM residual is only computed when ALL of:
      * `scorer` is provided,
      * haplotype is missense-only (no indels / LoF / skipped),
      * haplotype length matches the reference,
      * at least `min_residual_variants` applied substitutions.
    """
    logger = logger or logging.getLogger(__name__)
    vlist = list(variants)

    hp: HaplotypeProtein = apply_variants(ref_seq, vlist)

    # Classify counts by variant kind (among applied).
    def _count(kind: VariantKind) -> int:
        return sum(1 for v in hp.applied if v.kind == kind)

    n_missense = _count(VariantKind.MISSENSE)
    n_nonsense = _count(VariantKind.NONSENSE)
    n_frameshift = _count(VariantKind.FRAMESHIFT)
    n_indel = _count(VariantKind.DEL) + _count(VariantKind.DELINS)
    n_synonymous = _count(VariantKind.SYNONYMOUS)
    n_skipped = len(hp.skipped)

    identity = _count_identity(ref_seq, hp.seq)

    score = HaplotypeScore(
        variants_raw=variants_raw,
        variants=vlist,
        hap_seq=hp.seq,
        category=hp.overall_category,
        hap_length=len(hp.seq),
        length_changed=hp.length_changed,
        premature_stop_at=hp.premature_stop_at,
        identity=identity,
        n_missense=n_missense,
        n_nonsense=n_nonsense,
        n_frameshift=n_frameshift,
        n_indel=n_indel,
        n_synonymous=n_synonymous,
        n_skipped=n_skipped,
    )

    # ESM2 residual — only for missense-only, equal-length, multi-variant haplotypes.
    if (
        scorer is not None
        and hp.overall_category == "missense_multi"
        and not hp.length_changed
        and n_missense >= min_residual_variants
    ):
        from .epistasis import compute_residuals_for_diff
        r = compute_residuals_for_diff(
            scorer,
            transcript_id="haplotype",
            haplotype=1,
            ref_seq=ref_seq,
            hap_seq=hp.seq,
            min_substitutions=min_residual_variants,
            logger=logger,
        )
        score.s_additive = r.s_additive
        score.s_joint = r.s_joint
        score.residual = r.residual
        score.flagged = r.flagged

    return score


def score_haplotypes_batch(
    ref_seq: str,
    haplotypes: Sequence[tuple[str, Iterable[ProteinVariant]]],
    *,
    scorer=None,
    logger: Optional[logging.Logger] = None,
    **kwargs,
) -> list[HaplotypeScore]:
    """Score many haplotypes against one reference.

    Each element of `haplotypes` is `(raw_hgvs_string, parsed_variants)`.
    Returns one HaplotypeScore per input haplotype in the same order.
    """
    logger = logger or logging.getLogger(__name__)
    return [
        score_haplotype(ref_seq, raw, vs, scorer=scorer, logger=logger, **kwargs)
        for raw, vs in haplotypes
    ]


__all__ = [
    "HaplotypeScore",
    "score_haplotype",
    "score_haplotypes_batch",
]
