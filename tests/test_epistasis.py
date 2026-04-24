"""
Tests for the epistasis residual — the Phase 3 research contribution.

The heavy end-to-end tests are guarded behind an `esm` marker and auto-skip
when fair-esm/torch aren't installed or the ESM2 checkpoint isn't cached.
The pure-Python helpers (substitution_positions, decision rules, etc.) are
unit-tested without touching the model.
"""

from __future__ import annotations

import importlib.util

import pytest

from hapli.interpretation.epistasis import (
    EpistasisInput,
    compute_residual,
    compute_residuals_for_diff,
    substitution_positions,
)


HAS_ESM = (
    importlib.util.find_spec("torch") is not None
    and importlib.util.find_spec("esm") is not None
)
needs_esm = pytest.mark.skipif(not HAS_ESM, reason="fair-esm/torch not installed")


# ─────────────────────────────────────────────────────────────────────────────
# Pure-python helpers (no model)
# ─────────────────────────────────────────────────────────────────────────────
def test_substitution_positions_basic():
    assert substitution_positions("MAAAA", "MAVAA") == [3]
    assert substitution_positions("MAAAA", "MAAAA") == []
    assert substitution_positions("MAAAA", "MAVBA") == [3, 4]


def test_substitution_positions_rejects_length_mismatch():
    with pytest.raises(ValueError):
        substitution_positions("MAAA", "MAAAA")


# ─────────────────────────────────────────────────────────────────────────────
# Model-backed tests. All share a session-scoped scorer so the ESM2 checkpoint
# loads once.
# ─────────────────────────────────────────────────────────────────────────────
@pytest.fixture(scope="session")
def scorer():
    if not HAS_ESM:
        pytest.skip("fair-esm/torch not installed")
    from hapli.interpretation.esm_scoring import ESM2Scorer
    return ESM2Scorer()                              # default: tiny 8M model, CPU


REF = "MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"   # 40 aa


@needs_esm
def test_frameshift_rescue_window_has_large_positive_residual(scorer):
    """Hero case: 6 residues in a contiguous rescue window.

    Each variant scored in reference context looks catastrophic (ESM's prior
    for an R/C in a poly-Ala context is tiny). Jointly, the cluster is much
    more plausible. Residual must be large + positive.
    """
    hap = "MAAARCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    r = compute_residuals_for_diff(scorer, "T1", 1, REF, hap)
    assert r.n_variants == 6
    assert r.s_joint > r.s_additive, "joint should be less deleterious than additive"
    assert r.residual > 20.0, f"residual unexpectedly small: {r.residual}"
    assert r.flagged is True


@needs_esm
def test_single_substitution_does_not_compute_residual(scorer):
    hap = "MAVAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    r = compute_residuals_for_diff(scorer, "T1", 1, REF, hap)
    assert r.n_variants == 1
    assert r.flagged is False
    assert r.residual == 0.0    # untouched default


@needs_esm
def test_zero_substitutions_is_trivially_zero(scorer):
    r = compute_residuals_for_diff(scorer, "T1", 1, REF, REF)
    assert r.n_variants == 0
    assert r.flagged is False
    assert r.s_joint == 0.0


@needs_esm
def test_length_mismatch_is_skipped_not_errored(scorer):
    hap = REF[:20]
    r = compute_residuals_for_diff(scorer, "T1", 1, REF, hap)
    assert r.n_variants == 0
    assert r.flagged is False


@needs_esm
def test_scorer_caching_returns_identical_values_on_second_call(scorer):
    """Cached per-sequence log-probs should be byte-identical across calls."""
    seq = "MACDEFGHIKLMNPQRSTVWY"
    a = scorer.per_position_log_probs(seq)
    b = scorer.per_position_log_probs(seq)
    assert (a == b).all()


@needs_esm
def test_additive_score_is_strongly_negative_for_rare_substitutions(scorer):
    """Sanity check: for the 6-rare-residue cluster, the additive score alone
    (without the joint rescue) should look *very* deleterious — this is
    precisely the failure mode of variant-by-variant scoring that the joint
    term rescues in the full residual."""
    hap = "MAAARCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    r = compute_residuals_for_diff(scorer, "T1", 1, REF, hap)
    # Additive alone says "very deleterious" (dozens of log-units negative).
    assert r.s_additive < -20.0
    # Joint alone says "mildly deleterious" (under 10 log-units negative).
    assert r.s_joint > -10.0


@needs_esm
def test_explicit_input_accepted_by_compute_residual(scorer):
    hap = "MAAARCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    positions = substitution_positions(REF, hap)
    r = compute_residual(
        scorer,
        EpistasisInput(
            transcript="T1",
            haplotype=1,
            ref_seq=REF,
            hap_seq=hap,
            substitution_positions=positions,
        ),
    )
    assert r.n_variants == 6
    assert r.flagged is True
