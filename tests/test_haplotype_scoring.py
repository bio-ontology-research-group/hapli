"""Tests for the unified haplotype scorer."""

from __future__ import annotations

import importlib.util

import pytest

from hapli.core.hgvs import parse_hgvs_pro
from hapli.interpretation.haplotype_scoring import (
    HaplotypeScore,
    score_haplotype,
    score_haplotypes_batch,
)


HAS_ESM = (
    importlib.util.find_spec("torch") is not None
    and importlib.util.find_spec("esm") is not None
)
needs_esm = pytest.mark.skipif(not HAS_ESM, reason="ml extra not installed")


REF = "MAAAAAAAAAARCCCCCAAAAAAAAAAAAAAAAAAAAAAAA"   # 40 aa


# ─────────────────────────────────────────────────────────────────────────────
# No ESM — pure shape routing
# ─────────────────────────────────────────────────────────────────────────────
def test_single_missense_without_esm():
    s = score_haplotype(REF, "p.Arg12Gly", parse_hgvs_pro("p.Arg12Gly"))
    assert s.category == "missense_single"
    assert s.n_missense == 1
    assert s.hap_length == len(REF)
    assert not s.length_changed
    assert s.residual is None


def test_double_missense_without_esm_still_categorises():
    s = score_haplotype(REF, "p.[Arg12Gly;Cys13Trp]",
                        parse_hgvs_pro("p.[Arg12Gly;Cys13Trp]"))
    assert s.category == "missense_multi"
    assert s.n_missense == 2
    assert s.residual is None      # no scorer → no residual


def test_nonsense_routes_to_lof_category():
    s = score_haplotype(REF, "p.Arg12Ter", parse_hgvs_pro("p.Arg12Ter"))
    assert s.category == "lof"
    assert s.n_nonsense == 1
    assert s.premature_stop_at == 12
    assert s.length_changed


def test_frameshift_routes_to_lof_category():
    s = score_haplotype(REF, "p.Cys13GlyfsTer5",
                        parse_hgvs_pro("p.Cys13GlyfsTer5"))
    assert s.category == "lof"
    assert s.n_frameshift == 1


def test_range_del_routes_to_indel_category():
    s = score_haplotype(REF, "p.Arg12_Cys17del",
                        parse_hgvs_pro("p.Arg12_Cys17del"))
    assert s.category == "indel"
    assert s.n_indel == 1
    assert s.length_changed


def test_nterminal_del_simulating_sv_is_flagged_indel():
    """Large N-terminal range deletion (SV-like). Haplotype is shorter; category=indel."""
    s = score_haplotype(REF, "p.Met1_Ala11del",
                        parse_hgvs_pro("p.Met1_Ala11del"))
    assert s.category == "indel"
    assert s.hap_length == len(REF) - 11


def test_synonymous_unchanged():
    s = score_haplotype(REF, "p.Arg12=", parse_hgvs_pro("p.Arg12="))
    assert s.category == "synonymous"
    assert s.hap_seq == REF


def test_skipped_variant_tracked():
    # Reference mismatch — the M1 position is 'M' not 'W'
    s = score_haplotype(REF, "p.Trp1Gly", parse_hgvs_pro("p.Trp1Gly"))
    assert s.n_skipped == 1
    assert s.category == "unmodified"


def test_batch_wrapper_preserves_order():
    raws = ["p.Arg12Gly", "p.Arg12Ter", "p.[Arg12Gly;Cys13Trp]"]
    parsed = [parse_hgvs_pro(r) for r in raws]
    out = score_haplotypes_batch(REF, list(zip(raws, parsed)))
    assert [s.category for s in out] == ["missense_single", "lof", "missense_multi"]


# ─────────────────────────────────────────────────────────────────────────────
# With ESM — residual path
# ─────────────────────────────────────────────────────────────────────────────
@pytest.fixture(scope="module")
def scorer():
    if not HAS_ESM:
        pytest.skip("ml extra not installed")
    from hapli.interpretation.esm_scoring import ESM2Scorer
    return ESM2Scorer()


@needs_esm
def test_double_missense_gets_residual(scorer):
    s = score_haplotype(
        REF, "p.[Arg12Gly;Cys13Trp]",
        parse_hgvs_pro("p.[Arg12Gly;Cys13Trp]"),
        scorer=scorer,
    )
    assert s.category == "missense_multi"
    assert s.s_additive is not None
    assert s.s_joint is not None
    assert s.residual is not None


@needs_esm
def test_single_missense_still_no_residual_even_with_scorer(scorer):
    s = score_haplotype(
        REF, "p.Arg12Gly", parse_hgvs_pro("p.Arg12Gly"),
        scorer=scorer,
    )
    assert s.category == "missense_single"
    assert s.residual is None         # need ≥2 substitutions


@needs_esm
def test_lof_haplotype_skips_residual_even_with_scorer(scorer):
    """Length-changing haplotypes don't get an ESM residual — it's undefined
    for them (one sequence is shorter). But the shape signals (category,
    premature_stop_at, length_changed, identity) are still populated."""
    s = score_haplotype(
        REF, "p.Arg12Ter", parse_hgvs_pro("p.Arg12Ter"),
        scorer=scorer,
    )
    assert s.category == "lof"
    assert s.residual is None
    assert s.length_changed
    assert s.premature_stop_at == 12


@needs_esm
def test_indel_haplotype_skips_residual_even_with_scorer(scorer):
    s = score_haplotype(
        REF, "p.Arg12_Cys17del", parse_hgvs_pro("p.Arg12_Cys17del"),
        scorer=scorer,
    )
    assert s.category == "indel"
    assert s.residual is None
    assert s.length_changed
