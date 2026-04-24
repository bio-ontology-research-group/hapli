"""Tests for hapli.core.hgvs — HGVS parser + haplotype-protein constructor."""

from __future__ import annotations

import pytest

from hapli.core.hgvs import (
    HaplotypeProtein,
    ProteinVariant,
    VariantKind,
    apply_variants,
    parse_hgvs_pro,
)


REF = "MAAAAAAAAAARCCCCCAAAAAAAAAAAAAAAAAAAAAAAA"   # 40 aa


# ─────────────────────────────────────────────────────────────────────────────
# parse_hgvs_pro — every supported shape
# ─────────────────────────────────────────────────────────────────────────────
def test_parse_single_missense_three_letter():
    [v] = parse_hgvs_pro("p.Thr167Cys")
    assert v.kind == VariantKind.MISSENSE
    assert (v.start, v.end) == (167, 167)
    assert (v.ref_aa, v.alt_aa) == ("T", "C")


def test_parse_single_missense_one_letter():
    [v] = parse_hgvs_pro("p.T167C")
    assert v.kind == VariantKind.MISSENSE
    assert (v.ref_aa, v.alt_aa) == ("T", "C")


def test_parse_synonymous_equals_notation():
    [v] = parse_hgvs_pro("p.Arg200=")
    assert v.kind == VariantKind.SYNONYMOUS


def test_parse_synonymous_self_substitution():
    # p.R200R is technically valid HGVS-ish shorthand; we collapse to synonymous
    [v] = parse_hgvs_pro("p.Arg200Arg")
    assert v.kind == VariantKind.SYNONYMOUS


def test_parse_nonsense_Ter_and_star():
    [v1] = parse_hgvs_pro("p.Arg100Ter")
    [v2] = parse_hgvs_pro("p.Arg100*")
    for v in (v1, v2):
        assert v.kind == VariantKind.NONSENSE
        assert v.ref_aa == "R"
        assert v.alt_aa == "*"


def test_parse_del_single_and_range():
    [vs] = parse_hgvs_pro("p.Arg100del")
    [vr] = parse_hgvs_pro("p.Arg100_Lys110del")
    assert vs.kind == VariantKind.DEL and (vs.start, vs.end) == (100, 100)
    assert vr.kind == VariantKind.DEL and (vr.start, vr.end) == (100, 110)


def test_parse_delins_single_and_range():
    [vs] = parse_hgvs_pro("p.Arg100delinsGlyLeu")
    [vr] = parse_hgvs_pro("p.Arg100_Lys110delinsGL")
    assert vs.kind == VariantKind.DELINS and vs.alt_seq == "GL"
    assert vr.kind == VariantKind.DELINS and vr.alt_seq == "GL"


def test_parse_frameshift_with_offset():
    [v] = parse_hgvs_pro("p.Arg100GlyfsTer10")
    assert v.kind == VariantKind.FRAMESHIFT
    assert v.fs_first_alt_aa == "G"
    assert v.fs_ter_offset == 10


def test_parse_frameshift_minimal():
    [v] = parse_hgvs_pro("p.Arg100fs")
    assert v.kind == VariantKind.FRAMESHIFT
    assert v.fs_first_alt_aa is None


def test_parse_multi_variant_bracketed():
    vs = parse_hgvs_pro("p.[Thr10Gly;Val38Tyr]")
    assert len(vs) == 2
    assert vs[0].kind == VariantKind.MISSENSE and vs[0].start == 10
    assert vs[1].kind == VariantKind.MISSENSE and vs[1].start == 38


def test_parse_multi_variant_mixed_kinds():
    vs = parse_hgvs_pro("p.[Arg5Gly;Arg10Ter;Lys20_Leu30del]")
    kinds = [v.kind for v in vs]
    assert kinds == [VariantKind.MISSENSE, VariantKind.NONSENSE, VariantKind.DEL]


def test_parse_na_returns_empty():
    assert parse_hgvs_pro("NA") == []
    assert parse_hgvs_pro("") == []


def test_parse_malformed_returns_empty():
    assert parse_hgvs_pro("not-an-hgvs") == []
    assert parse_hgvs_pro("p.Xxx123Yyy") == []


def test_parse_parens_wrapper_stripped():
    # HGVS sometimes uses p.(X) for predicted-only notation
    vs = parse_hgvs_pro("p.(Thr167Cys)")
    assert len(vs) == 1


# ─────────────────────────────────────────────────────────────────────────────
# apply_variants — every supported shape
# ─────────────────────────────────────────────────────────────────────────────
def test_apply_single_missense():
    hp = apply_variants(REF, parse_hgvs_pro("p.Arg12Gly"))
    assert hp.seq[11] == "G"
    assert len(hp.seq) == len(REF)
    assert hp.overall_category == "missense_single"
    assert not hp.length_changed


def test_apply_double_missense_both_on_haplotype():
    hp = apply_variants(REF, parse_hgvs_pro("p.[Arg12Gly;Cys13Trp]"))
    assert hp.seq[11] == "G"
    assert hp.seq[12] == "W"
    assert hp.overall_category == "missense_multi"
    assert len(hp.seq) == len(REF)


def test_apply_synonymous_is_noop():
    hp = apply_variants(REF, parse_hgvs_pro("p.Arg12="))
    assert hp.seq == REF
    assert hp.overall_category == "synonymous"


def test_apply_nonsense_truncates():
    hp = apply_variants(REF, parse_hgvs_pro("p.Arg12Ter"))
    assert len(hp.seq) == 11          # everything at and after position 12 dropped
    assert hp.premature_stop_at == 12
    assert hp.overall_category == "lof"


def test_apply_single_residue_del_shortens_by_one():
    hp = apply_variants(REF, parse_hgvs_pro("p.Arg12del"))
    assert len(hp.seq) == len(REF) - 1
    assert hp.overall_category == "indel"
    assert hp.length_changed


def test_apply_range_del_removes_span():
    hp = apply_variants(REF, parse_hgvs_pro("p.Arg12_Cys17del"))
    assert len(hp.seq) == len(REF) - 6      # 12..17 inclusive = 6 residues
    assert hp.overall_category == "indel"


def test_apply_nterminal_del_simulates_whole_gene_loss():
    """Large N-terminal deletion, analogous to an SV that removes the start."""
    hp = apply_variants(REF, parse_hgvs_pro("p.Met1_Thr25del"))
    assert len(hp.seq) == len(REF) - 25
    assert hp.overall_category == "indel"


def test_apply_frameshift_truncates_and_flags_ptc():
    hp = apply_variants(REF, parse_hgvs_pro("p.Cys13GlyfsTer5"))
    # Anchor pos 13 truncated; the frameshift tail is 5 placeholder residues.
    assert hp.overall_category == "lof"
    assert hp.premature_stop_at == 18


def test_apply_ref_mismatch_is_skipped_not_errored():
    """If the reference amino acid in the HGVS doesn't match, we skip rather
    than corrupt the haplotype. Downstream consumers can look at `skipped`
    to decide whether to flag the row."""
    hp = apply_variants(REF, parse_hgvs_pro("p.Trp12Gly"))   # pos 12 is R not W
    assert len(hp.applied) == 0
    assert len(hp.skipped) == 1
    assert any("ref_mismatch" in c for c in hp.categories)


def test_apply_out_of_range_is_skipped():
    hp = apply_variants(REF, parse_hgvs_pro("p.Arg999Gly"))
    assert len(hp.applied) == 0
    assert any("out_of_range" in c for c in hp.categories)


def test_apply_compound_multi_variant_preserves_both():
    """Variants are applied in reverse-coordinate order so an upstream
    deletion doesn't shift the downstream missense position."""
    hp = apply_variants(REF, parse_hgvs_pro("p.[Ala2Trp;Arg12Gly]"))
    assert hp.seq[1] == "W"                                  # A2W applied
    assert hp.seq[11] == "G"                                 # R12G applied
    assert hp.overall_category == "missense_multi"


def test_apply_delins_replaces_with_specified_sequence():
    hp = apply_variants(REF, parse_hgvs_pro("p.Arg12delinsValLeu"))
    # Range [12,12] replaced with 'VL' → length increases by 1
    assert hp.seq[11] == "V"
    assert hp.seq[12] == "L"
    assert len(hp.seq) == len(REF) + 1


def test_apply_empty_variant_list_returns_unchanged():
    hp = apply_variants(REF, [])
    assert hp.seq == REF
    assert hp.overall_category == "unmodified"
    assert len(hp.applied) == 0
