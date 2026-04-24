"""Tests for hapli.core.protein_diff."""

from __future__ import annotations

from hapli.core.protein_diff import diff_proteins


REF = "MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"   # 40 aa


def test_identical_sequences_report_full_identity():
    d = diff_proteins("T1", 1, REF, REF)
    assert d.identity == 1.0
    assert d.substitutions == []
    assert d.premature_stop_at is None
    assert d.frameshift_region is None


def test_single_missense_is_tracked_without_frameshift_flag():
    hap = "M" + "A" * 10 + "V" + "A" * 28                # single A→V at aa 12
    d = diff_proteins("T1", 1, REF, hap)
    assert d.identity == 39 / 40
    assert len(d.substitutions) == 1
    assert d.substitutions[0]["ref_pos"] == 12
    assert d.substitutions[0]["ref_aa"] == "A"
    assert d.substitutions[0]["hap_aa"] == "V"
    assert d.frameshift_region is None


def test_frameshift_rescue_window_is_flagged_as_restored():
    """The +1 / -1 hero case. Same length, dense substitution cluster."""
    hap = "MAAARCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    d = diff_proteins("T1", 1, REF, hap)
    assert d.ref_length == d.hap_length == 40
    assert d.frameshift_region is not None
    assert d.frameshift_region["restored"] is True
    assert d.frameshift_region["start"] == 5
    assert d.frameshift_region["end"] == 10
    assert d.frameshift_region["n_substitutions"] >= 3


def test_premature_stop_in_truncated_hap_is_flagged():
    hap = "MAAARCCCCC*"                                 # 10 aa + stop, much shorter than ref
    d = diff_proteins("T1", 1, REF, hap)
    assert d.premature_stop_at == 11


def test_natural_terminal_stop_on_full_length_hap_is_not_premature():
    hap = REF + "*"
    d = diff_proteins("T1", 1, REF, hap)
    assert d.premature_stop_at is None


def test_mid_sequence_stop_is_premature_even_if_sequence_continues():
    hap = "MAAA*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"     # stop at position 5 but residues continue
    d = diff_proteins("T1", 1, REF, hap)
    assert d.premature_stop_at == 5


def test_empty_haplotype_identity_is_zero():
    d = diff_proteins("T1", 1, REF, "")
    assert d.identity == 0.0
    assert d.hap_length == 0


def test_both_empty_is_trivially_identical():
    d = diff_proteins("T1", 1, "", "")
    assert d.identity == 1.0
