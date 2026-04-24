"""Unit tests for the Phase 4 diploid aggregator + constraint loader."""

from __future__ import annotations

from pathlib import Path

import pytest

from hapli.core.schema import (
    ConstraintPriors,
    GeneEvidence,
    PresenceCall,
    ProteinDiff,
)
from hapli.external.constraint import (
    ConstraintLookup,
    load_clingen_dosage,
    load_gnomad_constraint,
)
from hapli.interpretation.diploid import (
    LOF_THRESHOLD,
    build_diploid_report,
    is_lof,
    score_haplotype,
)


# ─────────────────────────────────────────────────────────────────────────────
# score_haplotype
# ─────────────────────────────────────────────────────────────────────────────
def test_deleted_presence_scores_zero():
    p = PresenceCall(status="deleted", source="liftoff")
    assert score_haplotype(p, []) == 0.0


def test_not_run_scores_none():
    p = PresenceCall(status="not_run", source="liftoff")
    assert score_haplotype(p) is None


def test_intact_with_perfect_protein_scores_one():
    p = PresenceCall(status="intact", source="liftoff")
    diffs = [ProteinDiff(transcript="T1", haplotype=1, ref_length=40, hap_length=40, identity=1.0)]
    assert score_haplotype(p, diffs) == 1.0


def test_intact_with_partial_identity_scores_proportionally():
    p = PresenceCall(status="intact", source="liftoff")
    diffs = [ProteinDiff(transcript="T1", haplotype=1, ref_length=100, hap_length=100, identity=0.85)]
    assert score_haplotype(p, diffs) == 0.85


def test_nterm_premature_stop_scores_zero():
    """PTC in first 10% of the protein → effectively no functional product."""
    p = PresenceCall(status="intact", source="liftoff")
    diffs = [ProteinDiff(
        transcript="T1", haplotype=1,
        ref_length=100, hap_length=10, identity=0.10,
        premature_stop_at=5,
    )]
    assert score_haplotype(p, diffs) == 0.0


def test_mid_protein_premature_stop_scales_with_position():
    p = PresenceCall(status="intact", source="liftoff")
    diffs = [ProteinDiff(
        transcript="T1", haplotype=1,
        ref_length=100, hap_length=50, identity=0.50,
        premature_stop_at=50,
    )]
    # identity=0.5, ptc at position 50 out of 100 → 0.5 * 0.5 = 0.25
    s = score_haplotype(p, diffs)
    assert s is not None
    assert abs(s - 0.25) < 1e-9


def test_multi_transcript_takes_min():
    """Conservative: gene is only as functional as its worst isoform."""
    p = PresenceCall(status="intact", source="liftoff")
    diffs = [
        ProteinDiff(transcript="T1", haplotype=1, ref_length=100, hap_length=100, identity=0.99),
        ProteinDiff(transcript="T2", haplotype=1, ref_length=100, hap_length=100, identity=0.50),
    ]
    assert score_haplotype(p, diffs) == 0.50


def test_low_identity_presence_without_diffs_falls_back_to_sequence_id():
    p = PresenceCall(
        status="low_identity", source="liftoff", sequence_identity=0.73,
    )
    assert score_haplotype(p, []) == 0.73


def test_low_identity_presence_caps_protein_only_score():
    """Inversion regression: Liftoff reports the gene as low_identity (the
    inverted region couldn't be mapped cleanly) but the protein-extraction
    step translated only the unaffected flank and returned identity=1.0.
    The haplotype score must be dampened by presence.sequence_identity —
    otherwise a structurally disrupted gene gets a misleading score of 1.0.
    """
    p = PresenceCall(
        status="low_identity", source="liftoff", sequence_identity=0.40,
    )
    d = ProteinDiff(
        transcript="T1", haplotype=1, ref_length=100, hap_length=60,
        identity=1.0, premature_stop_at=None, frameshift_region=None,
        per_position=[],
    )
    # Protein looks perfect on its subset, but Liftoff's identity caps it:
    assert score_haplotype(p, [d]) == 0.40


def test_partial_presence_also_caps_protein_only_score():
    p = PresenceCall(
        status="partial", source="liftoff", sequence_identity=0.55,
    )
    d = ProteinDiff(
        transcript="T1", haplotype=1, ref_length=100, hap_length=100,
        identity=1.0, premature_stop_at=None, frameshift_region=None,
        per_position=[],
    )
    assert score_haplotype(p, [d]) == 0.55


def test_intact_presence_does_not_cap_protein_score():
    """The cap only fires on low_identity / partial; an 'intact' presence
    with sequence_identity=0.99 must NOT drag a protein identity of 1.0 down.
    """
    p = PresenceCall(
        status="intact", source="liftoff", sequence_identity=0.99,
    )
    d = ProteinDiff(
        transcript="T1", haplotype=1, ref_length=100, hap_length=100,
        identity=1.0, premature_stop_at=None, frameshift_region=None,
        per_position=[],
    )
    assert score_haplotype(p, [d]) == 1.0


# ─────────────────────────────────────────────────────────────────────────────
# is_lof
# ─────────────────────────────────────────────────────────────────────────────
def test_is_lof_threshold_behaviour():
    assert is_lof(0.0)
    assert is_lof(LOF_THRESHOLD)
    assert not is_lof(LOF_THRESHOLD + 0.01)
    assert not is_lof(None)
    # Custom threshold
    assert is_lof(0.3, threshold=0.4)
    assert not is_lof(0.5, threshold=0.4)


# ─────────────────────────────────────────────────────────────────────────────
# build_diploid_report
# ─────────────────────────────────────────────────────────────────────────────
def _evidence_with(
    hap1_presence: str,
    hap2_presence: str,
    hap1_identity: float | None = None,
    hap2_identity: float | None = None,
) -> GeneEvidence:
    ev = GeneEvidence(gene="G1")
    ev.presence["hap1"] = PresenceCall(status=hap1_presence, source="liftoff")
    ev.presence["hap2"] = PresenceCall(status=hap2_presence, source="liftoff")
    if hap1_identity is not None:
        ev.protein.append(ProteinDiff(
            transcript="T1", haplotype=1,
            ref_length=100, hap_length=100, identity=hap1_identity,
        ))
    if hap2_identity is not None:
        ev.protein.append(ProteinDiff(
            transcript="T1", haplotype=2,
            ref_length=100, hap_length=100, identity=hap2_identity,
        ))
    return ev


def test_compound_het_lof_both_haps_below_threshold():
    """Classic compound-het LoF: both haplotypes independently non-functional."""
    ev = _evidence_with("deleted", "deleted")
    rep = build_diploid_report(ev)
    assert rep.hap1_score == 0.0
    assert rep.hap2_score == 0.0
    assert rep.compound_het_lof is True


def test_not_compound_het_when_one_hap_functional():
    """Case 08-like: one deleted, one intact → not compound het LoF."""
    ev = _evidence_with("deleted", "intact", None, 1.0)
    rep = build_diploid_report(ev)
    assert rep.hap1_score == 0.0
    assert rep.hap2_score == 1.0
    assert rep.compound_het_lof is False


def test_rescue_case_is_not_compound_het_lof():
    """Case 05-like: hap1 has rescued frameshift → still scores high; hap2 reference."""
    ev = _evidence_with("intact", "intact", 0.85, 1.0)
    rep = build_diploid_report(ev)
    assert rep.hap1_score == 0.85
    assert rep.hap2_score == 1.0
    assert rep.compound_het_lof is False


def test_min_max_skip_none_values():
    ev = _evidence_with("not_run", "intact", None, 0.9)
    rep = build_diploid_report(ev)
    assert rep.hap1_score is None
    assert rep.hap2_score == 0.9
    assert rep.min_score == 0.9
    assert rep.max_score == 0.9
    # Can't call compound_het_lof when one hap is unknown
    assert rep.compound_het_lof is False


def test_constraints_propagate_verbatim():
    ev = _evidence_with("intact", "intact", 1.0, 1.0)
    priors = ConstraintPriors(pli=0.95, mis_z=2.1, clingen_haploinsufficient=True)
    rep = build_diploid_report(ev, constraints=priors)
    assert rep.constraints.pli == 0.95
    assert rep.constraints.mis_z == 2.1
    assert rep.constraints.clingen_haploinsufficient is True


# ─────────────────────────────────────────────────────────────────────────────
# constraint loader
# ─────────────────────────────────────────────────────────────────────────────
def test_gnomad_loader_picks_canonical_transcript(tmp_path):
    f = tmp_path / "gnomad.tsv"
    f.write_text(
        "gene\ttranscript\tcanonical\tlof.pLI\tmis.z_score\tlof.oe\tmis.oe\n"
        "BRCA1\tENST0001\ttrue\t0.998\t2.1\t0.12\t0.85\n"
        "BRCA1\tENST0002\tfalse\t0.50\t1.0\t0.5\t0.5\n"
    )
    out = load_gnomad_constraint(f)
    assert out["BRCA1"]["pli"] == 0.998
    assert out["BRCA1"]["mis_z"] == 2.1


def test_gnomad_loader_tolerates_missing_canonical_column(tmp_path):
    f = tmp_path / "gnomad.tsv"
    f.write_text(
        "gene\tpLI\tmis_z\n"
        "FOO\t0.1\t0.5\n"
    )
    out = load_gnomad_constraint(f)
    assert out["FOO"]["pli"] == 0.1
    assert out["FOO"]["oe_lof"] is None


def test_clingen_loader_interprets_scores(tmp_path):
    f = tmp_path / "clingen.tsv"
    f.write_text(
        "Gene Symbol\tHaploinsufficiency Score\tTriplosensitivity Score\n"
        "HI3\t3\t0\n"
        "TS3\t0\t3\n"
        "AR30\t30\t0\n"          # autosomal recessive, not HI
        "UNLIKELY40\t40\t40\n"   # dosage sensitivity unlikely
        "NODATA\tNA\tNA\n"
    )
    out = load_clingen_dosage(f)
    assert out["HI3"]["haploinsufficient"] is True
    assert out["HI3"]["triplosensitive"] is False
    assert out["TS3"]["haploinsufficient"] is False
    assert out["TS3"]["triplosensitive"] is True
    assert out["AR30"]["haploinsufficient"] is False
    assert out["UNLIKELY40"]["haploinsufficient"] is False
    assert out["NODATA"]["haploinsufficient"] is None


def test_constraint_lookup_unified_api(tmp_path):
    g = tmp_path / "g.tsv"
    g.write_text("gene\tpLI\n" "BRCA1\t0.998\n")
    c = tmp_path / "c.tsv"
    c.write_text("Gene Symbol\tHaploinsufficiency Score\tTriplosensitivity Score\n"
                 "BRCA1\t3\t0\n")
    lookup = ConstraintLookup.load(gnomad_path=g, clingen_path=c)
    priors = lookup.get("BRCA1")
    assert priors.pli == 0.998
    assert priors.clingen_haploinsufficient is True
    # Unknown gene → all None.
    missing = lookup.get("UNKNOWN")
    assert missing.pli is None
    assert missing.clingen_haploinsufficient is None


def test_constraint_lookup_with_no_sources_returns_empty_priors():
    lookup = ConstraintLookup()
    p = lookup.get("ANYTHING")
    assert p.pli is None and p.clingen_haploinsufficient is None
