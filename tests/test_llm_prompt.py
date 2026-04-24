"""
Unit tests for the LLM prompt builder.

The network call is out of scope (OpenRouter key required) but the prompt
construction must pull the schema-v2 evidence fields so the clinician-facing
report includes presence calls, compound-het LoF flags, epistasis residuals,
and ClinVar annotations — not just the schema-v1 minimap2 alignment block.
"""
from __future__ import annotations

from hapli.interpretation.llm import LLMInterpreter


def _render(data):
    itp = LLMInterpreter(api_key=None)
    return itp._construct_prompt(data)


def test_prompt_includes_compound_het_lof_flag():
    prompt = _render({
        "gene": "BRCA1", "sample": "S1",
        "evidence": {
            "presence": {
                "hap1": {"status": "uncertain", "sequence_identity": 0.95},
                "hap2": {"status": "uncertain", "sequence_identity": 0.98},
            },
            "diploid": {
                "hap1_score": 0.0, "hap2_score": 0.3,
                "compound_het_lof": True,
            },
        },
    })
    assert "compound_het_lof: TRUE" in prompt
    assert "hap1_score: 0.000" in prompt
    assert "hap2_score: 0.300" in prompt
    assert "BRCA1" in prompt


def test_prompt_includes_epistasis_residual():
    prompt = _render({
        "gene": "G1", "sample": "S1",
        "evidence": {
            "epistasis": [
                {"haplotype": 1, "residual": 57.05,
                 "s_additive": -58.48, "s_joint": -1.43,
                 "flagged": True},
            ],
        },
    })
    assert "residual +57.05" in prompt
    assert "[FLAGGED]" in prompt
    assert "S_add=-58.48" in prompt


def test_prompt_surfaces_clnsig_annotations():
    prompt = _render({
        "gene": "G1", "sample": "S1",
        "evidence": {
            "consequence": [
                {"haplotype": 1, "consequence": "missense",
                 "amino_acid_change": "p.Ala10Ser", "clnsig": "pathogenic"},
            ],
        },
    })
    assert "p.Ala10Ser" in prompt
    assert "ClinVar: pathogenic" in prompt


def test_prompt_reports_presence_deleted():
    prompt = _render({
        "gene": "G1", "sample": "S1",
        "evidence": {
            "presence": {
                "hap1": {"status": "deleted"},
                "hap2": {"status": "intact", "sequence_identity": 1.0},
            },
        },
    })
    assert "hap1: deleted" in prompt
    assert "hap2: intact" in prompt


def test_prompt_falls_back_to_legacy_schema_when_no_evidence():
    """Backward compatibility: older analysis JSONs with only the schema-v1
    `transcripts[].alignments` block must still produce a usable prompt."""
    prompt = _render({
        "gene": "G1", "sample": "S1",
        "transcripts": [
            {"id": "T1", "alignments": {
                "hap1": {"identity": 0.98, "nm": 5, "cigar": "100M"},
                "hap2": {"identity": 1.0, "nm": 0, "cigar": "100M"},
            }},
        ],
    })
    assert "Transcript T1" in prompt
    assert "hap1: Identity 98" in prompt


def test_prompt_reports_frameshift_restored():
    prompt = _render({
        "gene": "G1", "sample": "S1",
        "evidence": {
            "protein": [
                {"haplotype": 1, "identity": 0.85,
                 "premature_stop_at": None,
                 "frameshift_region": {"restored": True}},
                {"haplotype": 2, "identity": 1.0,
                 "premature_stop_at": None,
                 "frameshift_region": None},
            ],
        },
    })
    assert "frameshift restored" in prompt
