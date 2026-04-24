"""Tests for the population-level aggregator."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from hapli.workflow.aggregate import (
    GeneAggregate,
    SampleGeneRow,
    aggregate_files,
    aggregate_glob,
    parse_filename,
    write_per_gene_tsv,
    write_per_sample_tsv,
)


def _write_json(path: Path, data: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return path


def _payload(
    sample: str,
    gene: str,
    *,
    h1: float | None = 1.0,
    h2: float | None = 1.0,
    chl: bool = False,
    pres1: str = "intact",
    pres2: str = "intact",
    consequences: list[dict] | None = None,
    epistasis: list[dict] | None = None,
    pli: float | None = None,
    clingen: bool | None = None,
) -> dict:
    return {
        "schema_version": "2.0",
        "gene": gene,
        "sample": sample,
        "region": "chr1:1-100",
        "transcripts": [],
        "evidence": {
            "presence": {
                "hap1": {"status": pres1},
                "hap2": {"status": pres2},
            },
            "consequence": consequences or [],
            "epistasis": epistasis or [],
            "diploid": {
                "hap1_score": h1,
                "hap2_score": h2,
                "compound_het_lof": chl,
                "constraints": {
                    "pli": pli,
                    "clingen_haploinsufficient": clingen,
                },
            },
        },
    }


def test_parse_filename_extracts_sample_and_gene(tmp_path):
    p = tmp_path / "S1_BRCA1_analysis.json"
    sample, gene = parse_filename(p)
    assert sample == "S1"
    assert gene == "BRCA1"


def test_parse_filename_handles_underscored_genes(tmp_path):
    """Greedy split — 'S1_PAX_6_analysis.json' splits as sample=S1, gene=PAX_6."""
    p = tmp_path / "S1_PAX_6_analysis.json"
    sample, gene = parse_filename(p)
    assert sample == "S1"
    assert gene == "PAX_6"


def test_aggregate_basic(tmp_path):
    f1 = _write_json(tmp_path / "HG002_BRCA1_analysis.json", _payload(
        "HG002", "BRCA1", h1=0.0, h2=0.0, chl=True,
        pres1="deleted", pres2="intact",
        consequences=[
            {"haplotype": 1, "consequence": "stop_gained"},
            {"haplotype": 2, "consequence": "frameshift"},
        ],
    ))
    f2 = _write_json(tmp_path / "HG003_BRCA1_analysis.json", _payload(
        "HG003", "BRCA1", h1=1.0, h2=1.0,
    ))
    f3 = _write_json(tmp_path / "HG002_TP53_analysis.json", _payload(
        "HG002", "TP53", h1=0.7, h2=1.0,
        epistasis=[{"residual": -8.5, "flagged": True}],
    ))
    rows, aggs = aggregate_files([f1, f2, f3])
    by_key = {(r.sample, r.gene): r for r in rows}
    assert by_key[("HG002", "BRCA1")].compound_het_lof is True
    assert by_key[("HG002", "BRCA1")].n_lof_consequences == 2
    assert by_key[("HG002", "BRCA1")].presence_hap1 == "deleted"
    assert by_key[("HG003", "BRCA1")].compound_het_lof is False
    assert by_key[("HG002", "TP53")].epistasis_flagged is True
    assert by_key[("HG002", "TP53")].max_residual_abs == 8.5

    # Per-gene aggregates
    by_gene = {a.gene: a for a in aggs}
    brca = by_gene["BRCA1"]
    assert brca.n_samples == 2
    assert brca.n_compound_het_lof_samples == 1
    assert brca.n_lof_alleles == 2          # both haps of HG002 below 0.5
    assert brca.n_hap_records == 4          # 2 samples × 2 haps each
    assert brca.n_deleted_alleles == 1      # HG002 hap1=deleted
    tp53 = by_gene["TP53"]
    assert tp53.n_samples == 1
    assert tp53.n_epistasis_flagged_samples == 1


def test_aggregate_handles_missing_evidence(tmp_path):
    """A v1-shaped file (no `evidence` block) shouldn't crash the aggregator."""
    f = _write_json(tmp_path / "S1_FOO_analysis.json", {
        "gene": "FOO", "sample": "S1", "region": "chr1:1-100", "transcripts": [],
    })
    rows, aggs = aggregate_files([f])
    assert len(rows) == 1
    r = rows[0]
    assert r.hap1_score is None and r.hap2_score is None
    assert r.compound_het_lof is False
    assert r.n_consequence == 0


def test_glob_pattern(tmp_path):
    _write_json(tmp_path / "S1_G1_analysis.json", _payload("S1", "G1"))
    _write_json(tmp_path / "S2_G1_analysis.json", _payload("S2", "G1", h1=0.0, h2=0.0, chl=True))
    _write_json(tmp_path / "S1_G2_analysis.json", _payload("S1", "G2"))
    (tmp_path / "ignore_me.txt").write_text("not json")
    rows, aggs = aggregate_glob(str(tmp_path / "*_analysis.json"))
    assert len(rows) == 3
    by_gene = {a.gene: a for a in aggs}
    assert by_gene["G1"].n_compound_het_lof_samples == 1
    assert by_gene["G2"].n_samples == 1


def test_tsv_writers_round_trip(tmp_path):
    f = _write_json(tmp_path / "S1_G1_analysis.json", _payload(
        "S1", "G1", h1=0.5, h2=0.9, chl=False, pli=0.95, clingen=True,
    ))
    rows, aggs = aggregate_files([f])
    out_s = tmp_path / "per_sample.tsv"
    out_g = tmp_path / "per_gene.tsv"
    write_per_sample_tsv(rows, out_s)
    write_per_gene_tsv(aggs, out_g)
    s_text = out_s.read_text().splitlines()
    g_text = out_g.read_text().splitlines()
    assert s_text[0].startswith("sample\tgene")
    assert "S1\tG1" in s_text[1]
    assert g_text[0].startswith("gene\tn_samples")
    assert "G1\t1" in g_text[1]
    assert "0.95" in g_text[1]      # pli forwarded
