"""
Population-level aggregation across many `hapli analyze` / `assess` outputs.

For an HPRC-scale run (200 samples × N genes) the natural downstream
consumer is "what's the population LoF allele frequency per gene, and
where are the compound-het-LoF / epistasis hotspots?". This module walks
a glob of analysis JSONs, projects each to a comparable per-(gene, sample)
row, and emits:

  * one TSV row per (gene, sample) with key signals from evidence.diploid
    + evidence.epistasis,
  * one TSV row per gene with population aggregates (LoF allele frequency,
    compound-het-LoF count, mean residual flag rate, etc.).

Designed to be cheap (one streaming pass), robust to missing fields
(uses `read_evidence` defensive lookups), and forward-compatible with the
schema-v2 `evidence.*` block.
"""

from __future__ import annotations

import csv
import glob
import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from ..core.schema import read_evidence


# Standard hapli output filename: <sample>_<gene>_analysis.json
# (gene names can contain '_' so the split is rightmost on '_analysis.json')
_FNAME_RE = re.compile(r"^(?P<sample>[^/]+?)_(?P<gene>[^/]+)_analysis\.json$")


@dataclass
class SampleGeneRow:
    """One row per (sample, gene) in the per-sample TSV."""

    sample: str
    gene: str
    file: str
    hap1_score: float | None = None
    hap2_score: float | None = None
    min_score: float | None = None
    compound_het_lof: bool = False
    presence_hap1: str | None = None
    presence_hap2: str | None = None
    n_consequence: int = 0
    n_lof_consequences: int = 0
    n_epistasis_records: int = 0
    epistasis_flagged: bool = False
    max_residual_abs: float | None = None
    pli: float | None = None
    clingen_haploinsufficient: bool | None = None


@dataclass
class GeneAggregate:
    """Per-gene population summary collapsed across all samples."""

    gene: str
    n_samples: int = 0
    n_compound_het_lof_samples: int = 0
    n_epistasis_flagged_samples: int = 0
    # n_lof_alleles / 2*n_samples = LoF allele frequency
    n_lof_alleles: int = 0
    # n_hap_records counts how many haplotypes we actually scored — paired with
    # n_lof_alleles to give an allele frequency that ignores not_run haps.
    n_hap_records: int = 0
    n_deleted_alleles: int = 0
    pli: float | None = None
    clingen_haploinsufficient: bool | None = None
    samples: list[str] = field(default_factory=list)


_LOF_TAGS = {
    "stop_gained", "stop_lost", "start_lost", "frameshift",
    "splice_acceptor", "splice_donor",
}


def _is_lof_score(score: float | None, threshold: float = 0.5) -> bool:
    return score is not None and score <= threshold


def parse_filename(path: Path) -> tuple[str | None, str | None]:
    m = _FNAME_RE.match(path.name)
    if not m:
        return None, None
    return m["sample"], m["gene"]


def aggregate_files(
    paths: Iterable[Path],
    *,
    lof_threshold: float = 0.5,
    logger: logging.Logger | None = None,
) -> tuple[list[SampleGeneRow], list[GeneAggregate]]:
    """Stream-read all `*_analysis.json` paths and produce both tables."""
    logger = logger or logging.getLogger(__name__)
    rows: list[SampleGeneRow] = []
    by_gene: dict[str, GeneAggregate] = {}

    for path in paths:
        try:
            data = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning("could not read %s: %s", path, exc)
            continue
        # Prefer the JSON's own sample/gene fields (always present in schema-v2);
        # the filename split is ambiguous when sample names contain '_'
        # (e.g. "sample_001_GENE00001_analysis.json" splits as
        # sample='sample', gene='001_GENE00001' with a generic regex).
        sample = read_evidence(data, "sample")
        gene = read_evidence(data, "gene")
        if not sample or not gene:
            fallback_sample, fallback_gene = parse_filename(path)
            sample = sample or fallback_sample or path.stem
            gene = gene or fallback_gene or "?"

        # ── Per-sample row ──
        h1 = read_evidence(data, "evidence", "diploid", "hap1_score")
        h2 = read_evidence(data, "evidence", "diploid", "hap2_score")
        chl = bool(read_evidence(data, "evidence", "diploid", "compound_het_lof", default=False))
        pres1 = read_evidence(data, "evidence", "presence", "hap1", "status")
        pres2 = read_evidence(data, "evidence", "presence", "hap2", "status")
        consequences = data.get("evidence", {}).get("consequence", []) or []
        n_lof_csq = sum(1 for c in consequences if c.get("consequence") in _LOF_TAGS)
        epi = data.get("evidence", {}).get("epistasis", []) or []
        epi_flagged = any(bool(e.get("flagged")) for e in epi)
        max_abs_res = None
        if epi:
            try:
                max_abs_res = max(abs(float(e.get("residual", 0))) for e in epi)
            except (TypeError, ValueError):
                max_abs_res = None
        pli = read_evidence(data, "evidence", "diploid", "constraints", "pli")
        clingen_hi = read_evidence(data, "evidence", "diploid", "constraints", "clingen_haploinsufficient")

        min_score = None
        scored = [s for s in (h1, h2) if s is not None]
        if scored:
            min_score = min(scored)

        row = SampleGeneRow(
            sample=str(sample), gene=str(gene), file=str(path),
            hap1_score=h1, hap2_score=h2, min_score=min_score,
            compound_het_lof=chl,
            presence_hap1=pres1, presence_hap2=pres2,
            n_consequence=len(consequences), n_lof_consequences=n_lof_csq,
            n_epistasis_records=len(epi),
            epistasis_flagged=epi_flagged,
            max_residual_abs=max_abs_res,
            pli=pli,
            clingen_haploinsufficient=clingen_hi,
        )
        rows.append(row)

        # ── Per-gene aggregate accumulator ──
        agg = by_gene.setdefault(gene, GeneAggregate(gene=gene, pli=pli,
                                                    clingen_haploinsufficient=clingen_hi))
        agg.n_samples += 1
        agg.samples.append(str(sample))
        if chl:
            agg.n_compound_het_lof_samples += 1
        if epi_flagged:
            agg.n_epistasis_flagged_samples += 1
        for s in (h1, h2):
            if s is None:
                continue
            agg.n_hap_records += 1
            if _is_lof_score(s, lof_threshold):
                agg.n_lof_alleles += 1
        for p in (pres1, pres2):
            if p == "deleted":
                agg.n_deleted_alleles += 1
        # Forward priors from any sample that has them; they should be identical
        # across samples for the same gene.
        if agg.pli is None and pli is not None:
            agg.pli = pli
        if agg.clingen_haploinsufficient is None and clingen_hi is not None:
            agg.clingen_haploinsufficient = clingen_hi

    return rows, list(by_gene.values())


def write_per_sample_tsv(rows: list[SampleGeneRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "sample", "gene",
            "hap1_score", "hap2_score", "min_score",
            "compound_het_lof",
            "presence_hap1", "presence_hap2",
            "n_consequence", "n_lof_consequences",
            "n_epistasis_records", "epistasis_flagged", "max_residual_abs",
            "pli", "clingen_haploinsufficient",
            "file",
        ])
        for r in rows:
            w.writerow([
                r.sample, r.gene,
                _f(r.hap1_score), _f(r.hap2_score), _f(r.min_score),
                int(r.compound_het_lof),
                r.presence_hap1 or "", r.presence_hap2 or "",
                r.n_consequence, r.n_lof_consequences,
                r.n_epistasis_records, int(r.epistasis_flagged),
                _f(r.max_residual_abs),
                _f(r.pli),
                "" if r.clingen_haploinsufficient is None else int(r.clingen_haploinsufficient),
                r.file,
            ])


def write_per_gene_tsv(aggregates: list[GeneAggregate], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "gene", "n_samples", "n_hap_records", "n_lof_alleles",
            "lof_allele_freq", "n_deleted_alleles",
            "n_compound_het_lof_samples", "compound_het_lof_freq",
            "n_epistasis_flagged_samples", "epistasis_flagged_freq",
            "pli", "clingen_haploinsufficient",
        ])
        for a in sorted(aggregates, key=lambda x: x.gene):
            af = a.n_lof_alleles / a.n_hap_records if a.n_hap_records else 0.0
            chf = a.n_compound_het_lof_samples / a.n_samples if a.n_samples else 0.0
            ef = a.n_epistasis_flagged_samples / a.n_samples if a.n_samples else 0.0
            w.writerow([
                a.gene, a.n_samples, a.n_hap_records, a.n_lof_alleles,
                f"{af:.6f}", a.n_deleted_alleles,
                a.n_compound_het_lof_samples, f"{chf:.6f}",
                a.n_epistasis_flagged_samples, f"{ef:.6f}",
                _f(a.pli),
                "" if a.clingen_haploinsufficient is None else int(a.clingen_haploinsufficient),
            ])


def _f(v: float | None) -> str:
    if v is None:
        return ""
    return f"{v:.6g}"


def aggregate_glob(
    pattern: str,
    *,
    lof_threshold: float = 0.5,
    logger: logging.Logger | None = None,
) -> tuple[list[SampleGeneRow], list[GeneAggregate]]:
    """Resolve a glob pattern and aggregate the matching files."""
    paths = sorted(Path(p) for p in glob.glob(pattern, recursive=True))
    return aggregate_files(paths, lof_threshold=lof_threshold, logger=logger)


__all__ = [
    "GeneAggregate",
    "SampleGeneRow",
    "aggregate_files",
    "aggregate_glob",
    "parse_filename",
    "write_per_gene_tsv",
    "write_per_sample_tsv",
]
