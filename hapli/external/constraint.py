"""
Constraint-prior lookups for per-gene gnomAD + ClinGen metrics.

The production gnomAD v4 constraint table lives at

    https://gnomad.broadinstitute.org/downloads#v4-constraint

and ClinGen dosage sensitivity at

    https://search.clinicalgenome.org/kb/gene-dosage

Both are TSV-shaped. We ingest them into two in-memory dictionaries keyed by
HGNC gene symbol, then expose a single `ConstraintLookup.get(gene)` that
returns a populated `ConstraintPriors`.

The loader is forgiving:
  * Either column source may be absent → lookups still succeed but the
    corresponding fields are None.
  * gnomAD's per-transcript rows are collapsed to the canonical transcript
    (canonical=True) when the `canonical` column is present; otherwise the
    first row per gene wins.
  * Column names that vary across gnomAD releases are normalised via a small
    alias table.
"""

from __future__ import annotations

import csv
import gzip
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, IO

from ..core.schema import ConstraintPriors


# ─────────────────────────────────────────────────────────────────────────────
# gnomAD column name aliases. Canonical name → list of aliases observed across
# v2.1, v3, v4.1 release TSVs. When adding a new column, extend here rather
# than inserting branch-by-branch in the code.
# ─────────────────────────────────────────────────────────────────────────────
_GNOMAD_COLUMN_ALIASES = {
    "gene": ["gene", "gene_symbol", "gene_name"],
    "transcript": ["transcript", "transcript_id"],
    "canonical": ["canonical", "is_canonical", "mane_select"],
    "pli": ["pLI", "lof.pLI", "lof_pLI", "pli"],
    "mis_z": ["mis_z", "mis.z_score", "mis_z_score"],
    "oe_lof": ["oe_lof", "lof.oe", "lof_oe"],
    "oe_mis": ["oe_mis", "mis.oe", "mis_oe"],
}


# ─────────────────────────────────────────────────────────────────────────────
# ClinGen columns. Dosage TSV canonical columns (as of 2024 release):
#   Gene Symbol, Haploinsufficiency Score, Triplosensitivity Score
# Scores are integers:
#   3: sufficient evidence
#   2: some evidence
#   1: little evidence
#   0: no evidence
#   30: autosomal recessive
#   40: dosage sensitivity unlikely
# We treat scores >=2 as positive (evidence of dosage sensitivity).
# ─────────────────────────────────────────────────────────────────────────────
_CLINGEN_COLUMN_ALIASES = {
    "gene": ["Gene Symbol", "gene_symbol", "gene"],
    "hi": ["Haploinsufficiency Score", "haploinsufficiency_score", "hi_score"],
    "ts": ["Triplosensitivity Score", "triplosensitivity_score", "ts_score"],
}


def _open_table(path: Path) -> IO[str]:
    """Open a (possibly gzipped) TSV in text mode."""
    if str(path).endswith(".gz") or str(path).endswith(".bgz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def _resolve_columns(header: list[str], aliases: dict[str, list[str]]) -> dict[str, int]:
    """Return {canonical_name: column_index}. Missing keys are simply absent
    from the result — the loader treats them as optional."""
    lowered = {h.lower(): i for i, h in enumerate(header)}
    direct = {h: i for i, h in enumerate(header)}
    out: dict[str, int] = {}
    for canonical, names in aliases.items():
        for n in names:
            if n in direct:
                out[canonical] = direct[n]; break
            if n.lower() in lowered:
                out[canonical] = lowered[n.lower()]; break
    return out


def _parse_float(raw: str) -> float | None:
    if raw is None or raw == "" or raw.upper() == "NA":
        return None
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


# ─────────────────────────────────────────────────────────────────────────────
# Loaders
# ─────────────────────────────────────────────────────────────────────────────
def load_gnomad_constraint(path: Path | str) -> dict[str, dict[str, Any]]:
    """Return {gene_symbol: {pli, mis_z, oe_lof, oe_mis}}.

    When multiple transcripts are present per gene we prefer the canonical one.
    """
    path = Path(path)
    out: dict[str, dict[str, Any]] = {}
    with _open_table(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        cols = _resolve_columns(header, _GNOMAD_COLUMN_ALIASES)
        if "gene" not in cols:
            raise ValueError(f"{path}: no gene column (have: {header})")
        for row in reader:
            if not row:
                continue
            gene = row[cols["gene"]].strip()
            if not gene:
                continue
            # Prefer the canonical transcript when the column exists.
            is_canonical = True
            if "canonical" in cols:
                val = row[cols["canonical"]].strip().lower()
                is_canonical = val in {"true", "yes", "1", "canonical"}
            if gene in out and not is_canonical:
                continue
            out[gene] = {
                "pli": _parse_float(row[cols["pli"]]) if "pli" in cols else None,
                "mis_z": _parse_float(row[cols["mis_z"]]) if "mis_z" in cols else None,
                "oe_lof": _parse_float(row[cols["oe_lof"]]) if "oe_lof" in cols else None,
                "oe_mis": _parse_float(row[cols["oe_mis"]]) if "oe_mis" in cols else None,
            }
    return out


def load_clingen_dosage(path: Path | str) -> dict[str, dict[str, bool | None]]:
    """Return {gene_symbol: {haploinsufficient, triplosensitive}}.

    Scores ≥ 2 are interpreted as "evidence-backed dosage sensitivity"; 0, 1,
    30 (autosomal recessive), 40 (dosage sensitivity unlikely) and missing
    values collapse to False (or None if the column was absent).
    """
    path = Path(path)
    out: dict[str, dict[str, bool | None]] = {}
    with _open_table(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        cols = _resolve_columns(header, _CLINGEN_COLUMN_ALIASES)
        if "gene" not in cols:
            raise ValueError(f"{path}: no gene column (have: {header})")
        for row in reader:
            if not row:
                continue
            gene = row[cols["gene"]].strip()
            if not gene:
                continue
            hi = _interpret_clingen_score(row[cols["hi"]]) if "hi" in cols else None
            ts = _interpret_clingen_score(row[cols["ts"]]) if "ts" in cols else None
            out[gene] = {"haploinsufficient": hi, "triplosensitive": ts}
    return out


def _interpret_clingen_score(raw: str) -> bool | None:
    v = _parse_float(raw)
    if v is None:
        return None
    return v >= 2 and v < 30


# ─────────────────────────────────────────────────────────────────────────────
# Unified lookup
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class ConstraintLookup:
    """Unified gene → ConstraintPriors lookup."""

    gnomad: dict[str, dict[str, Any]] | None = None
    clingen: dict[str, dict[str, bool | None]] | None = None

    @classmethod
    def load(
        cls,
        gnomad_path: Path | str | None = None,
        clingen_path: Path | str | None = None,
        logger: logging.Logger | None = None,
    ) -> "ConstraintLookup":
        logger = logger or logging.getLogger(__name__)
        gnomad = None
        clingen = None
        if gnomad_path is not None:
            gnomad = load_gnomad_constraint(gnomad_path)
            logger.info("Loaded gnomAD constraint for %d genes from %s", len(gnomad), gnomad_path)
        if clingen_path is not None:
            clingen = load_clingen_dosage(clingen_path)
            logger.info("Loaded ClinGen dosage for %d genes from %s", len(clingen), clingen_path)
        return cls(gnomad=gnomad, clingen=clingen)

    def get(self, gene: str) -> ConstraintPriors:
        """Return priors for a gene. Missing entries yield all-None fields."""
        g = (self.gnomad or {}).get(gene, {})
        c = (self.clingen or {}).get(gene, {})
        return ConstraintPriors(
            pli=g.get("pli"),
            mis_z=g.get("mis_z"),
            oe_lof=g.get("oe_lof"),
            oe_mis=g.get("oe_mis"),
            clingen_haploinsufficient=c.get("haploinsufficient"),
            clingen_triplosensitive=c.get("triplosensitive"),
        )


__all__ = [
    "ConstraintLookup",
    "load_clingen_dosage",
    "load_gnomad_constraint",
]
