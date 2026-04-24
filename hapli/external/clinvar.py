"""
ClinVar cross-reference for hapli.

Tabix-indexed lookup over the ClinVar variant_summary.txt.gz (or the standard
ClinVar VCF), keyed by (CHROM, POS, REF, ALT). Returns the curated clinical
significance label so downstream code can annotate per-variant ConsequenceCalls
with `pathogenic` / `likely_pathogenic` / `benign` / `vus` / etc.

Canonical download:
  ClinVar VCF (GRCh38):
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
  ClinVar variant summary TSV:
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

Both are bgzipped and tabix-indexable. We expose two ingestion paths:

  ClinVarLookup.from_vcf(path) — uses pysam.VariantFile.fetch(); CHROM/POS/REF/ALT
                                 lookup gives back CLNSIG from INFO.
  ClinVarLookup.from_tsv(path) — uses tabix on a (CHROM, POS) prefix; designed
                                 for variant_summary.txt.gz.

For users who want a smaller subset, both ingestion paths accept any custom
file with the same column structure (subset by gene / chromosome before tabix).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pysam


# ClinVar's Clinical Significance vocabulary. Normalize to a small canonical set.
_CLNSIG_NORMAL = {
    "pathogenic": "pathogenic",
    "likely_pathogenic": "likely_pathogenic",
    "pathogenic/likely_pathogenic": "pathogenic",
    "uncertain_significance": "vus",
    "likely_benign": "likely_benign",
    "benign": "benign",
    "benign/likely_benign": "benign",
    "drug_response": "drug_response",
    "association": "association",
    "risk_factor": "risk_factor",
    "conflicting_interpretations_of_pathogenicity": "conflicting",
    "not_provided": "not_provided",
    "no_assertion_provided": "not_provided",
    "no_assertion_criteria_provided": "not_provided",
    "no_classifications_from_unflagged_records": "not_provided",
}


def _normalise_clnsig(raw: str | None) -> str | None:
    if not raw:
        return None
    # ClinVar uses both '_' and ' ' as separators; collapse and lowercase.
    s = raw.lower().replace(" ", "_").strip().rstrip(",")
    # Multi-significance entries (comma-separated): pick the most-pathogenic.
    if "," in s:
        parts = [p.strip() for p in s.split(",")]
    elif "|" in s:
        parts = [p.strip() for p in s.split("|")]
    else:
        parts = [s]
    # Priority: pathogenic > likely_pathogenic > vus > likely_benign > benign
    priority = ("pathogenic", "likely_pathogenic", "vus",
                "likely_benign", "benign", "drug_response", "association",
                "risk_factor", "conflicting", "not_provided")
    candidates = [_CLNSIG_NORMAL.get(p, p) for p in parts]
    for level in priority:
        if level in candidates:
            return level
    return candidates[0] if candidates else None


@dataclass(frozen=True)
class ClinVarHit:
    chrom: str
    pos: int
    ref: str
    alt: str
    clnsig: str | None
    raw_clnsig: str | None
    review_status: str | None = None
    gene_symbol: str | None = None
    rcv: str | None = None


class ClinVarLookup:
    """Random-access ClinVar lookup by (chrom, pos, ref, alt)."""

    def __init__(
        self,
        path: Path | str,
        kind: str,                                  # "vcf" or "tsv"
        chr_prefix: str | None = None,
        logger: logging.Logger | None = None,
    ):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(self.path)
        self.kind = kind
        self.chr_prefix = chr_prefix
        self.logger = logger or logging.getLogger(__name__)
        if kind == "vcf":
            self._vf = pysam.VariantFile(str(self.path))
            self._tabix = None
        elif kind == "tsv":
            self._vf = None
            self._tabix = pysam.TabixFile(str(self.path))
        else:
            raise ValueError(f"unknown kind {kind!r}; use 'vcf' or 'tsv'")

    @classmethod
    def from_vcf(cls, path: Path | str, **kwargs) -> "ClinVarLookup":
        return cls(path, kind="vcf", **kwargs)

    @classmethod
    def from_tsv(cls, path: Path | str, **kwargs) -> "ClinVarLookup":
        return cls(path, kind="tsv", **kwargs)

    def close(self) -> None:
        if self._vf is not None: self._vf.close()
        if self._tabix is not None: self._tabix.close()

    def __enter__(self) -> "ClinVarLookup":
        return self

    def __exit__(self, *_exc) -> None:
        self.close()

    def _normalise_chrom(self, chrom: str) -> str:
        # ClinVar VCF uses '1' (no chr prefix); query may use 'chr1'.
        if self.chr_prefix is None:
            try:
                if self.kind == "vcf":
                    contigs = set(self._vf.header.contigs.keys())
                else:
                    contigs = set(self._tabix.contigs)
                self.chr_prefix = "chr" if any(c.startswith("chr") for c in contigs) else ""
            except Exception:
                self.chr_prefix = ""
        if self.chr_prefix == "chr":
            return chrom if chrom.startswith("chr") else f"chr{chrom}"
        return chrom[3:] if chrom.startswith("chr") else chrom

    def get(self, chrom: str, pos: int, ref: str, alt: str) -> ClinVarHit | None:
        norm = self._normalise_chrom(chrom)
        if self.kind == "vcf":
            return self._get_vcf(norm, pos, ref, alt)
        return self._get_tsv(norm, pos, ref, alt)

    def _get_vcf(self, chrom: str, pos: int, ref: str, alt: str) -> ClinVarHit | None:
        try:
            recs = self._vf.fetch(chrom, pos - 1, pos)
        except (ValueError, KeyError):
            return None
        for rec in recs:
            if rec.pos != pos or rec.ref != ref:
                continue
            if alt not in (rec.alts or ()):
                continue
            info = rec.info
            def _info(key, default=None):
                """pysam .info.get() raises ValueError when the key isn't in the
                VCF header — common on subset/synthetic ClinVar files. Wrap so
                missing fields just return None."""
                try:
                    return info.get(key, default)
                except (ValueError, KeyError):
                    return default
            raw = _info("CLNSIG")
            if isinstance(raw, tuple):
                raw = ",".join(str(x) for x in raw)
            review = _info("CLNREVSTAT")
            if isinstance(review, tuple):
                review = ",".join(str(x) for x in review)
            gi = _info("GENEINFO")
            gene = gi.split(":", 1)[0] if gi else None
            rcv = _info("CLNACC")
            return ClinVarHit(
                chrom=chrom, pos=pos, ref=ref, alt=alt,
                clnsig=_normalise_clnsig(raw),
                raw_clnsig=raw,
                review_status=str(review) if review else None,
                gene_symbol=gene,
                rcv=str(rcv) if rcv else None,
            )
        return None

    def _get_tsv(self, chrom: str, pos: int, ref: str, alt: str) -> ClinVarHit | None:
        # variant_summary.txt.gz is tab-delimited; columns we care about are
        # GeneSymbol(4), ClinicalSignificance(6), ReviewStatus(24),
        # Chromosome(18), Start(19), ReferenceAlleleVCF(31), AlternateAlleleVCF(32),
        # Assembly(16). Indices are 0-based.
        try:
            rows = self._tabix.fetch(chrom, pos - 1, pos)
        except (ValueError, KeyError):
            return None
        for row in rows:
            f = row.split("\t")
            if len(f) < 33:
                continue
            try:
                if int(f[19]) != pos:
                    continue
            except (ValueError, IndexError):
                continue
            if f[31] != ref or f[32] != alt:
                continue
            return ClinVarHit(
                chrom=chrom, pos=pos, ref=ref, alt=alt,
                clnsig=_normalise_clnsig(f[6]),
                raw_clnsig=f[6],
                review_status=f[24],
                gene_symbol=f[4],
                rcv=None,
            )
        return None


__all__ = ["ClinVarHit", "ClinVarLookup"]
