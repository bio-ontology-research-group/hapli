"""
Lookup into the precomputed AlphaMissense per-missense-variant score table
(Cheng et al. 2023).

The table is distributed as a bgzipped + tabix-indexed TSV keyed on
(CHROM, POS, REF, ALT). One lookup per missense variant; pure I/O, no model.
Canonical download: https://zenodo.org/records/8208688 (hg38 scores ~11 GB).

The table schema we consume (columns after `#CHROM`):
    CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id  protein_variant
    am_pathogenicity  am_class

`am_pathogenicity` is a float in [0, 1]. `am_class` is one of
`likely_benign`, `ambiguous`, `likely_pathogenic`.

We intentionally support *any* table with that column layout so that
resource-constrained users can ship a gene-subset or transcript-subset
instead of the full 11 GB file.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pysam


DEFAULT_COLUMNS = (
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "genome",
    "uniprot_id",
    "transcript_id",
    "protein_variant",
    "am_pathogenicity",
    "am_class",
)


@dataclass(frozen=True)
class AlphaMissenseHit:
    chrom: str
    pos: int
    ref: str
    alt: str
    uniprot_id: str
    transcript_id: str
    protein_variant: str
    score: float
    am_class: str


class AlphaMissenseLookup:
    """Random-access AlphaMissense lookup over a tabix-indexed TSV."""

    def __init__(
        self,
        table_path: Path | str,
        columns: Iterable[str] = DEFAULT_COLUMNS,
        chr_prefix: str | None = None,
        logger: logging.Logger | None = None,
    ):
        """
        Parameters
        ----------
        table_path   : bgzipped TSV (with a `.tbi` sidecar).
        columns      : column names in the order they appear in the file. Must include
                       CHROM, POS, REF, ALT, am_pathogenicity, am_class.
        chr_prefix   : optional prefix to add/strip from the query CHROM to match the
                       table. Typical values: 'chr' (query chr1 → file chr1) or ''
                       (query chr1 → file 1). Auto-detected on first lookup if None.
        """
        self.table_path = Path(table_path)
        if not self.table_path.exists():
            raise FileNotFoundError(self.table_path)
        self.columns = list(columns)
        self._col_ix = {name: i for i, name in enumerate(self.columns)}
        for req in ("CHROM", "POS", "REF", "ALT", "am_pathogenicity", "am_class"):
            if req not in self._col_ix:
                raise ValueError(f"columns missing required field {req!r}: {self.columns}")
        self.chr_prefix = chr_prefix
        self.logger = logger or logging.getLogger(__name__)
        self._tabix = pysam.TabixFile(str(self.table_path))

    def close(self) -> None:
        self._tabix.close()

    def __enter__(self) -> "AlphaMissenseLookup":
        return self

    def __exit__(self, *_exc) -> None:
        self.close()

    def _normalise_chrom(self, chrom: str) -> str:
        if self.chr_prefix is None:
            # Auto-detect on first call by peeking at the tabix header contigs.
            try:
                contigs = set(self._tabix.contigs)
                has_chr = any(c.startswith("chr") for c in contigs)
                self.chr_prefix = "chr" if has_chr else ""
            except Exception:
                self.chr_prefix = ""
        if self.chr_prefix == "chr":
            return chrom if chrom.startswith("chr") else f"chr{chrom}"
        # prefix-stripping mode
        return chrom[3:] if chrom.startswith("chr") else chrom

    def get(self, chrom: str, pos: int, ref: str, alt: str) -> AlphaMissenseHit | None:
        """Return the AlphaMissense entry for this (chrom, pos, ref, alt), or None."""
        norm_chrom = self._normalise_chrom(chrom)
        try:
            rows = self._tabix.fetch(norm_chrom, pos - 1, pos)
        except ValueError:
            return None  # contig absent
        for row in rows:
            fields = row.rstrip("\n").split("\t")
            if len(fields) < len(self.columns):
                continue
            if int(fields[self._col_ix["POS"]]) != pos:
                continue
            if fields[self._col_ix["REF"]] != ref:
                continue
            if fields[self._col_ix["ALT"]] != alt:
                continue
            try:
                score = float(fields[self._col_ix["am_pathogenicity"]])
            except ValueError:
                continue
            return AlphaMissenseHit(
                chrom=fields[self._col_ix["CHROM"]],
                pos=int(fields[self._col_ix["POS"]]),
                ref=fields[self._col_ix["REF"]],
                alt=fields[self._col_ix["ALT"]],
                uniprot_id=fields[self._col_ix.get("uniprot_id", 0)] if "uniprot_id" in self._col_ix else "",
                transcript_id=fields[self._col_ix.get("transcript_id", 0)] if "transcript_id" in self._col_ix else "",
                protein_variant=fields[self._col_ix.get("protein_variant", 0)] if "protein_variant" in self._col_ix else "",
                score=score,
                am_class=fields[self._col_ix["am_class"]],
            )
        return None

    def get_many(
        self,
        variants: Iterable[tuple[str, int, str, str]],
    ) -> dict[tuple[str, int, str, str], AlphaMissenseHit | None]:
        """Batch lookup. Input: iterable of (chrom, pos, ref, alt) tuples."""
        out: dict[tuple[str, int, str, str], AlphaMissenseHit | None] = {}
        for v in variants:
            out[v] = self.get(*v)
        return out


def download_hint() -> str:
    """Human-readable message telling the user where to get the real AlphaMissense table."""
    return (
        "AlphaMissense hg38 precomputed scores are a ~11 GB bgzipped TSV hosted on "
        "Zenodo: https://zenodo.org/records/8208688 (see AlphaMissense_hg38.tsv.gz). "
        "After download, run `tabix -s 1 -b 2 -e 2 -S 1 AlphaMissense_hg38.tsv.gz` "
        "to create the .tbi sidecar. Point hapli at the file via --am-table or the "
        "HAPLI_ALPHAMISSENSE environment variable."
    )
