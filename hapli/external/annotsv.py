"""
Wrapper around AnnotSV (Geoffroy 2018, https://lbgi.fr/AnnotSV/) for
SV-aware gene annotation. Built for the `hapli` head-to-head benchmark:
hapli emits per-haplotype protein-level presence/LoF calls; AnnotSV
emits per-SV gene-impact + ACMG-class calls. Joining them gene-by-gene
on the same VCF is the comparison we publish.

AnnotSV reads a VCF and writes a TSV with one row per (SV, overlapping
gene). It does not run a per-haplotype protein extraction — its
"impact" call is a coordinate-based overlap classification (full / split
read-through / partial CDS / etc.) plus an ACMG class.

Install:
    conda install -c bioconda annotsv
    # then download the annotation bundle once:
    AnnotSV -annotationsDir $ANNOTSV/share/AnnotSV/Annotations_Human/
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional


DEFAULT_ANNOTSV = "AnnotSV"


class AnnotSVNotAvailable(RuntimeError):
    pass


@dataclass
class AnnotSVCall:
    """One row of the AnnotSV per-SV-per-gene TSV.

    Field names match the AnnotSV TSV columns. The full TSV has 100+
    columns; we surface only what the head-to-head needs. `extra`
    keeps the rest verbatim so reviewer requests for additional
    columns don't require re-running AnnotSV.
    """
    sv_chrom: str
    sv_start: int
    sv_end: int
    sv_type: str                     # DEL / DUP / INV / INS / BND / TRA / CNV
    sv_id: Optional[str]             # VCF ID column or AnnotSV-assigned
    gene_name: str                   # AnnotSV's `Gene_name` column
    location: str                    # full CDS / exonic / intronic / etc.
    annotation_mode: str             # "full" or "split" (per-gene split rows)
    acmg_class: Optional[int]        # 1-5 (Riggs 2020 SV ACMG classes)
    extra: dict[str, str] = field(default_factory=dict)


@dataclass
class AnnotSVResult:
    annotated_tsv: Path
    per_sv_calls: list[AnnotSVCall]


def run_annotsv(
    vcf_path: Path,
    *,
    genome_build: str = "GRCh38",
    annotation_dir: Optional[Path] = None,
    out_tsv: Optional[Path] = None,
    annotsv_path: str = DEFAULT_ANNOTSV,
    logger: Optional[logging.Logger] = None,
) -> AnnotSVResult:
    """Run AnnotSV on `vcf_path`, parse the resulting TSV, return calls.

    Parameters
    ----------
    vcf_path : input VCF (bgzip + tabix recommended; AnnotSV will
                accept plain VCF too).
    genome_build : "GRCh38" or "GRCh37".
    annotation_dir : optional override for $ANNOTSV/share/AnnotSV/Annotations_Human.
    out_tsv : where to write the output. Defaults to vcf_path with `.annotsv.tsv`.
    annotsv_path : binary name (resolved on PATH).

    Raises
    ------
    AnnotSVNotAvailable : if `AnnotSV` is not on PATH.
    RuntimeError        : if the AnnotSV invocation exits non-zero.
    """
    log = logger or logging.getLogger(__name__)
    if shutil.which(annotsv_path) is None:
        raise AnnotSVNotAvailable(
            f"{annotsv_path!r} not on PATH. Install via "
            f"`conda install -c bioconda annotsv` and ensure the binary is reachable."
        )
    if out_tsv is None:
        out_tsv = vcf_path.with_suffix(vcf_path.suffix + ".annotsv.tsv")
    cmd = [
        annotsv_path,
        "-SVinputFile", str(vcf_path),
        "-outputFile", str(out_tsv),
        "-genomeBuild", genome_build,
    ]
    if annotation_dir:
        cmd += ["-annotationsDir", str(annotation_dir)]
    log.info("Running AnnotSV: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"AnnotSV failed (rc={proc.returncode}):\n"
            f"stderr:\n{proc.stderr}\nstdout:\n{proc.stdout}"
        )
    calls = list(parse_annotsv_tsv(out_tsv))
    log.info("AnnotSV: parsed %d (SV, gene) rows from %s", len(calls), out_tsv.name)
    return AnnotSVResult(annotated_tsv=out_tsv, per_sv_calls=calls)


# Required AnnotSV columns we surface as typed fields. The rest go into
# `extra`. AnnotSV column names are case-sensitive and stable across
# v3.x; we tolerate either underscore or space variants.
_REQUIRED_COLS = {
    "sv_chrom":     ("SV_chrom",),
    "sv_start":     ("SV_start",),
    "sv_end":       ("SV_end",),
    "sv_type":      ("SV_type",),
    "sv_id":        ("SV_ID", "AnnotSV_ID"),
    "gene_name":    ("Gene_name",),
    "location":     ("Location",),
    "annotation_mode": ("Annotation_mode",),
    "acmg_class":   ("ACMG_class", "AnnotSV_ranking_score"),
}


def parse_annotsv_tsv(tsv_path: Path) -> Iterable[AnnotSVCall]:
    """Parse the AnnotSV TSV — one row per (SV, gene) for split mode,
    or one row per SV for full mode. Yields AnnotSVCall instances.

    Tolerates both v3.3 and v3.4 column naming. Unknown columns go
    into `extra` so downstream consumers can still access them.
    """
    with tsv_path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        col_idx = {name: i for i, name in enumerate(header)}
        # Resolve required-field column indices from the alias map
        idx_map: dict[str, Optional[int]] = {}
        for field_name, aliases in _REQUIRED_COLS.items():
            idx_map[field_name] = next(
                (col_idx[a] for a in aliases if a in col_idx), None
            )
        if idx_map["gene_name"] is None or idx_map["sv_chrom"] is None:
            raise ValueError(
                f"AnnotSV TSV {tsv_path} missing essential columns "
                f"(Gene_name, SV_chrom). Found: {header[:10]}…"
            )
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            try:
                sv_start = int(parts[idx_map["sv_start"]]) if idx_map["sv_start"] is not None else 0
                sv_end = int(parts[idx_map["sv_end"]]) if idx_map["sv_end"] is not None else 0
            except ValueError:
                continue
            acmg_raw = parts[idx_map["acmg_class"]] if idx_map["acmg_class"] is not None else ""
            try:
                acmg = int(acmg_raw) if acmg_raw and acmg_raw != "NA" else None
            except ValueError:
                acmg = None
            extra = {
                header[i]: parts[i]
                for i in range(len(header))
                if header[i] not in {a for aliases in _REQUIRED_COLS.values() for a in aliases}
            }
            yield AnnotSVCall(
                sv_chrom=parts[idx_map["sv_chrom"]] if idx_map["sv_chrom"] is not None else "",
                sv_start=sv_start,
                sv_end=sv_end,
                sv_type=parts[idx_map["sv_type"]] if idx_map["sv_type"] is not None else "",
                sv_id=(parts[idx_map["sv_id"]] or None) if idx_map["sv_id"] is not None else None,
                gene_name=parts[idx_map["gene_name"]],
                location=parts[idx_map["location"]] if idx_map["location"] is not None else "",
                annotation_mode=parts[idx_map["annotation_mode"]] if idx_map["annotation_mode"] is not None else "",
                acmg_class=acmg,
                extra=extra,
            )
