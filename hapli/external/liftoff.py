"""
Wrapper around the Liftoff binary (Shumate & Salzberg 2021).

Liftoff transfers a reference GFF3 onto a target assembly FASTA via minimap2.
We call it once per haplotype to get:

  (a) a per-haplotype GFF3 with valid coordinates on the haplotype (so Phase 2
      protein extraction can slice directly, bypassing the reference-coordinate-
      shift problem introduced by length-changing upstream variants), and
  (b) a per-gene presence call derived from Liftoff's own QC metadata
      (coverage, sequence_ID, valid_ORFs, matches_ref_protein, extra_copy_number)
      plus the unmapped_features.txt file that lists genes absent from the
      haplotype (e.g., removed by an SV DEL).

The binary is installed via `uv tool install git+https://github.com/agshumate/Liftoff`;
the wrapper shells out to the CLI rather than importing the module, so the
project env stays decoupled from Liftoff's pysam pin.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from ..core.schema import PresenceCall


_COPY_SUFFIX = re.compile(r"_\d+$")


DEFAULT_LIFTOFF = "liftoff"


@dataclass
class LiftoffResult:
    """Structured result of one Liftoff run on one haplotype."""

    lifted_gff: Path
    unmapped_txt: Path
    intermediate_dir: Path | None
    # Per-gene presence derived from the lifted GFF + unmapped list.
    presence: dict[str, PresenceCall] = field(default_factory=dict)
    # Convenience: which genes are present on this haplotype
    present_gene_ids: set[str] = field(default_factory=set)
    absent_gene_ids: set[str] = field(default_factory=set)


class LiftoffNotAvailable(RuntimeError):
    pass


def _parse_gff_attrs(attr_str: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in attr_str.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        k, v = item.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def _classify_gene(attrs: dict[str, str]) -> PresenceCall:
    """Translate Liftoff's lifted-GFF3 attributes into a PresenceCall."""
    coverage = float(attrs["coverage"]) if "coverage" in attrs else None
    seq_id = float(attrs["sequence_ID"]) if "sequence_ID" in attrs else None
    valid_orfs = int(attrs["valid_ORFs"]) if "valid_ORFs" in attrs else None
    extra_copy = int(attrs["extra_copy_number"]) if "extra_copy_number" in attrs else None
    matches_ref_protein = attrs.get("matches_ref_protein")
    if matches_ref_protein is not None:
        matches_ref_protein_b: bool | None = matches_ref_protein.lower() == "true"
    else:
        matches_ref_protein_b = None

    flags = {k: v for k, v in attrs.items() if k in ("partial_mapping", "low_identity")}

    if extra_copy is not None and extra_copy > 0:
        status = "duplicated"
    elif flags.get("partial_mapping") == "True" or (coverage is not None and coverage < 0.9):
        status = "partial"
    elif flags.get("low_identity") == "True" or (seq_id is not None and seq_id < 0.95):
        status = "low_identity"
    elif valid_orfs is not None and valid_orfs == 0:
        status = "uncertain"
    elif coverage is not None and coverage >= 0.9 and (seq_id is None or seq_id >= 0.95):
        status = "intact"
    else:
        status = "uncertain"

    return PresenceCall(
        status=status,
        source="liftoff",
        coverage=coverage,
        sequence_identity=seq_id,
        valid_orfs=valid_orfs,
        matches_ref_protein=matches_ref_protein_b,
        extra_copy_number=extra_copy,
        flags=flags,
    )


def parse_lifted_gff(gff_path: Path) -> dict[str, PresenceCall]:
    """Walk a Liftoff-produced GFF3 and emit one PresenceCall per `gene` record.

    Liftoff preserves the reference attributes (Name, biotype, ...) and appends
    `_0`, `_1`, ... to the ID= attribute to disambiguate tandem-duplicated copies.
    We register each call under *both* the base ID (with copy suffix stripped)
    and the Name attribute so downstream code can look up by whichever identifier
    it has. Tandem duplicates bump the call's `extra_copy_number` rather than
    overwriting.
    """
    result: dict[str, PresenceCall] = {}
    if not gff_path.exists():
        return result
    with gff_path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = _parse_gff_attrs(parts[8])
            call = _classify_gene(attrs)
            call.seqid = parts[0]
            call.start = int(parts[3])
            call.end = int(parts[4])
            call.strand = parts[6] if parts[6] in ("+", "-") else None

            gene_id = attrs.get("ID", "")
            base_id = _COPY_SUFFIX.sub("", gene_id) if gene_id else ""
            # GENCODE/Ensembl use `gene_name`; reference GFF3s use `Name`.
            # Register under both so downstream lookup-by-symbol works.
            name = attrs.get("Name", "")
            gene_name = attrs.get("gene_name", "")

            # All identifiers under which this record may be looked up.
            keys = {k for k in (base_id, gene_id, name, gene_name) if k}

            for key in keys:
                if key in result:
                    existing = result[key]
                    existing.extra_copy_number = (existing.extra_copy_number or 0) + 1
                    existing.status = "duplicated"
                else:
                    result[key] = call
    return result


def parse_unmapped(unmapped_path: Path) -> set[str]:
    if not unmapped_path.exists():
        return set()
    return {
        line.strip()
        for line in unmapped_path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }


def run_liftoff(
    haplotype_fasta: Path,
    reference_fasta: Path,
    reference_gff: Path,
    out_gff: Path,
    unmapped_txt: Path | None = None,
    intermediate_dir: Path | None = None,
    polish: bool = False,
    copies: bool = False,
    extra_args: list[str] | None = None,
    liftoff_path: str = DEFAULT_LIFTOFF,
    logger: logging.Logger | None = None,
) -> LiftoffResult:
    """Run Liftoff on one haplotype and parse the result.

    Parameters
    ----------
    haplotype_fasta : indexed FASTA of the target haplotype.
    reference_fasta : indexed reference FASTA that `reference_gff` is coordinated against.
    reference_gff   : reference GFF3 to lift.
    out_gff         : path for the lifted GFF3.
    unmapped_txt    : path for the unmapped-features list (defaults to `out_gff` + ".unmapped.txt").
    intermediate_dir: directory for Liftoff intermediate files (FASTA, SAM).
                      Defaults to `out_gff.parent / f"{out_gff.stem}.liftoff_tmp"`.
    polish          : pass `-polish` (cleans ORF-truncation noise).
    copies          : pass `-copies` (reports extra gene copies).
    extra_args      : extra CLI args passed verbatim.
    liftoff_path    : binary name / path.

    Returns
    -------
    LiftoffResult with parsed per-gene presence and the path to the lifted GFF3.
    """
    logger = logger or logging.getLogger(__name__)

    if shutil.which(liftoff_path) is None:
        raise LiftoffNotAvailable(
            f"{liftoff_path!r} not on PATH. Install via "
            f"`uv tool install git+https://github.com/agshumate/Liftoff`."
        )

    out_gff = Path(out_gff)
    unmapped_txt = Path(unmapped_txt) if unmapped_txt else out_gff.with_suffix(".unmapped.txt")
    intermediate_dir = Path(intermediate_dir) if intermediate_dir else out_gff.parent / f"{out_gff.stem}.liftoff_tmp"

    out_gff.parent.mkdir(parents=True, exist_ok=True)
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        liftoff_path,
        "-g", str(reference_gff),
        "-o", str(out_gff),
        "-u", str(unmapped_txt),
        "-dir", str(intermediate_dir),
    ]
    if polish:
        cmd.append("-polish")
    if copies:
        cmd.append("-copies")
    if extra_args:
        cmd.extend(extra_args)
    cmd.append(str(haplotype_fasta))
    cmd.append(str(reference_fasta))

    logger.debug("Running: %s", " ".join(cmd))
    run = subprocess.run(cmd, capture_output=True, text=True)
    if run.returncode != 0:
        raise RuntimeError(
            f"liftoff failed (rc={run.returncode}):\nstderr:\n{run.stderr}\nstdout:\n{run.stdout}"
        )

    # Liftoff's -polish writes an additional file named `<out_gff>_polished`
    # (the "_polished" is appended to the full filename, not before the extension).
    polished_gff = out_gff.parent / (out_gff.name + "_polished")
    effective_gff = polished_gff if polished_gff.exists() else out_gff

    presence = parse_lifted_gff(effective_gff)
    unmapped = parse_unmapped(unmapped_txt)

    # Mark unmapped genes as deleted.
    for gid in unmapped:
        presence[gid] = PresenceCall(status="deleted", source="liftoff")

    present = {gid for gid, p in presence.items() if p.status != "deleted"}
    absent = set(unmapped)

    return LiftoffResult(
        lifted_gff=effective_gff,
        unmapped_txt=unmapped_txt,
        intermediate_dir=intermediate_dir,
        presence=presence,
        present_gene_ids=present,
        absent_gene_ids=absent,
    )
