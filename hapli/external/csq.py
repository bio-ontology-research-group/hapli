"""
Wrapper around `bcftools csq` (Danecek 2017) for per-haplotype variant
consequence calling.

`bcftools csq -p a` emits a VCF with two annotations per record:

  INFO/BCSQ   — comma-separated list of consequence strings, each of the form
                Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change
                The Consequence may be a compound tag like "missense&inframe_altering"
                when bcftools joins two nearby variants on the same haplotype into
                a single functional unit (e.g. the frameshift + rescue pattern).
                Back-references "@POS" point to the compound record at an earlier
                position — we resolve them to the referenced BCSQ string.

  FORMAT/BCSQ — bitmask of indices into INFO/BCSQ, interleaved first / second
                haplotype. Bit 2k tells you whether INFO/BCSQ[k] applies to hap1;
                bit 2k+1 is the same for hap2.

This wrapper runs `bcftools csq`, parses the resulting records with pysam,
resolves back-references, and returns a list of ConsequenceCall objects keyed
by (sample, haplotype, transcript) that the pipeline can fold into the
per-gene evidence bundle.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pysam

from ..core.schema import ConsequenceCall


DEFAULT_BCFTOOLS = "bcftools"
BCSQ_FIELDS = 7   # Consequence | gene | transcript | biotype | strand | aa_change | dna_change


class BcftoolsNotAvailable(RuntimeError):
    pass


@dataclass
class CsqResult:
    annotated_vcf: Path
    consequences: list[ConsequenceCall]


def _parse_bcsq_entry(raw: str) -> dict | None:
    """Split one BCSQ entry into its 7 fields; return None on malformed input.

    '@POS' back-references are returned as {'backref': POS}.
    """
    if not raw:
        return None
    if raw.startswith("@"):
        try:
            return {"backref": int(raw[1:])}
        except ValueError:
            return None
    parts = raw.split("|")
    # csq emits 7 pipe-separated fields; pad if a tail field was empty.
    while len(parts) < BCSQ_FIELDS:
        parts.append("")
    return {
        "consequence": parts[0],
        "gene": parts[1],
        "transcript": parts[2],
        "biotype": parts[3],
        "strand": parts[4],
        "amino_acid_change": parts[5],
        "dna_change": parts[6],
    }


def _decode_bitmask(mask: int, n_entries: int) -> dict[int, set[int]]:
    """Return {bcsq_index_1based: {haplotype_idx_1based, ...}} from a FORMAT/BCSQ bitmask."""
    out: dict[int, set[int]] = {}
    for k in range(n_entries):
        for hap in (1, 2):
            bit = 2 * k + (hap - 1)
            if mask & (1 << bit):
                out.setdefault(k + 1, set()).add(hap)
    return out


def run_csq(
    vcf_path: Path,
    reference_fasta: Path,
    gff3: Path,
    sample: str,
    out_vcf: Path | None = None,
    bcftools_path: str = DEFAULT_BCFTOOLS,
    logger: logging.Logger | None = None,
) -> CsqResult:
    """Run `bcftools csq -p a` and parse the per-haplotype consequences.

    Parameters
    ----------
    vcf_path        : input bgzipped + tabix-indexed phased VCF.
    reference_fasta : indexed reference FASTA that the GFF coords are relative to.
    gff3            : reference GFF3 annotation.
    sample          : sample name to restrict output to (pulled from FORMAT/BCSQ).
    out_vcf         : where to write the annotated VCF (default: a tempfile).
    bcftools_path   : binary name.
    """
    logger = logger or logging.getLogger(__name__)

    if shutil.which(bcftools_path) is None:
        raise BcftoolsNotAvailable(f"{bcftools_path!r} not on PATH")

    if out_vcf is None:
        tmp = Path(tempfile.mkstemp(suffix=".csq.vcf.gz")[1])
        out_vcf = tmp

    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    # --force: tolerate GFFs with minor phase inconsistencies. bcftools' strict
    # `phase == cds_len % 3` check rejects many real-world GFFs (Ensembl is fine,
    # but hand-built / Liftoff / simulated GFFs commonly have phase=0 everywhere).
    cmd = [
        bcftools_path, "csq",
        "-p", "a",
        "--force",
        "-f", str(reference_fasta),
        "-g", str(gff3),
        "-s", sample,
        "-Oz",
        "-o", str(out_vcf),
        str(vcf_path),
    ]
    logger.debug("Running: %s", " ".join(cmd))
    run = subprocess.run(cmd, capture_output=True, text=True)
    if run.returncode != 0:
        raise RuntimeError(
            f"bcftools csq failed (rc={run.returncode}):\nstderr:\n{run.stderr}"
        )
    # Tabix-index the output so downstream code can do region queries.
    pysam.tabix_index(str(out_vcf), preset="vcf", force=True)

    consequences = _parse_annotated_vcf(out_vcf, sample=sample)
    return CsqResult(annotated_vcf=out_vcf, consequences=consequences)


def _parse_annotated_vcf(vcf: Path, sample: str) -> list[ConsequenceCall]:
    """Walk the csq-annotated VCF and emit one ConsequenceCall per
    (record, BCSQ-entry, haplotype) triple.

    Compound back-references ("@POS") are resolved to the primary entry; the
    `compound_id` field is set to the primary POS so downstream code can group.
    """
    calls: list[ConsequenceCall] = []
    # Primary BCSQ entries keyed by CHROM:POS, for back-ref resolution.
    primaries: dict[tuple[str, int], list[dict]] = {}

    with pysam.VariantFile(str(vcf)) as vf:
        if sample not in vf.header.samples:
            raise ValueError(f"Sample {sample!r} missing from csq output")
        for rec in vf:
            bcsq_info = rec.info.get("BCSQ")
            if not bcsq_info:
                continue
            entries = [_parse_bcsq_entry(e) for e in bcsq_info]
            entries = [e for e in entries if e is not None]

            # First record of a compound block: entries have no 'backref'.
            if all("backref" not in e for e in entries):
                primaries[(rec.chrom, rec.pos)] = entries
            # Back-references: replace with the primary entry.
            resolved = []
            for e in entries:
                if "backref" in e:
                    primary = primaries.get((rec.chrom, e["backref"]))
                    if primary:
                        resolved.extend(primary)
                else:
                    resolved.append(e)
            if not resolved:
                continue

            sample_rec = rec.samples[sample]
            mask = sample_rec.get("BCSQ") or 0
            if isinstance(mask, tuple):
                mask = mask[0] if mask else 0
            per_idx = _decode_bitmask(int(mask), len(resolved))

            for idx, entry in enumerate(resolved, start=1):
                haplotypes = per_idx.get(idx, set())
                if not haplotypes:
                    continue
                # dna_change often carries compound info like "1012T>TC+1029CT>C"
                compound_id = None
                dna = entry.get("dna_change", "")
                if "+" in dna:
                    compound_id = dna
                # Split composite consequences like "missense&inframe_altering"
                consequences = entry["consequence"].split("&") or [entry["consequence"]]
                # Amino-acid change encoding example: "4AAAAAAA>4ARCCCCC"
                aa = entry.get("amino_acid_change", "")
                hgvs_p = aa or None
                hgvs_c = dna or None
                for hap in haplotypes:
                    for cons in consequences:
                        calls.append(
                            ConsequenceCall(
                                chrom=rec.chrom,
                                pos=rec.pos,
                                ref=rec.ref,
                                alt=rec.alts[0] if rec.alts else "",
                                haplotype=hap,
                                transcript=entry.get("transcript", ""),
                                consequence=cons,
                                hgvs_c=hgvs_c,
                                hgvs_p=hgvs_p,
                                compound_id=compound_id,
                            )
                        )
    return calls


def group_by_haplotype(
    calls: Iterable[ConsequenceCall],
) -> dict[int, list[ConsequenceCall]]:
    """Convenience: {1: [...hap1 calls], 2: [...hap2 calls]}."""
    grouped: dict[int, list[ConsequenceCall]] = {1: [], 2: []}
    for c in calls:
        grouped.setdefault(c.haplotype, []).append(c)
    return grouped
