"""
Thin wrapper around `samtools faidx <region> | bcftools consensus` for
producing per-haplotype sequences over a region.

Replaces the hand-rolled `HaplotypeBuilder._apply_variants`, which silently
corrupted output on symbolic SVs (e.g. inserted the literal string "<DEL>"
into the haplotype FASTA).
"""

from __future__ import annotations

import logging
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pysam


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


REGION_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")


class UnsupportedVariantError(ValueError):
    """Raised when the region contains a variant class Mode A cannot handle."""


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive

    @classmethod
    def parse(cls, region: str) -> "Region":
        m = REGION_RE.match(region)
        if not m:
            raise ValueError(f"Region must look like 'chr:start-end', got {region!r}")
        return cls(m["chrom"], int(m["start"]), int(m["end"]))

    def __str__(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


def _is_breakend_alt(alt: str) -> bool:
    # Breakend ALTs look like "N[chr2:12345[" or "]chr2:12345]N" or "<BND>".
    return alt == "<BND>" or "[" in alt or "]" in alt


def _preflight(vcf_path: Path, sample: str, region: Region) -> None:
    """Raise UnsupportedVariantError if the region contains records Mode A can't materialise.

    Catches: BND breakends, CNV records, and symbolic INS without explicit ALT sequence.
    Everything else (SNV, small indel, <DEL>, <INV>, <DUP>) bcftools consensus handles.
    """
    with pysam.VariantFile(str(vcf_path)) as vf:
        if sample not in vf.header.samples:
            raise ValueError(
                f"Sample {sample!r} not in VCF {vcf_path}; "
                f"available samples: {list(vf.header.samples)}"
            )
        for rec in vf.fetch(region.chrom, region.start - 1, region.end):
            for alt in rec.alts or ():
                if _is_breakend_alt(alt):
                    raise UnsupportedVariantError(
                        f"Record at {rec.chrom}:{rec.pos} has breakend ALT {alt!r}. "
                        f"Mode A cannot materialise breakends in a per-chromosome consensus; "
                        f"use Mode B (pre-assembled haplotypes)."
                    )
                if alt == "<CNV>":
                    raise UnsupportedVariantError(
                        f"Record at {rec.chrom}:{rec.pos} is a <CNV>. "
                        f"Copy-number is allele-count, not allele-state; not representable in Mode A."
                    )
                if alt == "<INS>":
                    raise UnsupportedVariantError(
                        f"Record at {rec.chrom}:{rec.pos} is a symbolic <INS> without explicit sequence. "
                        f"Resolve the inserted sequence upstream (pangenome graph / long-read caller) before running hapli."
                    )


def _is_inv_alt(alt: str) -> bool:
    return alt == "<INV>"


def _is_dup_alt(alt: str) -> bool:
    # <DUP>, <DUP:TANDEM>, <DUP:INT> all treated as tandem-duplicating.
    return alt == "<DUP>" or alt.startswith("<DUP:")


def _svtype(rec: "pysam.VariantRecord") -> str | None:
    try:
        t = rec.info.get("SVTYPE")
    except (ValueError, KeyError):
        return None
    if isinstance(t, tuple):
        return t[0] if t else None
    return t


def _resolve_symbolic_svs(
    vcf_path: Path,
    reference_fasta: Path,
    region: "Region",
    out_dir: Path,
) -> Path | None:
    """If the region contains `<INV>` or `<DUP*>` records, write a copy of the VCF
    where those records are rewritten into explicit REF/ALT sequence alleles.
    bcftools consensus (through at least v1.21) does not handle these symbolic
    ALTs, but it handles the resolved equivalents fine.

    Returns the path to the rewritten, bgzipped+tabixed VCF, or None if no
    resolution was needed.

    Convention for resolution: the affected range is 1-based inclusive [POS..END],
    as our synthetic generator and common long-read callers produce. We rewrite as:
      <INV>:  REF = ref[POS-1:END], ALT = reverse_complement(REF)    (same length)
      <DUP*>: REF = ref[POS-1:END], ALT = REF + REF                  (tandem)
    """
    fa = pysam.FastaFile(str(reference_fasta))
    needs = False
    records_out: list = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        header = vf.header.copy()
        for rec in vf.fetch(region.chrom, region.start - 1, region.end):
            alts = rec.alts or ()
            if any(_is_inv_alt(a) or _is_dup_alt(a) for a in alts):
                needs = True
                break
    if not needs:
        return None

    resolved_path = out_dir / "resolved.vcf.gz"
    uncompressed = out_dir / "resolved.vcf"
    with pysam.VariantFile(str(vcf_path)) as vf:
        # Materialise a rewritten VCF by copying the header + all records,
        # replacing symbolic INV/DUP ALTs with explicit sequences.
        with uncompressed.open("w") as fout:
            # header
            for line in str(vf.header).splitlines():
                if not line:
                    continue
                fout.write(line + "\n")
            for rec in vf.fetch():
                alts = rec.alts or ()
                # Non-symbolic records: emit unchanged.
                if not any(_is_inv_alt(a) or _is_dup_alt(a) for a in alts):
                    fout.write(str(rec))
                    continue
                # Need END from INFO.
                try:
                    end = int(rec.info.get("END") or rec.stop)
                except (ValueError, TypeError, KeyError):
                    end = rec.stop
                pos = rec.pos
                ref_seq = fa.fetch(rec.chrom, pos - 1, end)
                # Emit one resolved record per ALT (preserves genotype indexing
                # because we replace all symbolic ALTs at the same position).
                new_alts = []
                for a in alts:
                    if _is_inv_alt(a):
                        new_alts.append(_reverse_complement(ref_seq))
                    elif _is_dup_alt(a):
                        new_alts.append(ref_seq + ref_seq)
                    else:
                        # Mixed symbolic + non-symbolic ALTs at one site are
                        # uncommon; pass through verbatim.
                        new_alts.append(a)
                # Build a minimal VCF line (tab-separated) preserving sample data.
                sample_fields = []
                for sname in vf.header.samples:
                    s = rec.samples[sname]
                    gt = s.get("GT")
                    if gt is None:
                        sample_fields.append(".")
                    else:
                        sep = "|" if s.phased else "/"
                        sample_fields.append(sep.join(str(g) if g is not None else "." for g in gt))
                # Drop END from INFO — the resolved record no longer needs it
                # and leaving END in place confuses bcftools about the new length.
                info_parts = []
                for k, v in rec.info.items():
                    if k == "END":
                        continue
                    if isinstance(v, tuple):
                        v = ",".join(str(x) for x in v)
                    if v is True:
                        info_parts.append(k)
                    else:
                        info_parts.append(f"{k}={v}")
                info_str = ";".join(info_parts) if info_parts else "."
                cols = [
                    rec.chrom, str(pos), rec.id or ".",
                    ref_seq, ",".join(new_alts),
                    str(rec.qual) if rec.qual is not None else ".",
                    ";".join(rec.filter.keys()) or ".",
                    info_str,
                    "GT",
                    *sample_fields,
                ]
                fout.write("\t".join(cols) + "\n")
    # bgzip + tabix
    pysam.tabix_compress(str(uncompressed), str(resolved_path), force=True)
    pysam.tabix_index(str(resolved_path), preset="vcf", force=True)
    uncompressed.unlink()
    return resolved_path


def consensus_region(
    reference_fasta: Path,
    vcf_path: Path,
    sample: str,
    region: str,
) -> dict[str, str]:
    """Materialise hap1/hap2 sequences over the given region via bcftools consensus.

    Parameters
    ----------
    reference_fasta : indexed reference FASTA (`.fai` sidecar required).
    vcf_path        : bgzipped, tabix-indexed, phased VCF.
    sample          : sample name present in the VCF.
    region          : 1-based inclusive region string "chr:start-end".

    Returns
    -------
    {'hap1': <seq>, 'hap2': <seq>}

    Raises
    ------
    UnsupportedVariantError : if the region contains BND, <CNV>, or symbolic <INS>.
    RuntimeError            : if either bcftools invocation fails.
    """
    logger = logging.getLogger(__name__)
    parsed = Region.parse(region)

    _preflight(vcf_path, sample, parsed)

    # bcftools consensus (≤1.21) does not handle symbolic <INV> / <DUP*> ALTs.
    # Rewrite any such records into explicit REF/ALT sequence form first.
    with tempfile.TemporaryDirectory(prefix="hapli_svresolve_") as tmp:
        resolved = _resolve_symbolic_svs(vcf_path, reference_fasta, parsed, Path(tmp))
        effective_vcf = resolved if resolved is not None else vcf_path
        if resolved is not None:
            logger.info("Resolved symbolic <INV>/<DUP*> records into explicit ALTs: %s", resolved)

        haps: dict[str, str] = {}
        for hap_idx in (1, 2):
            faidx = subprocess.Popen(
                ["samtools", "faidx", str(reference_fasta), str(parsed)],
                stdout=subprocess.PIPE,
            )
            consensus = subprocess.run(
                [
                    "bcftools", "consensus",
                    "--haplotype", str(hap_idx),
                    "--samples", sample,
                    str(effective_vcf),
                ],
                stdin=faidx.stdout,
                capture_output=True,
                text=True,
            )
            faidx.stdout.close() if faidx.stdout else None
            faidx_rc = faidx.wait()
            if faidx_rc != 0:
                raise RuntimeError(f"samtools faidx failed (rc={faidx_rc}) for region {parsed}")
            if consensus.returncode != 0:
                raise RuntimeError(
                    f"bcftools consensus failed (rc={consensus.returncode}) for hap{hap_idx}: "
                    f"{consensus.stderr.strip()}"
                )

            seq = _fasta_to_seq(consensus.stdout)
            logger.debug("hap%d length %d for %s sample=%s", hap_idx, len(seq), parsed, sample)
            haps[f"hap{hap_idx}"] = seq

        return haps


def _fasta_to_seq(fasta_text: str) -> str:
    """Concatenate all non-header lines of a single-record FASTA into one string."""
    return "".join(line for line in fasta_text.splitlines() if line and not line.startswith(">"))
