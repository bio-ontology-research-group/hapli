"""
Thin wrapper around `samtools faidx <region> | bcftools consensus` for
producing per-haplotype sequences over a region.

Replaces the hand-rolled `HaplotypeBuilder._apply_variants`, which silently
corrupted output on symbolic SVs (e.g. inserted the literal string "<DEL>"
into the haplotype FASTA).
"""

from __future__ import annotations

import contextlib
import fcntl
import logging
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pysam


@contextlib.contextmanager
def _file_lock(lock_path: Path):
    """Exclusive flock-based mutex around shared on-disk caches.

    Required because Snakemake fans hapli analyses out at (sample, gene)
    granularity but the per-haplotype consensus + Liftoff outputs are
    keyed at (sample, hap). Without this, 16 parallel jobs for the same
    sample race-overwrite each other's whole-hap FASTA and lifted GFF,
    producing inconsistent per-gene calls (e.g. one job reading the
    polished GFF while another is still writing it).

    Lock files are tiny sentinels; intentionally not deleted on release
    so they can be re-acquired on subsequent runs without re-creation.
    """
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "w") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


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


def find_overlapping_svs(
    vcf_path: Path,
    sample: str,
    hap_idx: int,
    region: "Region",
    cds_intervals: "list[tuple[int, int]] | None" = None,
) -> list[dict]:
    """Return SV records that the given haplotype carries an ALT for and that
    overlap the gene region (or the CDS intervals if provided).

    Used to build the `causing_svs` evidence chain on non-`intact`
    PresenceCalls so the SV story is machine-readable in downstream TSVs
    and joinable against AnnotSV/SvAnna outputs.

    `hap_idx` is 1 or 2 (matches bcftools consensus -H convention; the
    record's GT is read at position hap_idx-1).
    `cds_intervals` are 1-based inclusive [start, end] pairs in the same
    chromosome as `region`; if None we fall back to gene-region overlap.

    Returns a list of dicts with the CausingSVRecord field shape.
    """
    if hap_idx not in (1, 2):
        raise ValueError(f"hap_idx must be 1 or 2, got {hap_idx}")
    out: list[dict] = []
    cds_total = sum((e - s + 1) for s, e in cds_intervals) if cds_intervals else (region.end - region.start + 1)
    with pysam.VariantFile(str(vcf_path)) as vf:
        if sample not in vf.header.samples:
            raise ValueError(f"Sample {sample!r} not in VCF {vcf_path}")
        for rec in vf.fetch(region.chrom, region.start - 1, region.end):
            alts = rec.alts or ()
            if not alts:
                continue
            sv_type = _svtype(rec) or _classify_alt_as_sv(alts[0])
            if sv_type is None:
                continue
            try:
                end = int(rec.info.get("END") or rec.stop)
            except (ValueError, TypeError, KeyError):
                end = rec.stop
            pos = rec.pos
            length = max(0, end - pos)
            # GT for this hap
            gt = rec.samples[sample].get("GT")
            if gt is None or len(gt) < hap_idx:
                continue
            allele_idx = gt[hap_idx - 1]
            if allele_idx is None or allele_idx == 0:
                # Reference allele on this hap — no SV carried
                continue
            # CDS-overlap fraction
            if cds_intervals:
                overlap = 0
                for s, e in cds_intervals:
                    overlap += max(0, min(e, end) - max(s, pos) + 1)
                ovf = overlap / cds_total if cds_total else 0.0
            else:
                overlap = max(0, min(region.end, end) - max(region.start, pos) + 1)
                ovf = overlap / cds_total if cds_total else 0.0
            if overlap == 0:
                continue
            out.append({
                "chrom": rec.chrom,
                "pos": pos,
                "end": end,
                "sv_type": sv_type,
                "sv_id": rec.id,
                "length": length,
                "overlap_fraction": round(ovf, 4),
                "haplotype": hap_idx,
            })
    return out


def _classify_alt_as_sv(alt: str) -> str | None:
    """Best-effort SVTYPE inference from an ALT string when INFO/SVTYPE is missing."""
    if not alt:
        return None
    if _is_breakend_alt(alt):
        return "BND"
    if _is_inv_alt(alt):
        return "INV"
    if _is_dup_alt(alt):
        return "DUP"
    if alt == "<DEL>" or alt.startswith("<DEL:"):
        return "DEL"
    if alt.startswith("<INS"):
        return "INS"
    if alt == "<CNV>":
        return "CNV"
    return None


def _resolve_symbolic_svs(
    vcf_path: Path,
    reference_fasta: Path,
    region: "Region | None",
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
        if region is None:
            iterator = vf.fetch()
        else:
            iterator = vf.fetch(region.chrom, region.start - 1, region.end)
        for rec in iterator:
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


def consensus_whole_genome(
    reference_fasta: Path,
    vcf_path: Path,
    sample: str,
    output_dir: Path,
) -> dict[str, Path]:
    """Materialise per-haplotype whole-genome FASTAs via `bcftools consensus -f`.

    Why whole-genome and not per-region: when Liftoff runs against a small
    region FASTA, it tries to lift every reference gene onto that fragment
    and silently produces low-confidence chimeric matches for genes that
    don't belong there. With a whole-haplotype FASTA each gene lifts onto
    its real chromosome and the per-(sample, hap) lift cache is correct.

    Output is written to `{output_dir}/{sample}_hap{1,2}.whole.fa(.fai)`
    and re-used on subsequent gene calls for the same (sample, hap).

    Caller must have pre-filtered the VCF to drop unsupported symbolic
    alleles (`<INS:*>`, `<BND>`, `<CNV>`, `<INV>`, `<DUP*>`, `<DEL:*>`).
    Only `<DEL>` symbolic records and explicit-ALT records are tolerated;
    everything else triggers UnsupportedVariantError.

    Returns {'hap1': Path, 'hap2': Path}.
    """
    logger = logging.getLogger(__name__)
    output_dir.mkdir(parents=True, exist_ok=True)

    _preflight_whole(vcf_path, sample)

    # Serialize concurrent Snakemake fan-outs that all want to materialise
    # the same {sample}.whole.fa pair. See _file_lock docstring.
    lock_path = output_dir / f".{sample}.whole.lock"
    with _file_lock(lock_path):
        return _consensus_whole_genome_inner(
            reference_fasta, vcf_path, sample, output_dir, logger
        )


def _consensus_whole_genome_inner(reference_fasta, vcf_path, sample, output_dir, logger):
    # Resolve symbolic <INV>/<DUP*> ALTs into explicit REF/ALT sequences,
    # since bcftools consensus does not handle the symbolic forms.
    with tempfile.TemporaryDirectory(prefix="hapli_svresolve_whole_") as tmp:
        resolved = _resolve_symbolic_svs(vcf_path, reference_fasta, None, Path(tmp))
        effective_vcf = resolved if resolved is not None else vcf_path
        if resolved is not None:
            logger.info("Resolved whole-VCF symbolic <INV>/<DUP*> records → %s", resolved.name)

        out: dict[str, Path] = {}
        for hap_idx in (1, 2):
            out_fa = output_dir / f"{sample}_hap{hap_idx}.whole.fa"
            out_fai = Path(str(out_fa) + ".fai")
            if out_fa.exists() and out_fa.stat().st_size > 0 and out_fai.exists():
                logger.info("Whole-hap FASTA already cached: %s", out_fa.name)
                out[f"hap{hap_idx}"] = out_fa
                continue
            logger.info("Building whole-hap FASTA for sample=%s hap=%d → %s", sample, hap_idx, out_fa.name)
            with out_fa.open("w") as outf:
                proc = subprocess.run(
                    [
                        "bcftools", "consensus",
                        "-f", str(reference_fasta),
                        "-s", sample,
                        "-H", str(hap_idx),
                        str(effective_vcf),
                    ],
                    stdout=outf,
                    stderr=subprocess.PIPE,
                    text=True,
                )
            if proc.returncode != 0:
                out_fa.unlink(missing_ok=True)
                raise RuntimeError(
                    f"bcftools consensus (whole-genome) failed for sample={sample} hap={hap_idx} "
                    f"(rc={proc.returncode}): {proc.stderr.strip()}"
                )
            pysam.faidx(str(out_fa))
            out[f"hap{hap_idx}"] = out_fa
    return out


def _preflight_whole(vcf_path: Path, sample: str) -> None:
    """Whole-VCF preflight — refuse if any unsupported symbolic ALT slipped through."""
    with pysam.VariantFile(str(vcf_path)) as vf:
        if sample not in vf.header.samples:
            raise ValueError(
                f"Sample {sample!r} not in VCF {vcf_path}; "
                f"available samples: {list(vf.header.samples)[:10]}…"
            )
        for rec in vf.fetch():
            for alt in rec.alts or ():
                if _is_breakend_alt(alt):
                    raise UnsupportedVariantError(
                        f"Record at {rec.chrom}:{rec.pos} has breakend ALT {alt!r}. "
                        f"Mode A cannot materialise breakends; pre-filter or use Mode B."
                    )
                if alt == "<CNV>":
                    raise UnsupportedVariantError(
                        f"Record at {rec.chrom}:{rec.pos} is <CNV>; pre-filter the VCF."
                    )
                if alt.startswith("<INS") and alt != "<INS>":
                    # <INS:ME:ALU> etc.
                    raise UnsupportedVariantError(
                        f"Record at {rec.chrom}:{rec.pos} is symbolic mobile-element insertion {alt!r}; "
                        f"pre-filter the VCF (e.g. `bcftools view -e 'ALT~\"<INS\"'`)."
                    )
                if alt == "<INS>":
                    raise UnsupportedVariantError(
                        f"Record at {rec.chrom}:{rec.pos} is symbolic <INS> without explicit sequence."
                    )
                if alt.startswith("<DEL:"):
                    raise UnsupportedVariantError(
                        f"Record at {rec.chrom}:{rec.pos} is parameterized symbolic deletion {alt!r}; "
                        f"pre-filter the VCF."
                    )
                # <DUP*> and <INV> are tolerated: consensus_whole_genome resolves
                # them into explicit REF/ALT sequence form before bcftools consensus.
