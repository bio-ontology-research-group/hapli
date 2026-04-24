"""
Haplotype-resolved protein extraction.

Given a haplotype FASTA (from bcftools consensus or a de novo assembly) and a
GFF3 lifted onto *that same haplotype* (from Liftoff), produce the translated
protein sequence of each transcript.

Critical behaviour encoded here:

* The reference-GFF3 feature hierarchy (gene → mRNA → CDS) is preserved by
  Liftoff, so we walk the same way on the lifted GFF. The coordinates are
  valid on the haplotype FASTA — no re-mapping needed.
* Reverse-strand transcripts are handled correctly: extract every CDS on the
  + strand of the haplotype, concatenate in genomic order, then reverse-
  complement the *whole* concatenated sequence. The old pipeline bypassed
  `SequenceExtractor.get_sequence` (which would RC each exon individually,
  producing the wrong spliced sequence) by inlining this logic; it now lives
  in `splice_cds()` once, with the Phase 0 fix preserved as the canonical
  implementation.
* Translation tolerates partial-codon tails (Liftoff's lifted CDS sometimes
  drops a trailing 1–2 bp when polish is enabled); the remainder is logged
  but does not abort.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import pysam
from Bio.Seq import Seq


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement; preserves case and leaves N/n alone."""
    return seq.translate(_COMPLEMENT)[::-1]


# ─────────────────────────────────────────────────────────────────────────────
# Lightweight GFF view
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class GffRecord:
    seqid: str
    source: str
    featuretype: str
    start: int                 # 1-based inclusive
    end: int                   # 1-based inclusive
    strand: str
    phase: str
    attributes: dict[str, str]

    @property
    def id(self) -> str:
        return self.attributes.get("ID", "")

    @property
    def parent(self) -> str:
        return self.attributes.get("Parent", "")


def _parse_attrs(attr_str: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in attr_str.strip().rstrip(";").split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def load_gff(gff_path: Path) -> list[GffRecord]:
    """Read a (reference or lifted) GFF3 into a list of GffRecord objects."""
    records: list[GffRecord] = []
    with Path(gff_path).open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            try:
                records.append(
                    GffRecord(
                        seqid=parts[0],
                        source=parts[1],
                        featuretype=parts[2],
                        start=int(parts[3]),
                        end=int(parts[4]),
                        strand=parts[6],
                        phase=parts[7],
                        attributes=_parse_attrs(parts[8]),
                    )
                )
            except ValueError:
                continue
    return records


def transcripts_of_gene(records: list[GffRecord], gene_id: str) -> list[GffRecord]:
    """Return mRNA records whose Parent is `gene_id` (or whose gene ID matches after
    stripping Liftoff's `_N` copy suffix).
    """
    import re
    gene_base = re.sub(r"_\d+$", "", gene_id)

    # Find the gene record(s) — either direct ID match or base-ID match.
    gene_ids_here = {
        r.id for r in records
        if r.featuretype == "gene"
        and (r.id == gene_id or re.sub(r"_\d+$", "", r.id) == gene_base
             or r.attributes.get("Name") == gene_id)
    }
    if not gene_ids_here:
        return []
    return [
        r for r in records
        if r.featuretype == "mRNA" and r.parent in gene_ids_here
    ]


def cds_of_transcript(records: list[GffRecord], transcript_id: str) -> list[GffRecord]:
    """Return the CDS records whose Parent is `transcript_id`, sorted by genomic start."""
    cds = [r for r in records if r.featuretype == "CDS" and r.parent == transcript_id]
    cds.sort(key=lambda r: r.start)
    return cds


# ─────────────────────────────────────────────────────────────────────────────
# Splicing + translation
# ─────────────────────────────────────────────────────────────────────────────
def splice_cds(
    fasta: pysam.FastaFile,
    cds_records: Iterable[GffRecord],
    strand: str,
) -> str:
    """Extract + concatenate CDS sequences from the FASTA in genomic order,
    reverse-complementing the concatenated result if the transcript is on the
    minus strand.

    This is the Phase-0-fixed splicing pattern, factored out once.
    """
    parts: list[str] = []
    for cds in cds_records:
        if cds.seqid not in set(fasta.references):
            raise ValueError(
                f"CDS {cds.id!r} references seqid {cds.seqid!r} not in FASTA "
                f"(have: {sorted(fasta.references)})"
            )
        # pysam.fetch is 0-based half-open
        parts.append(fasta.fetch(cds.seqid, cds.start - 1, cds.end))
    full = "".join(parts)
    if strand == "-":
        full = reverse_complement(full)
    return full


@dataclass
class ProteinRecord:
    """One translated protein sequence + QC flags."""

    gene_id: str
    transcript_id: str
    seqid: str
    strand: str
    start: int                                       # 1-based inclusive on the FASTA
    end: int
    cds_length: int                                  # total spliced CDS length in bp
    protein: str                                     # translated sequence, stop codons shown as '*'
    start_codon_ok: bool                             # first codon is ATG
    stop_codon_ok: bool                              # last codon is a stop (before trailing '*' is stripped)
    premature_stop_at: int | None                    # 1-based aa index of first in-frame '*' BEFORE the last codon
    trailing_partial_codon: int = 0                  # 0, 1 or 2 leftover bp after the last full codon
    warnings: list[str] = field(default_factory=list)


def translate_cds(cds_dna: str, phase: int = 0) -> tuple[str, int]:
    """Translate a CDS DNA sequence starting at `phase` into protein.

    Returns (protein, trailing_partial_codon). The protein includes terminal
    `*` stop codons; callers decide whether to strip them.
    """
    seq = cds_dna[phase:]
    trailing = len(seq) % 3
    if trailing:
        seq = seq[:-trailing]
    return str(Seq(seq).translate()), trailing


def translate_transcript(
    fasta: pysam.FastaFile,
    gff_records: list[GffRecord],
    transcript: GffRecord,
    gene_id: str,
    logger: logging.Logger | None = None,
) -> ProteinRecord:
    """Build a ProteinRecord for one mRNA by slicing its CDS children."""
    logger = logger or logging.getLogger(__name__)
    cds = cds_of_transcript(gff_records, transcript.id)
    if not cds:
        return ProteinRecord(
            gene_id=gene_id,
            transcript_id=transcript.id,
            seqid=transcript.seqid,
            strand=transcript.strand,
            start=transcript.start,
            end=transcript.end,
            cds_length=0,
            protein="",
            start_codon_ok=False,
            stop_codon_ok=False,
            premature_stop_at=None,
            warnings=["no CDS records"],
        )

    dna = splice_cds(fasta, cds, transcript.strand)
    # Use the phase of the first CDS (on the + strand it's the first; on -
    # strand genomic-first is also transcript-first when we RC the whole).
    first_phase_str = cds[0].phase if transcript.strand == "+" else cds[-1].phase
    try:
        first_phase = int(first_phase_str) if first_phase_str != "." else 0
    except ValueError:
        first_phase = 0

    protein, trailing = translate_cds(dna, phase=first_phase)
    start_ok = protein.startswith("M")
    stop_ok = protein.endswith("*")
    # Premature stop: any '*' that is NOT the final character.
    premature_at: int | None = None
    body = protein[:-1] if stop_ok else protein
    for i, aa in enumerate(body, start=1):
        if aa == "*":
            premature_at = i
            break

    return ProteinRecord(
        gene_id=gene_id,
        transcript_id=transcript.id,
        seqid=transcript.seqid,
        strand=transcript.strand,
        start=transcript.start,
        end=transcript.end,
        cds_length=len(dna),
        protein=protein,
        start_codon_ok=start_ok,
        stop_codon_ok=stop_ok,
        premature_stop_at=premature_at,
        trailing_partial_codon=trailing,
    )


def extract_proteins_for_gene(
    fasta_path: Path | str,
    gff_path: Path | str,
    gene_id: str,
    logger: logging.Logger | None = None,
) -> dict[str, ProteinRecord]:
    """Return {transcript_id: ProteinRecord} for every mRNA under `gene_id`."""
    logger = logger or logging.getLogger(__name__)
    records = load_gff(Path(gff_path))
    transcripts = transcripts_of_gene(records, gene_id)
    if not transcripts:
        logger.info("No mRNA records for gene %s in %s", gene_id, gff_path)
        return {}

    out: dict[str, ProteinRecord] = {}
    with pysam.FastaFile(str(fasta_path)) as fa:
        for tx in transcripts:
            out[tx.id] = translate_transcript(fa, records, tx, gene_id, logger=logger)
    return out


def write_protein_fasta(
    proteins: Iterable[ProteinRecord],
    out_path: Path | str,
    header_prefix: str = "",
) -> Path:
    """Write a multi-record protein FASTA. `header_prefix` is prepended to the
    transcript ID (e.g. 'hap1.' or 'ref.') to keep records disambiguable when
    multiple haplotypes are concatenated."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        for p in proteins:
            if not p.protein:
                continue
            f.write(f">{header_prefix}{p.transcript_id} gene={p.gene_id} "
                    f"strand={p.strand} cds_length={p.cds_length} "
                    f"start_ok={p.start_codon_ok} stop_ok={p.stop_codon_ok}"
                    f"{f' pts={p.premature_stop_at}' if p.premature_stop_at else ''}\n")
            # Wrap at 60 characters, standard FASTA
            seq = p.protein
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")
    return out_path
