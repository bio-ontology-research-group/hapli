"""
Reference vs. haplotype protein diff.

Given a reference protein sequence and a haplotype protein sequence, produce a
structured `ProteinDiff` (see `hapli.core.schema.ProteinDiff`) carrying:

  * overall identity (matches / max(ref_len, hap_len))
  * the list of single-residue substitutions
  * position of the first premature stop, if any
  * a frameshift-region summary when an aligned block has a cluster of
    offsets consistent with a frame shift (detected heuristically from runs
    of insertions/deletions around the same window)

The alignment is done with BioPython's pairwise2 (already a transitive dep
via biopython>=1.86). For production-scale MAVE benchmarks we may swap in
parasail later; the function signature is stable.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterable

from Bio import Align

from .schema import ProteinDiff


# BLOSUM62-like scoring tuned for global alignment of closely-related
# sequences (haplotype vs reference). pairwise2 is deprecated; the modern
# Bio.Align.PairwiseAligner is the supported path.
def _make_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    return aligner


def _walk_alignment(ref_aln: str, hap_aln: str) -> tuple[
    list[dict],                 # substitutions
    list[tuple[int, int]],      # insertion runs in hap (hap_pos_1based, length)
    list[tuple[int, int]],      # deletion runs in hap (ref_pos_1based, length)
    int,                        # matches
]:
    """Walk a pairwise alignment and classify positions.

    Parameters
    ----------
    ref_aln, hap_aln : strings of equal length, with '-' for gaps.

    Returns
    -------
    substitutions : list of {ref_pos, ref_aa, hap_pos, hap_aa} dicts
    insertions    : runs of consecutive '-' in ref_aln (i.e. extra aa on hap)
    deletions     : runs of consecutive '-' in hap_aln (i.e. missing aa on hap)
    matches       : count of exact-match residue columns
    """
    substitutions: list[dict] = []
    insertions: list[tuple[int, int]] = []
    deletions: list[tuple[int, int]] = []
    matches = 0

    ref_pos = 0   # 0-based running index into the ungapped reference
    hap_pos = 0   # 0-based running index into the ungapped haplotype

    ins_start: int | None = None
    del_start: int | None = None
    ins_run_hap_start: int | None = None
    del_run_ref_start: int | None = None
    ins_run: int = 0
    del_run: int = 0

    for rc, hc in zip(ref_aln, hap_aln):
        if rc == "-" and hc != "-":
            if ins_start is None:
                ins_start = hap_pos + 1
            ins_run += 1
            hap_pos += 1
        else:
            if ins_start is not None:
                insertions.append((ins_start, ins_run))
                ins_start = None
                ins_run = 0

        if hc == "-" and rc != "-":
            if del_start is None:
                del_start = ref_pos + 1
            del_run += 1
            ref_pos += 1
        else:
            if del_start is not None:
                deletions.append((del_start, del_run))
                del_start = None
                del_run = 0

        if rc != "-" and hc != "-":
            ref_pos += 1
            hap_pos += 1
            if rc == hc:
                matches += 1
            else:
                substitutions.append({
                    "ref_pos": ref_pos,
                    "ref_aa": rc,
                    "hap_pos": hap_pos,
                    "hap_aa": hc,
                })

    # Flush tail runs
    if ins_start is not None:
        insertions.append((ins_start, ins_run))
    if del_start is not None:
        deletions.append((del_start, del_run))

    return substitutions, insertions, deletions, matches


def _detect_frameshift_region(
    substitutions: list[dict],
    insertions: list[tuple[int, int]],
    deletions: list[tuple[int, int]],
    ref_length: int,
    hap_length: int,
) -> dict | None:
    """Heuristically detect a frameshift-window pattern.

    Two signatures we flag:

      * **Rescued** frameshift — a dense cluster of substitutions with no net
        length change at the protein level (total ref_length == hap_length, or
        flanking insertion + deletion that cancel at the nucleotide level).
        This is the "+1 insertion / -1 deletion on the same haplotype" pattern:
        the protein's reading frame is shifted between the two indels, then
        restored. At the protein-alignment level this manifests as a short
        cluster of substitutions with no surrounding gaps.

      * **Unrescued** frameshift — the alignment carries a terminal indel run
        beyond a substitution cluster, consistent with a frame-disrupting event
        that runs to the end of the protein.

    We don't try to be definitive; callers should treat `frameshift_region` as
    a flag to surface to humans / downstream function predictors.
    """
    if not substitutions:
        return None
    # Cluster consecutive substitutions within 3 aa of each other.
    clusters: list[list[dict]] = []
    cur: list[dict] = [substitutions[0]]
    for s in substitutions[1:]:
        if s["ref_pos"] - cur[-1]["ref_pos"] <= 3:
            cur.append(s)
        else:
            clusters.append(cur)
            cur = [s]
    clusters.append(cur)

    big = max(clusters, key=len)
    if len(big) < 3:
        return None

    first_pos = big[0]["ref_pos"]
    last_pos = big[-1]["ref_pos"]
    flanking_ins = [ix for ix in insertions if first_pos - 10 <= ix[0] <= last_pos + 10]
    flanking_del = [dx for dx in deletions if first_pos - 10 <= dx[0] <= last_pos + 10]

    # Rescued: protein lengths match exactly, so frame must have been restored
    # (or flanking indels cancel at the amino-acid level).
    length_rescued = ref_length == hap_length
    indel_rescued = (
        sum(l for _, l in flanking_ins) == sum(l for _, l in flanking_del)
        and bool(flanking_ins) and bool(flanking_del)
    )
    restored = length_rescued or indel_rescued

    return {
        "start": first_pos,
        "end": last_pos,
        "n_substitutions": len(big),
        "flanking_insertions": flanking_ins,
        "flanking_deletions": flanking_del,
        "restored": restored,
    }


def diff_proteins(
    transcript_id: str,
    haplotype: int,
    ref_seq: str,
    hap_seq: str,
    strip_trailing_stop: bool = True,
    logger: logging.Logger | None = None,
) -> ProteinDiff:
    """Align ref against hap, return a structured ProteinDiff."""
    logger = logger or logging.getLogger(__name__)
    # Detect the first in-frame stop in the haplotype BEFORE any stripping.
    # A stop is "premature" when:
    #   (a) residues follow it in the haplotype (mid-sequence '*'), OR
    #   (b) it's the last residue but the haplotype is materially shorter than
    #       the reference, meaning the protein was truncated.
    raw_hap = hap_seq
    raw_ref_len = len(ref_seq.rstrip("*"))
    premature_at: int | None = None
    for i, aa in enumerate(raw_hap, start=1):
        if aa == "*":
            followed_by_residues = i < len(raw_hap)
            truncated = (not followed_by_residues) and (i - 1) < (raw_ref_len - 5)
            if followed_by_residues or truncated:
                premature_at = i
            break

    if strip_trailing_stop:
        ref_seq = ref_seq.rstrip("*")
        hap_seq = hap_seq.rstrip("*")

    if not ref_seq and not hap_seq:
        return ProteinDiff(
            transcript=transcript_id,
            haplotype=haplotype,
            ref_length=0,
            hap_length=0,
            identity=1.0,
        )
    if not hap_seq:
        return ProteinDiff(
            transcript=transcript_id,
            haplotype=haplotype,
            ref_length=len(ref_seq),
            hap_length=0,
            identity=0.0,
        )

    aligner = _make_aligner()
    alignments = aligner.align(ref_seq, hap_seq)
    best = alignments[0]
    # Bio.Align.Alignment supports indexing: row 0 is the target (reference),
    # row 1 is the query (haplotype). Each row is the aligned sequence with '-'
    # placed at gap columns.
    ref_aln = best[0]
    hap_aln = best[1]

    substitutions, insertions, deletions, matches = _walk_alignment(ref_aln, hap_aln)
    denom = max(len(ref_seq), len(hap_seq))
    identity = matches / denom if denom > 0 else 1.0

    frameshift_region = _detect_frameshift_region(
        substitutions, insertions, deletions,
        ref_length=len(ref_seq), hap_length=len(hap_seq),
    )

    return ProteinDiff(
        transcript=transcript_id,
        haplotype=haplotype,
        ref_length=len(ref_seq),
        hap_length=len(hap_seq),
        identity=identity,
        substitutions=substitutions,
        premature_stop_at=premature_at,
        frameshift_region=frameshift_region,
    )


def diff_many(
    transcript_id: str,
    ref_seq: str,
    hap_seqs: dict[int, str],
    **kwargs,
) -> list[ProteinDiff]:
    """Diff the reference against each haplotype. `hap_seqs` maps 1/2 → seq."""
    return [
        diff_proteins(transcript_id, h, ref_seq, s, **kwargs)
        for h, s in hap_seqs.items()
    ]
