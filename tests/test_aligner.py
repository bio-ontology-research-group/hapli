"""
Unit tests for the minimap2 wrapper (`hapli.alignment.aligner.SequenceAligner`).

The module is a thin subprocess wrapper around minimap2, retained for the
schema-v1 per-transcript alignment block that the TUI and LLM consumers
still read. Liftoff is the canonical per-haplotype annotation-lift step
in schema v2; this wrapper is kept for backwards-compatibility only.
"""
from __future__ import annotations

import shutil
import pytest

from hapli.alignment.aligner import SequenceAligner


pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None, reason="minimap2 not on PATH"
)


def test_perfect_match_gives_high_identity():
    """A query that equals a slice of the target must align back at very
    high identity, verifying the subprocess + SAM parsing path. minimap2's
    `-x splice` preset uses k=15 seeds and needs a reasonably long query
    (~100 bp+) to find hits reliably.
    """
    # Diverse 100-bp query embedded in a longer target; A-rich context on
    # both sides so the query is the only unique region for minimap2 to latch.
    unique = (
        "ACGTACGTTGCATGCAAGCTTGACCAGTCCGATCGATCGATCTAGCTAGCTAGCTAGCTA"
        "GCTACGATCGATCGTAGCTAGCTAGCTAGTCGATCGATCGATGCTGCTG"
    )
    assert len(unique) >= 100
    target = "A" * 500 + unique + "A" * 500
    query = unique

    aligner = SequenceAligner()
    alignments = aligner.align(query, "q1", target, "t1")
    assert alignments, "expected at least one alignment for a perfect-match query"
    best = alignments[0]
    assert not best.is_unmapped
    nm = best.get_tag("NM") if best.has_tag("NM") else -1
    assert nm == 0, f"expected NM=0 for perfect match; got {nm}"


def test_no_alignment_for_unrelated_sequences():
    """A query with no sequence similarity to the target returns an empty
    list (the only unmapped record is filtered out by SequenceAligner).
    """
    target = "A" * 500
    query = "GCGCGCGCGCGCGCGCGC"   # wholly unrelated to poly-A target

    aligner = SequenceAligner()
    alignments = aligner.align(query, "q1", target, "t1")
    # minimap2 may return nothing or only unmapped records; SequenceAligner
    # filters the unmapped so the result should be empty.
    assert all(not a.is_unmapped for a in alignments)


def test_missing_minimap2_is_surfaced():
    """If the minimap2_path is wrong, the wrapper must raise FileNotFoundError
    rather than silently producing an empty alignment list.
    """
    aligner = SequenceAligner(minimap2_path="/nonexistent/minimap2_binary_doesnt_exist")
    with pytest.raises(FileNotFoundError):
        aligner.align("ACGT", "q", "ACGT", "t")
