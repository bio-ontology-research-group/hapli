#!/usr/bin/env python3
"""
Unified MaveDB scoreset scorer — handles single, double, N-way, nonsense,
frameshift, indel, and SV-shaped HGVS variants in one pass.

For each row of the scoreset:

  1. Parse the `hgvs_pro` column into one or more ProteinVariants (bracketed
     multi-variant form supported).
  2. Apply them to the UniProt reference protein → synthetic haplotype protein.
  3. Score the haplotype with the unified scorer (hapli.interpretation.haplotype_scoring):
       * category (missense_single / missense_multi / lof / indel / synonymous / ...)
       * hap_length, length_changed, premature_stop_at, identity
       * For missense-only equal-length haplotypes with ≥2 substitutions:
         S_additive, S_joint, residual (the Phase 3 epistasis signal).

  4. Join with the MAVE experimental score column.

Output CSV is wide enough to support every downstream benchmark:
  - single-mutant AUROC of ESM log-odds (vs. MAVE LoF labels) — use rows
    where `category = missense_single`.
  - double-mutant epistasis residual analysis — use rows where
    `category = missense_multi` and `residual is not null`.
  - LoF/SV-like accuracy — use rows where `category ∈ {lof, indel}` and
    check the `hapli_score_heuristic` column.

One ESM2 forward pass per unique haplotype protein sequence (cached by
sha256 → disk). The reference forward pass is also cached. On the GB1
pairwise scan this amortises to ~hours on a CPU for 500K rows with the
tiny model; use `--max-rows` for quick smoke tests.

Usage:
  ./score_haplotypes.py --scores scores.csv --reference reference.fa --out scored.csv
                        [--with-esm] [--model esm2_t6_8M_UR50D] [--max-rows N]
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

from hapli.core.hgvs import parse_hgvs_pro
from hapli.interpretation.haplotype_scoring import score_haplotype


def load_reference_fasta(path: Path) -> str:
    seq: list[str] = []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)


_OUTPUT_COLUMNS = [
    "accession", "hgvs_pro", "mave_score", "n_replicates",
    "category", "hap_length", "length_changed",
    "premature_stop_at", "identity",
    "n_missense", "n_nonsense", "n_frameshift", "n_indel",
    "n_synonymous", "n_skipped",
    "s_additive", "s_joint", "residual", "flagged",
    # A cheap shape-only LoF heuristic that doesn't need ESM — mirrors what the
    # hapli pipeline's diploid aggregator would emit per haplotype.
    "hapli_score_heuristic",
]


def _heuristic_score(identity: float, ptc: int | None, ref_len: int) -> float:
    """Simple analogue of the Phase 4 diploid score for use as a LoF heuristic.
    Length-changing or truncating haplotypes get low scores; missense-only
    haplotypes pass through `identity` verbatim."""
    if ptc is not None:
        frac = ptc / ref_len if ref_len else 0
        if frac <= 0.1:
            return 0.0
        return identity * frac
    return identity


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--with-esm", action="store_true",
                    help="compute the ESM residual where applicable (requires `ml` extra)")
    ap.add_argument("--model", default="esm2_t6_8M_UR50D")
    ap.add_argument("--max-rows", type=int, default=None,
                    help="cap parsed rows for quick runs")
    ap.add_argument("--resume", action="store_true",
                    help="if --out exists, skip rows whose hgvs_pro is already written and append")
    args = ap.parse_args(argv)

    ref_seq = load_reference_fasta(args.reference)
    print(f"reference: {len(ref_seq)} aa", file=sys.stderr)

    scorer = None
    if args.with_esm:
        from hapli.interpretation.esm_scoring import ESM2Scorer
        scorer = ESM2Scorer(model_name=args.model)
        # Prime the cache with a reference pass so lookups are O(1) thereafter.
        t0 = time.time()
        _ = scorer.per_position_log_probs(ref_seq)
        print(f"ref pLL computed in {time.time() - t0:.1f}s (cached)", file=sys.stderr)

    cat_counts: dict[str, int] = {}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    n_rows = 0; n_skipped = 0

    done: set[str] = set()
    mode = "w"
    if args.resume and args.out.exists() and args.out.stat().st_size > 0:
        with args.out.open() as fprev:
            prev = csv.DictReader(fprev)
            for row in prev:
                h = (row.get("hgvs_pro") or "").strip()
                if h:
                    done.add(h)
        print(f"resume: {len(done)} rows already scored in {args.out}", file=sys.stderr)
        mode = "a"

    t0 = time.time()
    with args.scores.open() as fin, args.out.open(mode, newline="") as fout:
        reader = csv.DictReader(fin)
        writer = csv.writer(fout)
        if mode == "w":
            writer.writerow(_OUTPUT_COLUMNS)
        for r in reader:
            if args.max_rows is not None and n_rows >= args.max_rows:
                break
            hgvs = (r.get("hgvs_pro") or "").strip()
            if not hgvs or hgvs == "NA":
                n_skipped += 1; continue
            if hgvs in done:
                continue
            variants = parse_hgvs_pro(hgvs)
            if not variants:
                n_skipped += 1; continue

            try:
                mave_score = float(r["score"])
            except (TypeError, ValueError, KeyError):
                mave_score_str = r.get("score", "")
            else:
                mave_score_str = f"{mave_score:.6g}"

            s = score_haplotype(ref_seq, hgvs, variants, scorer=scorer)
            cat_counts[s.category] = cat_counts.get(s.category, 0) + 1

            writer.writerow([
                r.get("accession", ""),
                hgvs,
                mave_score_str,
                r.get("replicates", r.get("err", "")),
                s.category,
                s.hap_length,
                int(s.length_changed),
                s.premature_stop_at if s.premature_stop_at is not None else "",
                f"{s.identity:.6f}",
                s.n_missense, s.n_nonsense, s.n_frameshift, s.n_indel,
                s.n_synonymous, s.n_skipped,
                f"{s.s_additive:.6f}" if s.s_additive is not None else "",
                f"{s.s_joint:.6f}"    if s.s_joint    is not None else "",
                f"{s.residual:.6f}"   if s.residual   is not None else "",
                int(s.flagged),
                f"{_heuristic_score(s.identity, s.premature_stop_at, len(ref_seq)):.6f}",
            ])
            n_rows += 1
            if n_rows % 5000 == 0:
                fout.flush()
                print(f"  {n_rows} rows in {time.time() - t0:.0f}s", file=sys.stderr)

    print(
        f"wrote {args.out}\n"
        f"  rows kept:    {n_rows}\n"
        f"  rows skipped: {n_skipped}\n"
        f"  category distribution: {cat_counts}\n"
        f"  wall time:    {time.time() - t0:.1f}s",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
