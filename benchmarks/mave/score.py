#!/usr/bin/env python3
"""
Score a MaveDB scoreset with hapli's ESM2-based per-variant log-odds.

This is the single-mutant leg of the Phase 5 MAVE benchmark. For each row in
the scoreset we parse the HGVS protein notation, validate it against the
reference sequence, and compute:

  * esm_log_odds   — log p(alt | ref_context) - log p(ref | ref_context), using
                     a single ESM2 forward pass on the reference protein
  * esm_ref_lp     — log p(ref | ref_context), the model's prior at that site
  * esm_alt_lp     — log p(alt | ref_context)

Output CSV columns:

  accession, hgvs_pro, position, ref_aa, alt_aa, mave_score, n_replicates,
  esm_ref_lp, esm_alt_lp, esm_log_odds

One forward pass on the reference protein (~30 s on CPU for 1.8k aa with
the 8M-param model) serves all variants — the heavy cost is amortised.
Double-mutant benchmarks require a separate script (one forward pass per
unique haplotype protein) — documented in the Phase 5 README.

Usage:
  ./score.py --scores urn_mavedb_...scores.csv
             --reference urn_mavedb_...reference.fa
             --out scored.csv
             [--model esm2_t6_8M_UR50D]
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import time
from pathlib import Path

import numpy as np


# three-letter → one-letter amino-acid code (HGVS protein uses three-letter)
_AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}
_HGVS_PRO_SUB = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|=|Ter)$")


def parse_hgvs_pro(hgvs: str) -> tuple[str, int, str] | None:
    """Parse a simple HGVS protein substitution. Returns (ref_aa, pos, alt_aa)
    or None if the notation is something we don't score (synonymous, indel,
    frameshift, nonsense with a stop, multi-variant). Uses 1-based `pos`."""
    if not hgvs or hgvs == "NA":
        return None
    m = _HGVS_PRO_SUB.match(hgvs.strip())
    if not m:
        return None
    ref_aa3, pos_str, alt_aa3 = m.groups()
    ref_aa = _AA3.get(ref_aa3)
    if ref_aa is None:
        return None
    if alt_aa3 == "=":        # synonymous - skip for missense benchmark
        return None
    alt_aa = _AA3.get(alt_aa3) if alt_aa3 in _AA3 else alt_aa3
    if not alt_aa:
        return None
    # Skip nonsense / stop-gain (Ter, *) - ESM2 has no amino-acid token for them.
    if alt_aa in ("*", "Ter"):
        return None
    return ref_aa, int(pos_str), alt_aa


def load_reference_fasta(path: Path) -> str:
    seq: list[str] = []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--model", default="esm2_t6_8M_UR50D")
    ap.add_argument("--max-rows", type=int, default=None, help="cap rows for quick runs")
    args = ap.parse_args(argv)

    ref_seq = load_reference_fasta(args.reference)
    print(f"reference: {len(ref_seq)} aa", file=sys.stderr)

    # Lazy import so other scripts that just want the parser can run without torch.
    from hapli.interpretation.esm_scoring import ESM2Scorer
    scorer = ESM2Scorer(model_name=args.model)

    # One forward pass on the reference → log-prob matrix over all positions × tokens
    # (cached via ESM2Scorer's mechanism — downstream alt_log_probs calls are O(n)).
    t0 = time.time()
    per_pos = scorer.per_position_log_probs(ref_seq)
    print(f"ref pLL per-position computed in {time.time() - t0:.1f}s "
          f"(cached by sequence sha256)", file=sys.stderr)

    # Gather rows + parsed HGVS
    rows: list[dict[str, str]] = []
    skipped = 0
    with args.scores.open() as f:
        reader = csv.DictReader(f)
        for r in reader:
            if args.max_rows is not None and len(rows) >= args.max_rows:
                break
            hgvs = r.get("hgvs_pro") or ""
            parsed = parse_hgvs_pro(hgvs)
            if parsed is None:
                skipped += 1
                continue
            ref_aa, pos, alt_aa = parsed
            if pos < 1 or pos > len(ref_seq):
                skipped += 1; continue
            if ref_seq[pos - 1] != ref_aa:
                skipped += 1; continue      # HGVS doesn't match reference
            try:
                mave_score = float(r["score"])
            except (ValueError, KeyError, TypeError):
                skipped += 1; continue
            rows.append({
                "accession": r.get("accession", ""),
                "hgvs_pro": hgvs,
                "position": pos,
                "ref_aa": ref_aa,
                "alt_aa": alt_aa,
                "mave_score": mave_score,
                "n_replicates": r.get("replicates", ""),
            })
    print(f"kept {len(rows)} variants, skipped {skipped}", file=sys.stderr)

    # Batch: one alt_log_probs call per unique alt in its position — but simpler
    # to just call alt_log_probs on the whole list (one forward pass already done).
    positions = [r["position"] for r in rows]
    alts = [r["alt_aa"] for r in rows]
    # alt_log_probs does one forward on the reference (cached) and then per-position
    # lookups into the full logit matrix.
    alt_lps = scorer.alt_log_probs(ref_seq, positions, alts)
    ref_lps = per_pos[[p - 1 for p in positions]]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "accession", "hgvs_pro", "position", "ref_aa", "alt_aa",
            "mave_score", "n_replicates",
            "esm_ref_lp", "esm_alt_lp", "esm_log_odds",
        ])
        for r, rlp, alp in zip(rows, ref_lps, alt_lps):
            w.writerow([
                r["accession"], r["hgvs_pro"], r["position"], r["ref_aa"], r["alt_aa"],
                f"{r['mave_score']:.6g}", r["n_replicates"],
                f"{float(rlp):.6f}", f"{float(alp):.6f}",
                f"{float(alp - rlp):.6f}",
            ])
    print(f"wrote {args.out} ({len(rows)} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
