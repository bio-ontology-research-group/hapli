#!/usr/bin/env python3
"""
Position-stratified epistasis analysis on a pairwise-scan MAVE dataset.

Epistasis theory predicts that residue pairs close in sequence (or in 3D
space) carry more epistatic signal than distant pairs — contact-mediated
interactions concentrate near each other. A protein language model that
has learned real co-evolutionary structure should reflect this: its
haplotype-level residual should be strongest on close pairs.

This script bins the hapli-scored double-mutants by `|pos1 − pos2|` and
reports, per bin:

  * n_rows
  * mean |residual|
  * Spearman(residual, mave_score)
  * Spearman(residual, E_exp) if the raw scores CSV is provided

Ideal result: closer bins have larger |residual| and stronger Spearman
correlations with E_exp, tapering off for distant bins.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path


_RE_DOUBLE_BRACKET = re.compile(
    r"^p\.\["
    r"(?:[A-Z][a-z]{2}|[A-Z])(\d+)(?:[A-Z][a-z]{2}|[A-Z])"
    r";"
    r"(?:[A-Z][a-z]{2}|[A-Z])(\d+)(?:[A-Z][a-z]{2}|[A-Z])"
    r"\]$"
)
_RE_ADJ_DELINS = re.compile(
    r"^p\.(?:[A-Z][a-z]{2}|[A-Z])(\d+)_(?:[A-Z][a-z]{2}|[A-Z])(\d+)delins"
)


def position_pair(hgvs: str) -> tuple[int, int] | None:
    m = _RE_DOUBLE_BRACKET.match(hgvs)
    if m:
        p1, p2 = map(int, m.groups())
        return (min(p1, p2), max(p1, p2))
    m = _RE_ADJ_DELINS.match(hgvs)
    if m:
        p1, p2 = map(int, m.groups())
        return (min(p1, p2), max(p1, p2))
    return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scored", type=Path, required=True,
                    help="hapli-scored CSV with missense_multi rows")
    ap.add_argument("--bins", default="1,5,10,30,200",
                    help="comma-separated upper-distance thresholds (default: 1,5,10,30,200)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--figure", type=Path, default=None)
    args = ap.parse_args(argv)

    bin_edges = [int(x) for x in args.bins.split(",")]
    bin_labels = []
    prev = 0
    for b in bin_edges:
        bin_labels.append(f"{prev + 1}-{b}")
        prev = b

    # Collect per-bin
    bucket: dict[str, list[dict]] = defaultdict(list)
    distances: list[int] = []

    with args.scored.open() as f:
        for r in csv.DictReader(f):
            if r.get("category") != "missense_multi":
                continue
            pair = position_pair(r["hgvs_pro"])
            if pair is None:
                continue
            d = pair[1] - pair[0]
            # pick bin
            label = None
            prev = 0
            for i, edge in enumerate(bin_edges):
                if d <= edge:
                    label = bin_labels[i]; break
                prev = edge
            if label is None:
                label = f">{bin_edges[-1]}"
            try:
                residual = float(r["residual"])
                mave = float(r["mave_score"])
                s_add = float(r["s_additive"])
                s_joint = float(r["s_joint"])
            except (ValueError, KeyError, TypeError):
                continue
            bucket[label].append({
                "pair": pair, "distance": d,
                "residual": residual, "mave_score": mave,
                "s_additive": s_add, "s_joint": s_joint,
            })
            distances.append(d)

    if not bucket:
        print("no rows collected", file=sys.stderr)
        return 1

    import numpy as np
    from scipy.stats import spearmanr

    stats = {}
    for label in bin_labels + [f">{bin_edges[-1]}"]:
        rows = bucket.get(label, [])
        if not rows:
            continue
        res = np.array([r["residual"] for r in rows])
        mave = np.array([r["mave_score"] for r in rows])
        s_add = np.array([r["s_additive"] for r in rows])
        s_joint = np.array([r["s_joint"] for r in rows])
        stats[label] = {
            "n": len(rows),
            "mean_abs_residual": float(np.abs(res).mean()),
            "std_residual": float(res.std()),
            # Correlation of each predictor with observed fitness, within this bin.
            "spearman_residual_vs_mave":    float(spearmanr(res, mave)[0]) if len(res) >= 20 else None,
            "spearman_s_joint_vs_mave":     float(spearmanr(s_joint, mave)[0]) if len(res) >= 20 else None,
            "spearman_s_additive_vs_mave":  float(spearmanr(s_add, mave)[0]) if len(res) >= 20 else None,
        }

    summary = {
        "scored_file": str(args.scored),
        "total_rows": sum(len(v) for v in bucket.values()),
        "distance_bins": bin_labels + [f">{bin_edges[-1]}"],
        "per_bin": stats,
        "distance_range": {"min": int(min(distances)), "max": int(max(distances))},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    # Figure: mean |residual| by bin + correlation strength by bin
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = [l for l in bin_labels + [f">{bin_edges[-1]}"] if l in stats]
    x = range(len(labels))
    mean_abs = [stats[l]["mean_abs_residual"] for l in labels]
    corr_add = [abs(stats[l]["spearman_s_additive_vs_mave"] or 0) for l in labels]
    corr_joint = [abs(stats[l]["spearman_s_joint_vs_mave"] or 0) for l in labels]
    ns = [stats[l]["n"] for l in labels]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.bar(x, mean_abs, color="tab:green", alpha=0.7)
    ax1.set_xticks(x); ax1.set_xticklabels(labels)
    ax1.set_xlabel("|pos1 - pos2| bin")
    ax1.set_ylabel("mean |residual|")
    ax1.set_title("Residual magnitude by position-pair distance")
    for xi, n in zip(x, ns):
        ax1.text(xi, 0, f"n={n}", ha="center", va="bottom", fontsize=8, rotation=0)

    ax2.plot(x, corr_add, "o-", label="|Spearman(S_additive, MAVE)|", color="tab:red")
    ax2.plot(x, corr_joint, "o-", label="|Spearman(S_joint, MAVE)|", color="tab:blue")
    ax2.set_xticks(x); ax2.set_xticklabels(labels)
    ax2.set_xlabel("|pos1 - pos2| bin")
    ax2.set_ylabel("|Spearman correlation with MAVE|")
    ax2.set_title("Predictor strength by position-pair distance")
    ax2.legend()
    fig.tight_layout()
    fig_path = args.figure or args.out.with_suffix(".png")
    fig.savefig(fig_path, dpi=120)
    print(f"wrote figure {fig_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
