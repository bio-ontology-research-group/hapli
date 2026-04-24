#!/usr/bin/env python3
"""
Render the head-to-head JSON as a publication-quality figure.

Two panels:
  1. Per-case categorical summary: bcftools/csq view vs hapli view, with
     hapli-only-signal annotations.
  2. Tally: count of "hapli-only" signals across the case suite, broken
     out by signal type.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", type=Path,
                    default=Path("benchmarks/gnomad_flips/results/head_to_head.json"))
    ap.add_argument("--out", type=Path,
                    default=Path("benchmarks/gnomad_flips/results/head_to_head.png"))
    args = ap.parse_args(argv)

    if not args.summary.exists():
        print(f"summary missing: {args.summary} — run run_head_to_head.py first", file=sys.stderr)
        return 1

    summary = json.loads(args.summary.read_text())
    rows = [r for r in summary["rows"] if "error" not in r]
    if not rows:
        print("no rows", file=sys.stderr)
        return 1

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5),
                                    gridspec_kw={"width_ratios": [3, 1]})

    # -- Panel 1: per-case summary --
    case_names = [r["case"] for r in rows]
    n = len(case_names)

    # csq-side: list LoF tags + missense/synonymous + total record count
    csq_lof_class = []
    csq_other = []
    for r in rows:
        c = r["csq"]
        any_lof = any(h.get("lof_class") for h in c.values())
        any_records = sum(h.get("n_records", 0) for h in c.values())
        csq_lof_class.append(any_lof)
        csq_other.append(any_records - sum(1 for h in c.values() if h.get("lof_class")))

    # hapli-side: hap1_score, hap2_score, compound_het_lof, epistasis_flagged
    h1 = [r["hapli"].get("hap1_score") for r in rows]
    h2 = [r["hapli"].get("hap2_score") for r in rows]
    chl = [r["hapli"].get("compound_het_lof") for r in rows]
    epi = [r["hapli"].get("epistasis_flagged") for r in rows]

    y = np.arange(n)
    bar_h = 0.35
    h1_vals = [v if v is not None else 0 for v in h1]
    h2_vals = [v if v is not None else 0 for v in h2]
    ax1.barh(y - bar_h / 2, h1_vals, bar_h, label="hap1 score", color="tab:blue", alpha=0.7)
    ax1.barh(y + bar_h / 2, h2_vals, bar_h, label="hap2 score", color="tab:cyan", alpha=0.7)

    # Annotations
    for i, r in enumerate(rows):
        signals = r["hapli_only_signals"] or []
        ann = []
        if any(s.startswith("compound_het_lof") for s in signals):
            ann.append("CHL")
        if any(s.startswith("epistasis") for s in signals):
            ann.append("EPI")
        if any("deleted" in s for s in signals):
            ann.append("DEL")
        if csq_lof_class[i]:
            ann.append("csq.LoF")
        text = " | ".join(ann) if ann else "(no flips)"
        ax1.text(1.05, i, text, va="center", fontsize=10, family="monospace")

    ax1.set_yticks(y)
    ax1.set_yticklabels(case_names)
    ax1.set_xlabel("hapli per-haplotype function score (1.0 = intact, 0.0 = LoF)")
    ax1.set_xlim(0, 2.0)
    ax1.axvline(0.5, color="grey", linestyle="--", linewidth=0.8)
    ax1.text(0.5, -0.3, "LoF threshold = 0.5", fontsize=8, color="grey")
    ax1.set_title("hapli vs. bcftools/csq — per-case head-to-head\n"
                  "(bars = hapli scores; right column = signal labels)",
                  fontsize=11)
    ax1.legend(loc="upper right")

    # -- Panel 2: tally --
    tally = summary["tally"]
    bars = [
        ("hapli\ncompound_het_lof", tally["hapli_compound_het_lof_calls"]),
        ("hapli\nepistasis flagged", tally["hapli_epistasis_flagged_calls"]),
        ("csq LoF\nhap-records", tally["csq_lof_haplotype_records"]),
    ]
    labels, counts = zip(*bars)
    bar_colors = ["tab:green", "tab:orange", "tab:red"]
    ax2.bar(range(len(labels)), counts, color=bar_colors, alpha=0.8)
    ax2.set_xticks(range(len(labels)))
    ax2.set_xticklabels(labels, fontsize=9)
    ax2.set_ylabel("count across case suite")
    ax2.set_title("Signal tally\n(top two are hapli-only)", fontsize=11)
    for i, c in enumerate(counts):
        ax2.text(i, c + 0.05, str(c), ha="center", fontsize=11, fontweight="bold")

    fig.suptitle(
        f"hapli vs. bcftools/csq head-to-head — {len(rows)} curated cases",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote figure {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
