#!/usr/bin/env python3
"""
Analyse a hapli-scored MAVE scoreset: compute AUROC of ESM2 log-odds as a
predictor of experimental LoF, and emit a PR / ROC figure.

This is the single-mutant baseline for the Phase 5 MAVE benchmark. The
interpretation is: ESM2's per-variant log-odds (S_additive per variant) is
the baseline that the epistasis residual is measured *against* on the
double-mutant leg. Getting respectable single-mutant AUROC confirms the
integration works; the headline number for the paper is the double-mutant
residual vs. baseline delta (next benchmark sub-task).

Usage:
  ./analyze.py --scored scored.csv --lof-threshold 0.5 --out-dir results/
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scored", type=Path, required=True)
    ap.add_argument("--lof-threshold", type=float, default=0.5,
                    help="MAVE scores on the LoF side of this threshold are labelled LoF (default: 0.5)")
    ap.add_argument("--lof-is", choices=["low", "high"], default="low",
                    help="whether LOW or HIGH MAVE scores indicate LoF. "
                         "'low' (default) is the Findlay 2018 SGE / Giacomelli TP53 convention; "
                         "'high' is the Starita 2015 BRCA1 RING depletion convention.")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--figure", type=Path, default=None,
                    help="optional PNG path for the ROC+PR figure")
    args = ap.parse_args(argv)

    log_odds: list[float] = []
    mave: list[float] = []
    with args.scored.open() as f:
        for r in csv.DictReader(f):
            try:
                log_odds.append(float(r["esm_log_odds"]))
                mave.append(float(r["mave_score"]))
            except (KeyError, ValueError):
                continue

    if not log_odds:
        print("no rows to analyse", file=sys.stderr)
        return 1

    # Binary label: which side of the threshold is LoF depends on the scoreset.
    if args.lof_is == "low":
        labels = [1 if m < args.lof_threshold else 0 for m in mave]
    else:
        labels = [1 if m > args.lof_threshold else 0 for m in mave]
    # ESM2 predictor: more negative log_odds → more deleterious → higher
    # "LoF probability" after sign flip.
    pred = [-x for x in log_odds]

    import numpy as np
    from sklearn.metrics import (
        roc_auc_score, average_precision_score,
        roc_curve, precision_recall_curve,
    )
    y = np.array(labels)
    s = np.array(pred)
    if y.sum() == 0 or y.sum() == len(y):
        print(f"degenerate labels (all-{'LoF' if y.sum() else 'functional'}); "
              f"can't compute AUROC", file=sys.stderr)
        return 1

    roc_auc = roc_auc_score(y, s)
    pr_auc = average_precision_score(y, s)
    fpr, tpr, _ = roc_curve(y, s)
    pr_p, pr_r, _ = precision_recall_curve(y, s)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "scored_file": str(args.scored),
        "n_variants": int(len(y)),
        "n_lof": int(y.sum()),
        "n_functional": int((1 - y).sum()),
        "lof_threshold": args.lof_threshold,
        "lof_is": args.lof_is,
        "roc_auc": float(roc_auc),
        "pr_auc": float(pr_auc),
        "predictor": "esm_log_odds (inverted)",
    }
    summary_path = args.out_dir / (args.scored.stem + ".auroc.json")
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    fig_path = args.figure or args.out_dir / (args.scored.stem + ".auroc.png")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(fpr, tpr, label=f"ESM2 log-odds (AUROC={roc_auc:.3f})")
    ax1.plot([0, 1], [0, 1], "--", color="grey")
    ax1.set_xlabel("FPR"); ax1.set_ylabel("TPR"); ax1.set_title("ROC")
    ax1.legend()
    ax2.plot(pr_r, pr_p, label=f"ESM2 log-odds (AUPR={pr_auc:.3f})")
    ax2.set_xlabel("Recall"); ax2.set_ylabel("Precision"); ax2.set_title("PR")
    ax2.legend()
    fig.suptitle(f"{args.scored.stem} — {len(y)} variants, {int(y.sum())} LoF")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=120)
    print(f"wrote figure {fig_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
