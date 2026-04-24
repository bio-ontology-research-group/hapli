#!/usr/bin/env python3
"""
Analyse a hapli-scored double-mutant MAVE scoreset.

Compares three predictors of the experimental MAVE score:

  * S_additive — ESM per-variant log-odds summed in reference context.
                 This is the baseline a variant-by-variant tool (LOFTEE,
                 BCFtools/csq, raw AlphaMissense-aggregate) would produce.

  * S_joint    — ESM full-sequence pseudo-log-likelihood delta. Contains
                 the joint / haplotype context.

  * residual   — S_joint − S_additive. The hapli research contribution;
                 captures the part of the fitness signal that variant-by-
                 variant cannot reach.

Reports Spearman / Pearson correlations of each predictor against the MAVE
score, and their AUROC at a fitness-threshold-derived LoF label. Also
reports the delta between S_joint and S_additive for the same rows — the
headline "epistasis gain" number for the paper.
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
    ap.add_argument("--category", default="missense_multi",
                    help="restrict analysis to rows of this category (default: missense_multi)")
    ap.add_argument("--lof-threshold", type=float, default=0.1,
                    help="MAVE scores on the LoF side count as LoF (default: 0.1; GB1 fitness scale)")
    ap.add_argument("--lof-is", choices=["low", "high"], default="low",
                    help="whether LOW or HIGH MAVE scores indicate LoF (default: low)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--figure", type=Path, default=None)
    args = ap.parse_args(argv)

    rows: list[dict[str, float]] = []
    with args.scored.open() as f:
        for r in csv.DictReader(f):
            if r.get("category") != args.category:
                continue
            try:
                ms = float(r["mave_score"])
                sa = float(r["s_additive"])
                sj = float(r["s_joint"])
                res = float(r["residual"])
                hi = float(r["hapli_score_heuristic"])
            except (ValueError, KeyError, TypeError):
                continue
            rows.append({
                "mave_score": ms,
                "s_additive": sa,
                "s_joint": sj,
                "residual": res,
                "hapli_heuristic": hi,
            })

    if len(rows) < 50:
        print(f"only {len(rows)} rows pass filters; need more for statistics", file=sys.stderr)
        return 1

    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr, spearmanr
    from sklearn.metrics import roc_auc_score

    ms = np.array([r["mave_score"] for r in rows])
    sa = np.array([r["s_additive"] for r in rows])
    sj = np.array([r["s_joint"] for r in rows])
    res = np.array([r["residual"] for r in rows])
    hi = np.array([r["hapli_heuristic"] for r in rows])

    # Binary LoF label from the MAVE score
    if args.lof_is == "low":
        lof = (ms < args.lof_threshold).astype(int)
    else:
        lof = (ms > args.lof_threshold).astype(int)

    def _corr(name, pred, inverted=True):
        y = -pred if inverted else pred
        rho_s, p_s = spearmanr(y, ms)
        rho_p, p_p = pearsonr(y, ms)
        return {
            "predictor": name,
            "spearman": float(rho_s), "spearman_p": float(p_s),
            "pearson": float(rho_p),  "pearson_p": float(p_p),
        }

    def _auroc(name, pred):
        if lof.sum() == 0 or lof.sum() == len(lof):
            return None
        return roc_auc_score(lof, -pred)   # more negative pred → more LoF

    summary = {
        "scored_file": str(args.scored),
        "category": args.category,
        "n_rows": len(rows),
        "n_lof": int(lof.sum()),
        "n_functional": int((1 - lof).sum()),
        "lof_threshold": args.lof_threshold,
        "lof_is": args.lof_is,
        "correlations": [
            _corr("S_additive", sa),
            _corr("S_joint", sj),
            _corr("residual", res),
            _corr("hapli_heuristic", hi, inverted=False),
        ],
        "auroc": {
            "S_additive": _auroc("S_additive", sa),
            "S_joint":    _auroc("S_joint", sj),
            "residual":   _auroc("residual", res),
            # |residual| alone doesn't distinguish rescue (+) from synergy (−);
            # its AUROC against LoF is reported for completeness but not claimed.
            "abs_residual": _auroc("abs_residual", -np.abs(res)) if len(rows) > 0 else None,
        },
        # Headline: does S_joint correlate better with MAVE than S_additive does?
        # If yes, the haplotype-level view adds signal over the variant-by-variant baseline.
        "joint_vs_additive_delta": {
            "spearman_gain":  float(spearmanr(-sj, ms)[0] - spearmanr(-sa, ms)[0]),
            "mean_abs_residual": float(np.abs(res).mean()),
        },
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    # Four-panel figure
    fig_path = args.figure or args.out.with_suffix(".png")
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes[0, 0].scatter(-sa, ms, s=4, alpha=0.4)
    axes[0, 0].set_xlabel("-S_additive (high = model says deleterious)")
    axes[0, 0].set_ylabel("MAVE experimental score")
    axes[0, 0].set_title(f"S_additive  Spearman={spearmanr(-sa, ms)[0]:.3f}")
    axes[0, 1].scatter(-sj, ms, s=4, alpha=0.4, color="tab:orange")
    axes[0, 1].set_xlabel("-S_joint")
    axes[0, 1].set_ylabel("MAVE experimental score")
    axes[0, 1].set_title(f"S_joint     Spearman={spearmanr(-sj, ms)[0]:.3f}")
    axes[1, 0].scatter(res, ms, s=4, alpha=0.4, color="tab:green")
    axes[1, 0].set_xlabel("residual = S_joint - S_additive")
    axes[1, 0].set_ylabel("MAVE experimental score")
    axes[1, 0].set_title(f"residual    Spearman={spearmanr(res, ms)[0]:.3f}")
    axes[1, 1].hist(res, bins=40, color="tab:green", alpha=0.7)
    axes[1, 1].set_xlabel("residual"); axes[1, 1].set_ylabel("count")
    axes[1, 1].set_title(f"Residual distribution\n(mean |res|={np.abs(res).mean():.2f})")
    fig.suptitle(f"{args.scored.stem} — {len(rows)} {args.category}  (hapli vs. MAVE)")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=120)
    print(f"wrote figure {fig_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
