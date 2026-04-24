#!/usr/bin/env python3
"""
Compose the paper's headline figure from a completed score_haplotypes.py run.

Input: a `<scoreset>.scored.csv` produced by `score_haplotypes.py --with-esm`,
plus the matching raw `<urn>.scores.csv` for the experimental-epistasis
target (computed via marginal averaging à la Olson 2014).

Output: a 2×2 figure suitable for the paper's central claim:

  Top-left:    AUROC bars — S_additive vs S_joint vs residual
               (the per-LoF-threshold improvement of joint over additive)
  Top-right:   Scatter of hapli.residual vs experimental epistasis E_exp
               (the cleanest non-tautological epistasis correlation)
  Bottom-left: Scatter of hapli.S_joint vs MAVE fitness
  Bottom-right: residual distribution + |residual|>threshold count

Headline numbers are also printed to stdout so they can be pasted directly
into the paper's results section.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path


# Match either p.[X;Y] or p.X_Ydelins for double mutants in MaveDB
_RE_BRACKET = re.compile(
    r"^p\.\["
    r"(?:[A-Z][a-z]{2}|[A-Z])(\d+)(?:[A-Z][a-z]{2}|[A-Z])"
    r";"
    r"(?:[A-Z][a-z]{2}|[A-Z])(\d+)(?:[A-Z][a-z]{2}|[A-Z])"
    r"\]$"
)
_RE_DELINS = re.compile(
    r"^p\.(?:[A-Z][a-z]{2}|[A-Z])(\d+)_(?:[A-Z][a-z]{2}|[A-Z])(\d+)delins"
)
_AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
        "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
        "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}


def _aa(code: str) -> str | None:
    if code in _AA3: return _AA3[code]
    if len(code) == 1 and code in "ACDEFGHIKLMNPQRSTVWY": return code
    return None


def _parse_pair(hgvs: str) -> tuple[tuple[int, str], tuple[int, str]] | None:
    """Return ((pos1, alt_aa1), (pos2, alt_aa2)) or None.
    Skips synonymous positions (where ref == alt)."""
    m = _RE_BRACKET.match(hgvs)
    if m:
        # Extract using a more detailed regex
        m2 = re.match(
            r"^p\.\[([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z]);"
            r"([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z])\]$",
            hgvs,
        )
        if not m2:
            return None
        r1, p1, a1, r2, p2, a2 = m2.groups()
        a1_, r1_ = _aa(a1), _aa(r1)
        a2_, r2_ = _aa(a2), _aa(r2)
        if a1_ is None or a2_ is None:
            return None
        v1 = (int(p1), a1_) if a1_ != r1_ else None
        v2 = (int(p2), a2_) if a2_ != r2_ else None
        if v1 is None and v2 is None:
            return None
        # Pad missing side with self-pair so downstream code always gets two
        # positions; the marginal-averaging will fold synonymous positions in
        # naturally (they contribute s_hat = grand_mean).
        return (v1 or (int(p1), r1_)), (v2 or (int(p2), r2_))
    return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scored", type=Path, required=True,
                    help="output of score_haplotypes.py --with-esm")
    ap.add_argument("--scores", type=Path, required=True,
                    help="raw MaveDB scores CSV used for E_exp marginal averaging")
    ap.add_argument("--out", type=Path, required=True,
                    help="output PNG path")
    ap.add_argument("--summary", type=Path, default=None,
                    help="optional JSON summary path")
    ap.add_argument("--lof-threshold", type=float, default=0.1)
    ap.add_argument("--lof-is", choices=["low", "high"], default="low")
    args = ap.parse_args(argv)

    # ── E_exp marginal estimation from raw scores ─────────────────────────
    single_agg: dict[tuple[int, str], list[float]] = defaultdict(list)
    all_scores: list[float] = []
    with args.scores.open() as f:
        for r in csv.DictReader(f):
            hgvs = (r.get("hgvs_pro") or "").strip()
            try:
                s = float(r["score"])
            except (TypeError, ValueError, KeyError):
                continue
            pair = _parse_pair(hgvs)
            if pair is None:
                continue
            v1, v2 = pair
            all_scores.append(s)
            single_agg[v1].append(s)
            single_agg[v2].append(s)
    if not all_scores:
        print("no parseable doubles in scores file", file=sys.stderr)
        return 1
    grand_mean = sum(all_scores) / len(all_scores)
    s_hat = {k: sum(vs) / len(vs) for k, vs in single_agg.items()}

    # ── Join with hapli scored CSV ────────────────────────────────────────
    rows: list[dict] = []
    with args.scored.open() as f:
        for r in csv.DictReader(f):
            if r.get("category") != "missense_multi":
                continue
            try:
                ms = float(r["mave_score"])
                sa = float(r["s_additive"])
                sj = float(r["s_joint"])
                resid = float(r["residual"])
            except (TypeError, ValueError, KeyError):
                continue
            pair = _parse_pair(r["hgvs_pro"])
            if pair is None:
                continue
            v1, v2 = pair
            s1 = s_hat.get(v1)
            s2 = s_hat.get(v2)
            if s1 is None or s2 is None:
                continue
            e_exp = ms - s1 - s2 + grand_mean
            rows.append({
                "mave_score": ms, "S_additive": sa, "S_joint": sj,
                "residual": resid, "E_exp": e_exp,
            })
    if len(rows) < 50:
        print(f"only {len(rows)} usable rows; aborting", file=sys.stderr)
        return 1

    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr, pearsonr
    from sklearn.metrics import roc_auc_score

    ms   = np.array([r["mave_score"] for r in rows])
    sa   = np.array([r["S_additive"] for r in rows])
    sj   = np.array([r["S_joint"]    for r in rows])
    res  = np.array([r["residual"]   for r in rows])
    eexp = np.array([r["E_exp"]      for r in rows])

    # LoF labels
    lof = (ms < args.lof_threshold).astype(int) if args.lof_is == "low" \
          else (ms > args.lof_threshold).astype(int)
    auroc_add  = roc_auc_score(lof, -sa)
    auroc_joint = roc_auc_score(lof, -sj)
    auroc_res   = roc_auc_score(lof, -res)
    rho_add_eexp,   _ = spearmanr(sa, eexp)
    rho_joint_eexp, _ = spearmanr(sj, eexp)
    rho_res_eexp,   _ = spearmanr(res, eexp)
    rho_joint_mave, _ = spearmanr(-sj, ms)

    # ── Figure ────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.30)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    # Panel 1 — AUROC bars
    names = ["S_additive\n(per-variant)", "S_joint\n(haplotype)", "residual\n(hapli novelty)"]
    vals  = [auroc_add, auroc_joint, auroc_res]
    colors = ["tab:red", "tab:blue", "tab:green"]
    ax1.bar(range(3), vals, color=colors, alpha=0.85)
    ax1.set_xticks(range(3))
    ax1.set_xticklabels(names, fontsize=10)
    ax1.set_ylim(0.5, 1.0)
    ax1.set_ylabel("AUROC vs MAVE LoF label")
    ax1.axhline(0.5, color="grey", linestyle="--", linewidth=0.8)
    for i, v in enumerate(vals):
        ax1.text(i, v + 0.005, f"{v:.3f}", ha="center", fontweight="bold")
    delta = auroc_joint - auroc_add
    ax1.set_title(f"A. AUROC: joint beats additive by ΔAUROC = {delta:+.3f}",
                  fontsize=11, loc="left")

    # Panel 2 — residual vs E_exp scatter
    ax2.scatter(res, eexp, s=4, alpha=0.4, color="tab:green")
    ax2.set_xlabel("hapli residual = S_joint − S_additive")
    ax2.set_ylabel("experimental epistasis E_exp")
    ax2.axhline(0, color="grey", linestyle="--", linewidth=0.5)
    ax2.axvline(0, color="grey", linestyle="--", linewidth=0.5)
    ax2.set_title(
        f"B. residual vs experimental epistasis\n"
        f"   ρ(residual, E_exp) = {rho_res_eexp:+.3f}    "
        f"vs ρ(S_additive, E_exp) = {rho_add_eexp:+.3f}",
        fontsize=11, loc="left",
    )

    # Panel 3 — S_joint vs MAVE fitness
    ax3.scatter(-sj, ms, s=4, alpha=0.4, color="tab:blue")
    ax3.set_xlabel("−S_joint  (high = ESM predicts deleterious)")
    ax3.set_ylabel("MAVE fitness")
    ax3.set_title(
        f"C. fitness prediction\n   ρ(−S_joint, MAVE) = {rho_joint_mave:+.3f}",
        fontsize=11, loc="left",
    )

    # Panel 4 — residual distribution
    ax4.hist(res, bins=40, color="tab:green", alpha=0.7, edgecolor="black", linewidth=0.3)
    ax4.axvline(0, color="black", linestyle="--", linewidth=0.8)
    ax4.set_xlabel("residual = S_joint − S_additive")
    ax4.set_ylabel("variant count")
    n_strong = int(np.sum(np.abs(res) > 3))
    ax4.set_title(
        f"D. residual distribution\n"
        f"   mean |residual| = {np.abs(res).mean():.2f},  "
        f"|residual| > 3.0 in {n_strong}/{len(rows)} variants ({n_strong / len(rows):.1%})",
        fontsize=11, loc="left",
    )

    fig.suptitle(
        f"hapli vs variant-by-variant on real double-mutant MAVE data\n"
        f"({args.scored.stem} — {len(rows)} double mutants)",
        fontsize=13, fontweight="bold",
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote figure {args.out}")

    summary = {
        "scored_file": str(args.scored),
        "scores_file": str(args.scores),
        "n_rows": len(rows),
        "lof_threshold": args.lof_threshold,
        "auroc": {
            "S_additive": auroc_add,
            "S_joint":    auroc_joint,
            "residual":   auroc_res,
            "ΔAUROC_joint_minus_additive": auroc_joint - auroc_add,
        },
        "spearman_vs_E_exp": {
            "S_additive": rho_add_eexp,
            "S_joint":    rho_joint_eexp,
            "residual":   rho_res_eexp,
            "ratio_residual_over_additive": abs(rho_res_eexp) / max(abs(rho_add_eexp), 1e-9),
        },
        "spearman_S_joint_vs_mave": rho_joint_mave,
        "residual_distribution": {
            "mean_abs": float(np.abs(res).mean()),
            "n_with_abs_gt_3": int(np.sum(np.abs(res) > 3)),
        },
    }
    if args.summary:
        args.summary.write_text(json.dumps(summary, indent=2))
        print(f"wrote summary {args.summary}")
    print("\n=== HEADLINE ===")
    for k in ("auroc", "spearman_vs_E_exp"):
        print(f"  {k}:")
        for kk, vv in summary[k].items():
            print(f"    {kk}: {vv:+.4f}" if isinstance(vv, float) else f"    {kk}: {vv}")
    print(f"  spearman_S_joint_vs_mave: {summary['spearman_S_joint_vs_mave']:+.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
