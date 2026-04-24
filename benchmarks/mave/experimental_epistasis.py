#!/usr/bin/env python3
"""
Compute experimental epistasis on a pairwise DMS scoreset.

Following Olson 2014 (GB1 pairwise scan): since the scoreset contains no pure
single-residue mutants, we estimate the marginal fitness of each single
substitution by averaging over all pairs containing it. The "wildtype"
baseline is the grand mean.

For each double haplotype with substitutions (v1, v2) from the scoreset:

    s_hat(v) := mean over all observed pairs containing v of s(v, .)
    s_bar    := grand mean of all observed pairs

    E_exp(v1, v2) := s(v1, v2) - s_hat(v1) - s_hat(v2) + s_bar

E_exp > 0 → observed double is MORE fit than additive prediction (positive epistasis / rescue)
E_exp < 0 → observed double is LESS fit than additive prediction (negative / synergistic epistasis)
E_exp ≈ 0 → additive; no epistasis

This yields a ground-truth epistasis target for each haplotype. We then
correlate hapli's residual (from the unified scorer, `s_joint - s_additive`
in ESM log-odds) against E_exp. The paper's headline claim is that hapli
tracks this experimental epistasis signal — variant-by-variant predictors
(by construction) produce zero on this axis.

Usage:
  ./experimental_epistasis.py --scores raw-scores.csv --scored hapli-scored.csv --out analysis.json
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path


# Match p.[Thr10Gly;Val38Tyr]-style double mutants.
_RE_DOUBLE_BRACKET = re.compile(
    r"^p\.\["
    r"([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z])"
    r";"
    r"([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z])"
    r"\]$"
)

# Match p.Thr10_Leu11delinsXxxYyy-style adjacent pairs — Olson 2014's
# position-adjacent doubles are encoded this way.
_RE_ADJ_DELINS = re.compile(
    r"^p\.([A-Z][a-z]{2}|[A-Z])(\d+)"
    r"_([A-Z][a-z]{2}|[A-Z])(\d+)"
    r"delins([A-Za-z]+)$"
)

_AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def _aa1(code: str) -> str | None:
    if code in _AA3_TO_1:
        return _AA3_TO_1[code]
    if len(code) == 1 and code in "ACDEFGHIKLMNPQRSTVWY":
        return code
    return None


def parse_double(hgvs: str) -> tuple[tuple[int, str], tuple[int, str]] | None:
    """Parse a double mutant into two (position, alt_aa) pairs, or None.
    Accepts both bracket form and adjacent-delins form. Synonymous positions
    (ref == alt) are dropped (they're not substitutions)."""
    m = _RE_DOUBLE_BRACKET.match(hgvs)
    if m:
        r1, p1, a1, r2, p2, a2 = m.groups()
        a1_ = _aa1(a1); a2_ = _aa1(a2)
        r1_ = _aa1(r1); r2_ = _aa1(r2)
        if a1_ is None or a2_ is None:
            return None
        v1 = (int(p1), a1_) if a1_ != r1_ else None
        v2 = (int(p2), a2_) if a2_ != r2_ else None
        if v1 is None and v2 is None:
            return None
        return v1, v2  # one may be None if synonymous

    m = _RE_ADJ_DELINS.match(hgvs)
    if m:
        r1, p1, r2, p2, alt_text = m.groups()
        # Expect two residues in alt_text; support 3-letter (6 chars) or 1-letter (2 chars).
        alt = _decode_aa_string(alt_text)
        if alt is None or len(alt) != 2:
            return None
        r1_, r2_ = _aa1(r1), _aa1(r2)
        v1 = (int(p1), alt[0]) if alt[0] != r1_ else None
        v2 = (int(p2), alt[1]) if alt[1] != r2_ else None
        if v1 is None and v2 is None:
            return None
        return v1, v2

    return None


def _decode_aa_string(s: str) -> str | None:
    if all(c in "ACDEFGHIKLMNPQRSTVWY" for c in s):
        return s
    if len(s) % 3 != 0:
        return None
    out = []
    for i in range(0, len(s), 3):
        c = _aa1(s[i : i + 3])
        if c is None:
            return None
        out.append(c)
    return "".join(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", type=Path, required=True,
                    help="raw MaveDB scores CSV (every row is a (v1, v2) + experimental fitness)")
    ap.add_argument("--scored", type=Path, required=True,
                    help="hapli-scored CSV produced by score_haplotypes.py --with-esm")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--figure", type=Path, default=None)
    args = ap.parse_args(argv)

    # ── Pass 1 over raw scores: build marginal estimates ────────────────────
    # single_fitness[(pos, aa)] = (sum_of_doubles, count)
    single_agg: dict[tuple[int, str], list[float]] = defaultdict(list)
    all_scores: list[float] = []
    doubles_parsed: dict[str, tuple[float, tuple[int, str], tuple[int, str]]] = {}

    with args.scores.open() as f:
        for r in csv.DictReader(f):
            hgvs = (r.get("hgvs_pro") or "").strip()
            if not hgvs or hgvs == "NA":
                continue
            try:
                score = float(r["score"])
            except (ValueError, KeyError, TypeError):
                continue
            parsed = parse_double(hgvs)
            if parsed is None:
                continue
            v1, v2 = parsed
            all_scores.append(score)
            if v1 is not None:
                single_agg[v1].append(score)
            if v2 is not None:
                single_agg[v2].append(score)
            if v1 is not None and v2 is not None:
                doubles_parsed[hgvs] = (score, v1, v2)

    if not all_scores:
        print("no parseable rows", file=sys.stderr)
        return 1

    grand_mean = sum(all_scores) / len(all_scores)
    s_hat = {k: sum(vs) / len(vs) for k, vs in single_agg.items()}

    print(
        f"raw pass: {len(all_scores)} rows, {len(single_agg)} distinct single "
        f"positions covered, grand mean = {grand_mean:.4f}",
        file=sys.stderr,
    )

    # ── Pass 2 over hapli-scored CSV: join + compute E_exp ──────────────────
    rows_out: list[dict] = []
    with args.scored.open() as f:
        for r in csv.DictReader(f):
            hgvs = r["hgvs_pro"]
            if r.get("category") != "missense_multi":
                continue
            try:
                mave_score = float(r["mave_score"])
                residual = float(r["residual"])
                s_add = float(r["s_additive"])
                s_joint = float(r["s_joint"])
            except (ValueError, KeyError, TypeError):
                continue
            if hgvs not in doubles_parsed:
                continue
            _, v1, v2 = doubles_parsed[hgvs]
            s1 = s_hat.get(v1)
            s2 = s_hat.get(v2)
            if s1 is None or s2 is None:
                continue
            e_exp = mave_score - s1 - s2 + grand_mean
            # The "linear-additive prediction" — what variant-by-variant would produce
            linear_add = s1 + s2 - grand_mean
            rows_out.append({
                "hgvs_pro": hgvs,
                "mave_score": mave_score,
                "linear_additive_prediction": linear_add,
                "E_exp": e_exp,
                "hapli_s_additive": s_add,
                "hapli_s_joint": s_joint,
                "hapli_residual": residual,
            })

    if len(rows_out) < 20:
        print(f"only {len(rows_out)} rows joined; need more data", file=sys.stderr)
        return 1

    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr, spearmanr

    e_exp = np.array([r["E_exp"] for r in rows_out])
    residual = np.array([r["hapli_residual"] for r in rows_out])
    s_joint = np.array([r["hapli_s_joint"] for r in rows_out])
    s_add = np.array([r["hapli_s_additive"] for r in rows_out])
    linear = np.array([r["linear_additive_prediction"] for r in rows_out])
    mave = np.array([r["mave_score"] for r in rows_out])

    def corr(name, x, y):
        rs, ps = spearmanr(x, y); rp, pp = pearsonr(x, y)
        return {"name": name, "spearman": float(rs), "spearman_p": float(ps),
                "pearson": float(rp), "pearson_p": float(pp), "n": int(len(x))}

    summary = {
        "scores_file": str(args.scores),
        "scored_file": str(args.scored),
        "n_rows": len(rows_out),
        "grand_mean_fitness": grand_mean,
        "n_singles_covered": len(single_agg),
        "headline_correlations": {
            # The point of this benchmark: hapli.residual tracks experimental epistasis;
            # variant-by-variant (s_additive alone) does not.
            "residual_vs_E_exp":           corr("residual_vs_E_exp", residual, e_exp),
            "s_joint_vs_E_exp":            corr("s_joint_vs_E_exp", s_joint, e_exp),
            "s_additive_vs_E_exp":         corr("s_additive_vs_E_exp", s_add, e_exp),
            # Sanity: the linear-additive fitness prediction has NO information
            # about epistasis by construction — this should be ~0 correlation.
            "linear_additive_vs_E_exp":    corr("linear_additive_vs_E_exp", linear, e_exp),
            # Raw fitness vs E_exp — a positive correlation just reflects that
            # epistasis shifts measured fitness.
            "mave_score_vs_E_exp":         corr("mave_score_vs_E_exp", mave, e_exp),
        },
        "residual_distribution": {
            "mean_abs": float(np.abs(residual).mean()),
            "mean_abs_E_exp": float(np.abs(e_exp).mean()),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    fig_path = args.figure or args.out.with_suffix(".png")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    axes[0].scatter(residual, e_exp, s=4, alpha=0.4, color="tab:blue")
    axes[0].set_xlabel("hapli residual (S_joint − S_additive)")
    axes[0].set_ylabel("experimental epistasis E_exp")
    axes[0].set_title(
        f"residual vs E_exp   Spearman={spearmanr(residual, e_exp)[0]:.3f}\n"
        f"(haplotype-level signal, the paper's claim)"
    )
    axes[1].scatter(s_add, e_exp, s=4, alpha=0.4, color="tab:red")
    axes[1].set_xlabel("hapli S_additive (variant-by-variant baseline)")
    axes[1].set_ylabel("experimental epistasis E_exp")
    axes[1].set_title(
        f"S_additive vs E_exp   Spearman={spearmanr(s_add, e_exp)[0]:.3f}\n"
        f"(what variant-by-variant tools see)"
    )
    axes[2].hist(e_exp, bins=40, color="tab:grey", alpha=0.7)
    axes[2].axvline(0, color="black", linestyle="--")
    axes[2].set_xlabel("experimental epistasis E_exp")
    axes[2].set_ylabel("count")
    axes[2].set_title(f"E_exp distribution\nmean |E_exp| = {np.abs(e_exp).mean():.3f}")
    fig.suptitle(
        f"{args.scored.stem} — {len(rows_out)} doubles   "
        f"(residual-vs-E_exp is the paper's headline)"
    )
    fig.tight_layout()
    fig.savefig(fig_path, dpi=120)
    print(f"wrote figure {fig_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
