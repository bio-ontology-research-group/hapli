#!/usr/bin/env python3
"""
Synthetic SV-haplotype benchmark. Constructs a tiny reference protein with
known structure (start codon, poly-Ala body, domain-boundary motif, stop)
then runs the unified scorer on SV-shaped haplotype variants expressed
in HGVS protein notation.

Validates that the Phase 5 "works for any variant shape" claim holds:
single missense, double missense, nonsense (stop-gain), frameshift,
single-residue deletion, range deletion, large N-terminal deletion
(whole-gene-loss proxy), compound LoF+missense.

Each case defines an `expect` dict whose keys are any of the HaplotypeScore
fields; the benchmark emits PASS/FAIL per assertion and a final summary.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from hapli.core.hgvs import parse_hgvs_pro
from hapli.interpretation.haplotype_scoring import score_haplotype


REF = (
    "M"                          # pos 1 (start)
    + "A" * 9                    # pos 2-10
    + "RCVTQCK"                  # pos 11-17 (structured motif)
    + "A" * 27                   # pos 18-44
    + "GGGWFK"                   # pos 45-50 (C-terminal motif)
)
assert len(REF) == 50


@dataclass
class SVCase:
    name: str
    hgvs: str
    expect: dict[str, Any] = field(default_factory=dict)
    description: str = ""


CASES = [
    # --- baselines (should behave normally) ---
    SVCase(
        name="single_missense_noop",
        hgvs="p.Cys16Ser",
        expect={"category": "missense_single", "length_changed": False,
                "premature_stop_at": None, "n_missense": 1},
        description="Control: single missense, no SV shape.",
    ),
    SVCase(
        name="synonymous",
        hgvs="p.Ala4=",
        expect={"category": "synonymous", "hap_seq_equals_ref": True,
                "length_changed": False, "n_synonymous": 1},
        description="Control: synonymous — no-op.",
    ),
    # --- stop-gain / LoF ---
    SVCase(
        # REF position 13 is 'V' (positions 11-17 = 'RCVTQCK').
        name="early_stop_gain",
        hgvs="p.Val13Ter",
        expect={"category": "lof", "length_changed": True,
                "premature_stop_at": 13, "hap_length": 12,
                "score_heuristic_max": 0.5},
        description="Stop-gain at position 13 → LoF, heuristic score ≤ 0.5.",
    ),
    SVCase(
        # REF position 48 is 'W' (positions 45-50 = 'GGGWFK').
        name="late_stop_gain",
        hgvs="p.Trp48Ter",
        expect={"category": "lof", "length_changed": True,
                "premature_stop_at": 48, "hap_length": 47},
        description="Late stop — fraction-preserved heuristic is non-zero.",
    ),
    # --- frameshift ---
    SVCase(
        name="frameshift_with_ter_offset",
        hgvs="p.Thr14GlyfsTer5",
        expect={"category": "lof", "n_frameshift": 1,
                "premature_stop_at": 19, "length_changed": True},
        description="Frameshift — LoF category, PTC reported.",
    ),
    SVCase(
        name="frameshift_minimal",
        hgvs="p.Thr14fs",
        expect={"category": "lof", "n_frameshift": 1,
                "length_changed": True},
        description="Bare frameshift notation without Ter offset.",
    ),
    # --- in-frame deletions (indel category) ---
    SVCase(
        name="single_residue_del",
        hgvs="p.Ala4del",
        expect={"category": "indel", "n_indel": 1,
                "length_changed": True, "hap_length": 49},
        description="In-frame single-residue deletion.",
    ),
    SVCase(
        name="range_del_small",
        hgvs="p.Arg11_Cys16del",
        expect={"category": "indel", "n_indel": 1,
                "length_changed": True, "hap_length": 44,
                "score_heuristic_max": 0.95},
        description="6-residue in-frame deletion.",
    ),
    SVCase(
        name="range_del_large_nterm_like_sv_whole_gene_loss",
        hgvs="p.Met1_Ala10del",
        expect={"category": "indel",
                "hap_length": 40, "length_changed": True,
                "score_heuristic_max": 0.9},
        description="Large N-terminal deletion — SV-like whole-gene loss proxy.",
    ),
    # --- compound (LoF + missense) ---
    SVCase(
        name="compound_lof_and_missense",
        hgvs="p.[Ala5Val;Val13Ter]",
        expect={"category": "lof",
                "n_nonsense": 1, "n_missense": 1,
                "length_changed": True, "premature_stop_at": 13},
        description="Compound: stop-gain + missense on same hap → LoF dominates.",
    ),
    # --- out-of-range / skipped ---
    SVCase(
        name="out_of_range_skipped",
        hgvs="p.Ala999Val",
        expect={"category": "unmodified", "n_skipped": 1},
        description="Variant outside the reference — skipped, not errored.",
    ),
    # --- residual is None for all length-changing cases ---
    SVCase(
        name="residual_is_none_for_truncation",
        hgvs="p.Val13Ter",
        expect={"residual_is_none": True},
        description="ESM residual is undefined for length-changing haplotypes; the shape signals carry the information.",
    ),
]


def _check(score, case: SVCase) -> list[tuple[str, bool, str]]:
    out: list[tuple[str, bool, str]] = []
    for k, v in case.expect.items():
        if k == "score_heuristic_max":
            hs = _heuristic_score(score, len(REF))
            out.append((k, hs <= v, f"heuristic={hs:.3f} <= {v}"))
        elif k == "hap_seq_equals_ref":
            out.append((k, (score.hap_seq == REF) == v, f"hap_seq == REF = {score.hap_seq == REF}"))
        elif k == "residual_is_none":
            out.append((k, (score.residual is None) == v, f"residual={score.residual}"))
        else:
            actual = getattr(score, k)
            out.append((k, actual == v, f"{k}={actual} expected {v}"))
    return out


def _heuristic_score(score, ref_len: int) -> float:
    if score.premature_stop_at is not None:
        frac = score.premature_stop_at / ref_len if ref_len else 0
        if frac <= 0.1:
            return 0.0
        return score.identity * frac
    return score.identity


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    results: list[dict] = []
    fail = 0
    total_checks = 0
    for case in CASES:
        variants = parse_hgvs_pro(case.hgvs)
        score = score_haplotype(REF, case.hgvs, variants)
        checks = _check(score, case)
        case_fails = [c for c in checks if not c[1]]
        status = "PASS" if not case_fails else "FAIL"
        results.append({
            "name": case.name,
            "hgvs": case.hgvs,
            "category": score.category,
            "hap_length": score.hap_length,
            "status": status,
            "checks": [{"check": c[0], "ok": c[1], "detail": c[2]} for c in checks],
        })
        total_checks += len(checks)
        fail += len(case_fails)
        marker = "✓" if status == "PASS" else "✗"
        print(f"[{marker}] {case.name:45s}  {case.hgvs:30s}  → {score.category} "
              f"len={score.hap_length} pts={score.premature_stop_at}")
        for c_name, ok, detail in case_fails:
            print(f"        FAIL: {c_name}: {detail}")

    summary_path = args.out / "sv_benchmark.json"
    summary_path.write_text(json.dumps({
        "reference_length": len(REF),
        "n_cases": len(CASES),
        "n_assertions": total_checks,
        "n_failures": fail,
        "pass_rate": (total_checks - fail) / total_checks,
        "cases": results,
    }, indent=2))
    print(f"\n{len(CASES)} cases, {total_checks} assertions, {fail} failures")
    print(f"wrote {summary_path}")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
