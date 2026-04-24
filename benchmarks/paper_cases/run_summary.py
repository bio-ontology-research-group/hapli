#!/usr/bin/env python3
"""
Run hapli's Mode A pipeline on every paper case and emit a summary TSV.

Becomes paper Table 1: a per-case row showing what hapli's evidence bundle
looks like for each motivating variant class. Directly comparable to what
bcftools/csq would produce, which is why the TSV also captures per-variant
consequences when available.

Output columns:
  case, hap1_presence, hap2_presence, hap1_score, hap2_score,
  compound_het_lof, n_consequences, hap1_categories, hap2_categories,
  hap1_identity, hap2_identity, hap1_ptc, hap2_ptc, frameshift_restored,
  joint_residual, flagged_epistasis

Usage:
  uv run python3 benchmarks/paper_cases/run_summary.py \\
      --cases-dir data/paper_cases \\
      --out benchmarks/paper_cases/summary.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path


CASE_NAMES = [
    "01_synonymous_pair",
    "02_single_benign_missense",
    "03_single_stop_gain",
    "04_compound_het_lof",
    "05_frameshift_rescue",
    "06_compound_missense_pocket",
    "07_splice_donor_cryptic",
    "08_symbolic_del_whole_gene",
    "09_upstream_del_shift",
    "10_inversion_body",
    "11_tandem_dup_small_gene",
    "12_inframe_domain_indel",
]


def _cats_for_hap(consequences: list, hap: int) -> str:
    cats = sorted({c.get("consequence") for c in consequences if c.get("haplotype") == hap})
    return ",".join(c for c in cats if c) or "-"


def _protein_for_hap(proteins: list, hap: int) -> dict:
    for p in proteins:
        if p.get("haplotype") == hap:
            return p
    return {}


def _summarise(json_path: Path) -> dict:
    data = json.loads(json_path.read_text())
    ev = data.get("evidence", {})
    pres = ev.get("presence", {})
    dip = ev.get("diploid", {})
    cons = ev.get("consequence", [])
    prot = ev.get("protein", [])
    epi = ev.get("epistasis", [])

    p1 = _protein_for_hap(prot, 1)
    p2 = _protein_for_hap(prot, 2)
    frameshift_restored = bool(p1.get("frameshift_region") or {}) and \
                          (p1.get("frameshift_region") or {}).get("restored", False)

    joint_residual = None
    flagged = False
    for r in epi:
        if r.get("haplotype") == 1 and r.get("residual") is not None:
            joint_residual = r["residual"]
            flagged = bool(r.get("flagged"))
            break

    def _pres_status(h: str) -> str:
        p = pres.get(h)
        if isinstance(p, dict):
            return p.get("status", "?")
        return str(p) if p else "-"

    return {
        "hap1_presence": _pres_status("hap1"),
        "hap2_presence": _pres_status("hap2"),
        "hap1_score": dip.get("hap1_score"),
        "hap2_score": dip.get("hap2_score"),
        "compound_het_lof": dip.get("compound_het_lof"),
        "n_consequences": len(cons),
        "hap1_categories": _cats_for_hap(cons, 1),
        "hap2_categories": _cats_for_hap(cons, 2),
        "hap1_identity": p1.get("identity"),
        "hap2_identity": p2.get("identity"),
        "hap1_ptc": p1.get("premature_stop_at"),
        "hap2_ptc": p2.get("premature_stop_at"),
        "frameshift_restored": frameshift_restored,
        "joint_residual": joint_residual,
        "flagged_epistasis": flagged,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases-dir", type=Path, default=Path("data/paper_cases"))
    ap.add_argument("--out", type=Path, default=Path("benchmarks/paper_cases/summary.tsv"))
    ap.add_argument("--regenerate", action="store_true",
                    help="regenerate paper-case bundles first")
    args = ap.parse_args(argv)

    repo_root = Path(__file__).resolve().parent.parent.parent
    cases_dir = args.cases_dir if args.cases_dir.is_absolute() else repo_root / args.cases_dir
    out_path = args.out if args.out.is_absolute() else repo_root / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if args.regenerate or not cases_dir.exists():
        subprocess.run(
            [sys.executable, str(repo_root / "scripts" / "generate_paper_cases.py"),
             "--all", "--out", str(cases_dir), "--force"],
            check=True,
        )

    rows = []
    for case in CASE_NAMES:
        case_dir = cases_dir / case
        if not case_dir.exists():
            print(f"[warn] missing case bundle: {case_dir}", file=sys.stderr)
            continue
        with tempfile.TemporaryDirectory(prefix=f"hapli_{case}_") as td:
            out_dir = Path(td)
            result = subprocess.run(
                ["uv", "run", "main.py", "analyze",
                 "--gene", "G1", "--sample", "S1",
                 "--vcf", str(case_dir / "phased.vcf.gz"),
                 "--reference", str(case_dir / "reference.fa"),
                 "--gff", str(case_dir / "annotation.gff3"),
                 "--output-dir", str(out_dir)],
                cwd=repo_root, capture_output=True, text=True, timeout=240,
            )
            if result.returncode != 0:
                print(f"[fail] {case}: rc={result.returncode}\n{result.stderr}",
                      file=sys.stderr)
                rows.append({"case": case, "status": "FAILED"})
                continue
            json_path = out_dir / "S1_G1_analysis.json"
            if not json_path.exists():
                rows.append({"case": case, "status": "NO_JSON"})
                continue
            summary = _summarise(json_path)
            rows.append({"case": case, "status": "ok", **summary})

    columns = [
        "case", "status",
        "hap1_presence", "hap2_presence",
        "hap1_score", "hap2_score", "compound_het_lof",
        "hap1_identity", "hap2_identity",
        "hap1_ptc", "hap2_ptc",
        "frameshift_restored",
        "n_consequences", "hap1_categories", "hap2_categories",
        "joint_residual", "flagged_epistasis",
    ]
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"wrote {out_path} — {len(rows)} rows", file=sys.stderr)
    # Echo the table to stdout
    w2 = csv.DictWriter(sys.stdout, fieldnames=columns, delimiter="\t", extrasaction="ignore")
    w2.writeheader()
    for r in rows:
        w2.writerow(r)
    return 0


if __name__ == "__main__":
    sys.exit(main())
