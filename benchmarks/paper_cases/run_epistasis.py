#!/usr/bin/env python3
"""
Run hapli with --with-esm across paper cases that have ≥2 coding substitutions
on one haplotype, and record the additive-vs-joint epistasis residual for each.

These numbers land in paper Figure 2 / Table 2: controlled in-silico cases
whose residuals are large and predictably-signed by construction.

Output TSV columns:
  case, n_variants, S_additive, S_joint, residual, flagged, hap1_identity,
  hap2_identity, frameshift_restored

Usage:
  uv run python3 benchmarks/paper_cases/run_epistasis.py \\
      --cases-dir data/paper_cases \\
      --out benchmarks/paper_cases/epistasis.tsv \\
      --esm-model esm2_t6_8M_UR50D
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path


# Cases where we expect an epistasis residual signal (hap1 has ≥2 coding subs).
# 01 single missense, 02 single missense: nothing to aggregate (skipped).
# 03 stop-gain: PTC before residual calc (skipped).
# 04 compound-het LoF: variants on different haplotypes (skipped).
# 08 whole-gene DEL: no coding diff (skipped).
# 09 upstream DEL: no coding diff (skipped).
# 12 in-frame del: single event (skipped).
# 05, 06, 07, 10 are the interesting ones. 10 is a disruptive SV so may not
# produce an interpretable residual; we still run it and record whatever emits.
EPISTASIS_CASES = [
    "05_frameshift_rescue",
    "06_compound_missense_pocket",
    "07_splice_donor_cryptic",
    "10_inversion_body",
]


def _extract_epistasis(json_path: Path) -> dict:
    data = json.loads(json_path.read_text())
    ev = data.get("evidence", {})
    epi = ev.get("epistasis", [])
    prot = ev.get("protein", [])

    hap1_prot = next((p for p in prot if p.get("haplotype") == 1), {})
    hap2_prot = next((p for p in prot if p.get("haplotype") == 2), {})

    # Report the hap1 residual (that's where we put the epistasis signal by design).
    hap1_epi = next((r for r in epi if r.get("haplotype") == 1), {})

    return {
        "n_variants": hap1_epi.get("n_variants"),
        "S_additive": hap1_epi.get("s_additive"),
        "S_joint": hap1_epi.get("s_joint"),
        "residual": hap1_epi.get("residual"),
        "flagged": hap1_epi.get("flagged"),
        "hap1_identity": hap1_prot.get("identity"),
        "hap2_identity": hap2_prot.get("identity"),
        "frameshift_restored": bool((hap1_prot.get("frameshift_region") or {}).get("restored")),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases-dir", type=Path, default=Path("data/paper_cases"))
    ap.add_argument("--out", type=Path, default=Path("benchmarks/paper_cases/epistasis.tsv"))
    ap.add_argument("--esm-model", default="esm2_t6_8M_UR50D",
                    help="ESM2 checkpoint name (default: 8M small)")
    args = ap.parse_args(argv)

    repo_root = Path(__file__).resolve().parent.parent.parent
    cases_dir = args.cases_dir if args.cases_dir.is_absolute() else repo_root / args.cases_dir
    out_path = args.out if args.out.is_absolute() else repo_root / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for case in EPISTASIS_CASES:
        case_dir = cases_dir / case
        if not case_dir.exists():
            print(f"[skip] missing case bundle: {case_dir}", file=sys.stderr)
            continue
        with tempfile.TemporaryDirectory(prefix=f"hapli_epi_{case}_") as td:
            out_dir = Path(td)
            result = subprocess.run(
                ["uv", "run", "main.py", "analyze",
                 "--gene", "G1", "--sample", "S1",
                 "--vcf", str(case_dir / "phased.vcf.gz"),
                 "--reference", str(case_dir / "reference.fa"),
                 "--gff", str(case_dir / "annotation.gff3"),
                 "--output-dir", str(out_dir),
                 "--with-esm", "--esm-model", args.esm_model],
                cwd=repo_root, capture_output=True, text=True, timeout=600,
            )
            if result.returncode != 0:
                print(f"[fail] {case}: rc={result.returncode}\n{result.stderr[-500:]}",
                      file=sys.stderr)
                rows.append({"case": case, "status": "FAILED"})
                continue
            json_path = out_dir / "S1_G1_analysis.json"
            if not json_path.exists():
                rows.append({"case": case, "status": "NO_JSON"})
                continue
            summary = _extract_epistasis(json_path)
            rows.append({"case": case, "status": "ok", **summary})
            print(
                f"[ok] {case}: n={summary['n_variants']}, "
                f"S_add={summary['S_additive']}, S_joint={summary['S_joint']}, "
                f"residual={summary['residual']}, flagged={summary['flagged']}",
                file=sys.stderr,
            )

    columns = ["case", "status", "n_variants", "S_additive", "S_joint",
               "residual", "flagged", "hap1_identity", "hap2_identity",
               "frameshift_restored"]
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"\nwrote {out_path} — {len(rows)} rows", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
