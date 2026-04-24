#!/usr/bin/env python3
"""
1000G + ACMG-SF-v3.2 post-aggregation analysis.

Consumes the per-sample and per-gene TSVs produced by `hapli aggregate` at
the end of the Snakemake run, and produces the two paper deliverables:

  1. **Paper Table 2** — per-gene summary over the cohort:
       gene, n_samples, LoF_allele_freq, n_compound_het_lof_samples,
       n_svpresence_hap_records, pLI (if available), ACMG_inheritance,
       disagreement_class (axis-1 / axis-2 / both / none).

  2. **Paper Figure 3** — scatter of hapli's per-gene LoF allele
     frequency against gnomAD's published per-gene observed/expected
     LoF ratio (oe_lof), highlighting the ACMG SF v3.2 subset.

Claims the paper makes from these outputs:

  * "Across 2,504 1000G Phase 3 individuals over the 79 ACMG SF v3.2
     genes, hapli emits N compound-het-LoF genotype calls (axis-2)
     that bcftools/csq emits as 2N separate per-variant records
     without diploid aggregation."
  * "Across the same cohort, hapli's Liftoff-derived presence flag
     recovers M gene-level SV events (axis-1: deleted / low_identity /
     duplicated) that the per-variant tools are structurally silent on."
  * "Per-gene LoF allele frequency agrees with gnomAD's per-gene
     observed LoF count (Pearson r ≥ X), validating that the diploid
     aggregation is not systematically mis-calibrated."

Usage:
  uv run python3 benchmarks/1000g_acmg/analyze.py \\
      --per-sample results/1000g_acmg/aggregate/per_sample.tsv \\
      --per-gene   results/1000g_acmg/aggregate/per_gene.tsv \\
      --gene-list  benchmarks/1000g_acmg/acmg_sf_v3_2.tsv \\
      --gnomad-constraint data/gnomad/gnomad_v4_constraint_metrics.tsv  \\
      --out        benchmarks/1000g_acmg/results/
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def _load_acmg_list(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    with path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["hgnc_symbol"]] = {
                "disease_category": row.get("disease_category", ""),
                "inheritance": row.get("inheritance", ""),
            }
    return out


def _load_per_gene(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["gene"]] = row
    return out


def _load_gnomad_constraint(path: Path) -> dict[str, dict[str, float]]:
    """Pull per-gene pLI / oe_lof from a gnomAD constraint TSV.
    Columns: gene, canonical, pLI, oe_lof, mis_z (subset; gnomAD v4 format)."""
    out: dict[str, dict[str, float]] = {}
    if not path or not path.exists():
        return out
    with path.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            g = row.get("gene") or row.get("gene_symbol") or ""
            if not g:
                continue
            if row.get("canonical") not in ("true", "TRUE", "1", "True", None, ""):
                # gnomAD flags canonical transcripts; skip non-canonical rows.
                continue
            try:
                pli = float(row.get("pLI", "") or "nan")
                oe = float(row.get("oe_lof", "") or "nan")
            except (ValueError, TypeError):
                pli = float("nan")
                oe = float("nan")
            # If duplicated, take the first canonical row.
            if g not in out:
                out[g] = {"pLI": pli, "oe_lof": oe}
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-sample", type=Path, required=True)
    ap.add_argument("--per-gene", type=Path, required=True)
    ap.add_argument("--gene-list", type=Path, required=True,
                    help="TSV of ACMG-SF-v3.2 genes with inheritance / category columns.")
    ap.add_argument("--gnomad-constraint", type=Path, default=None,
                    help="Optional gnomAD v4 per-gene constraint TSV.")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output directory for table + figure.")
    args = ap.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)

    acmg = _load_acmg_list(args.gene_list)
    per_gene = _load_per_gene(args.per_gene)
    gnomad = _load_gnomad_constraint(args.gnomad_constraint) if args.gnomad_constraint else {}

    # ─ Build the per-gene enriched table ─────────────────────────────────────
    rows: list[dict] = []
    for sym, info in sorted(acmg.items()):
        g = per_gene.get(sym, {})
        # Axis-2 disagreement: samples with compound_het_lof=True on this gene
        n_chl = int(g.get("n_compound_het_lof_samples", 0) or 0)
        # Axis-1 (SV presence) disagreement: would need a separate count —
        # for now, read from the optional `n_sv_presence_hap_records` column
        # that the aggregator emits when it's present.
        n_sv = int(g.get("n_sv_presence_hap_records", 0) or 0)
        lof_af = float(g.get("lof_allele_freq", 0.0) or 0.0)
        n_samples = int(g.get("n_samples", 0) or 0)
        n_hap = int(g.get("n_hap_records", 0) or 0)
        n_lof_alleles = int(g.get("n_lof_alleles", 0) or 0)

        g4 = gnomad.get(sym, {})
        pli = g4.get("pLI")
        oe_lof = g4.get("oe_lof")

        if n_chl > 0 and n_sv > 0:
            disagreement = "axis-1+2"
        elif n_chl > 0:
            disagreement = "axis-2"
        elif n_sv > 0:
            disagreement = "axis-1"
        else:
            disagreement = "none"

        rows.append({
            "gene": sym,
            "disease_category": info["disease_category"],
            "inheritance": info["inheritance"],
            "n_samples": n_samples,
            "n_hap_records": n_hap,
            "n_lof_alleles": n_lof_alleles,
            "lof_allele_freq": lof_af,
            "n_compound_het_lof_samples": n_chl,
            "n_sv_presence_hap_records": n_sv,
            "disagreement_class": disagreement,
            "gnomad_pLI": f"{pli:.3f}" if pli is not None and pli == pli else "",
            "gnomad_oe_lof": f"{oe_lof:.3f}" if oe_lof is not None and oe_lof == oe_lof else "",
        })

    out_tbl = args.out / "per_gene_summary.tsv"
    columns = list(rows[0].keys()) if rows else ["gene"]
    with out_tbl.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {out_tbl}", file=sys.stderr)

    # ─ Tallies for the paper abstract ────────────────────────────────────────
    tally = {
        "n_genes_run": len(rows),
        "n_genes_axis2_disagreement": sum(1 for r in rows if r["n_compound_het_lof_samples"] > 0),
        "n_genes_axis1_disagreement": sum(1 for r in rows if r["n_sv_presence_hap_records"] > 0),
        "total_compound_het_lof_events": sum(r["n_compound_het_lof_samples"] for r in rows),
        "total_sv_presence_events": sum(r["n_sv_presence_hap_records"] for r in rows),
        "n_samples_max": max((r["n_samples"] for r in rows), default=0),
    }
    (args.out / "tally.json").write_text(json.dumps(tally, indent=2))
    print(json.dumps(tally, indent=2), file=sys.stderr)

    # ─ Paper Figure 3: per-gene LoF AF vs gnomAD oe_lof ──────────────────────
    if gnomad and any(r["gnomad_oe_lof"] for r in rows):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not available; skipping figure", file=sys.stderr)
        else:
            xs = [float(r["gnomad_oe_lof"]) for r in rows if r["gnomad_oe_lof"]]
            ys = [r["lof_allele_freq"] for r in rows if r["gnomad_oe_lof"]]
            labels = [r["gene"] for r in rows if r["gnomad_oe_lof"]]
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(xs, ys, s=24, alpha=0.7)
            ax.set_xlabel("gnomAD v4 oe_lof (per-gene)")
            ax.set_ylabel("hapli 1000G LoF allele freq")
            ax.set_title("hapli 1000G Phase 3 LoF AF vs. gnomAD oe_lof, ACMG-SF-v3.2 genes")
            # Annotate high-LoF genes
            for x, y, lab in zip(xs, ys, labels):
                if y > 0.001:
                    ax.annotate(lab, (x, y), fontsize=7, alpha=0.8)
            fig_path = args.out / "lof_af_vs_oe_lof.png"
            fig.tight_layout()
            fig.savefig(fig_path, dpi=120)
            print(f"wrote {fig_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
