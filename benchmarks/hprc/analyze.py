#!/usr/bin/env python3
"""
HPRC R2 × Mode B post-aggregation analysis.

Produces the two HPRC-specific paper deliverables:

  1. **Per-gene presence distribution** across the HPRC cohort: how many
     haplotype records per gene are classified intact / low_identity /
     partial / deleted / duplicated by Liftoff. This is the Axis-1 scale
     demo: gene-level structural signal across real assembled diploids.

  2. **Per-super-population LoF allele frequency**: group samples by
     AFR / AMR / EAS / EUR / SAS (from the sample manifest) and report
     per-(gene, super-pop) LoF AF. Validates that hapli's Mode B output
     is calibrated across ancestries.

Usage:
  uv run python3 benchmarks/hprc/analyze.py \\
      --per-sample results/hprc/aggregate/per_sample.tsv \\
      --per-gene   results/hprc/aggregate/per_gene.tsv \\
      --sample-manifest benchmarks/hprc/sample_list.tsv \\
      --out        benchmarks/hprc/results/
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


def _load_manifest(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "sample_id":
                continue
            if len(parts) < 3:
                continue
            out[parts[0]] = {
                "hap1_url": parts[1],
                "hap2_url": parts[2],
                "population": parts[3] if len(parts) > 3 else "",
                "super_population": parts[4] if len(parts) > 4 else "",
            }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-sample", type=Path, required=True)
    ap.add_argument("--per-gene",   type=Path, required=True)
    ap.add_argument("--sample-manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    manifest = _load_manifest(args.sample_manifest)

    with args.per_sample.open() as f:
        per_sample = list(csv.DictReader(f, delimiter="\t"))
    with args.per_gene.open() as f:
        per_gene = list(csv.DictReader(f, delimiter="\t"))

    # ─ Table A: per-gene presence-status distribution ────────────────────────
    presence_counts: dict[str, Counter] = defaultdict(Counter)
    for r in per_sample:
        gene = r["gene"]
        for hap_col in ("presence_hap1", "presence_hap2"):
            status = r.get(hap_col, "") or "unknown"
            presence_counts[gene][status] += 1

    table_a = args.out / "per_gene_presence.tsv"
    with table_a.open("w") as f:
        # Headers: gene, intact, low_identity, partial, deleted, duplicated, uncertain, not_run
        statuses = ["intact", "low_identity", "partial", "deleted", "duplicated",
                    "uncertain", "not_run"]
        f.write("gene\t" + "\t".join(statuses) + "\ttotal_hap_records\n")
        for gene in sorted(presence_counts):
            counts = presence_counts[gene]
            row = [gene]
            tot = sum(counts.values())
            for s in statuses:
                row.append(str(counts.get(s, 0)))
            row.append(str(tot))
            f.write("\t".join(row) + "\n")
    print(f"wrote {table_a}", file=sys.stderr)

    # ─ Table B: per-(gene, super-pop) LoF AF ─────────────────────────────────
    lof_alleles_by_gene_pop: dict[tuple[str, str], list[int]] = defaultdict(list)
    total_by_gene_pop: dict[tuple[str, str], int] = defaultdict(int)
    for r in per_sample:
        sample = r["sample"]
        gene = r["gene"]
        super_pop = manifest.get(sample, {}).get("super_population", "UNKNOWN") or "UNKNOWN"
        # Count each haplotype that's LoF (score ≤ 0.5) as one allele
        for hap_col, score_col in (("presence_hap1", "hap1_score"),
                                   ("presence_hap2", "hap2_score")):
            try:
                s = float(r.get(score_col, "") or "nan")
            except ValueError:
                s = float("nan")
            total_by_gene_pop[(gene, super_pop)] += 1
            if s == s and s <= 0.5:
                lof_alleles_by_gene_pop[(gene, super_pop)] = \
                    lof_alleles_by_gene_pop.get((gene, super_pop), []) + [1]

    table_b = args.out / "per_gene_per_superpop_lof_af.tsv"
    with table_b.open("w") as f:
        f.write("gene\tsuper_population\tn_lof_alleles\tn_total_alleles\tlof_allele_freq\n")
        for (gene, pop), tot in sorted(total_by_gene_pop.items()):
            n_lof = len(lof_alleles_by_gene_pop.get((gene, pop), []))
            af = n_lof / tot if tot else 0.0
            f.write(f"{gene}\t{pop}\t{n_lof}\t{tot}\t{af:.6f}\n")
    print(f"wrote {table_b}", file=sys.stderr)

    # ─ Tally for paper abstract ──────────────────────────────────────────────
    total_hap = sum(sum(c.values()) for c in presence_counts.values())
    sv_presence_records = sum(
        c.get("deleted", 0) + c.get("duplicated", 0) + c.get("low_identity", 0)
        + c.get("partial", 0)
        for c in presence_counts.values()
    )
    compound_het_samples = sum(
        1 for r in per_sample
        if r.get("compound_het_lof") in ("1", "True", "true")
    )
    tally = {
        "n_samples": len({r["sample"] for r in per_sample}),
        "n_genes": len({r["gene"] for r in per_sample}),
        "n_hap_records_total": total_hap,
        "n_sv_presence_records": sv_presence_records,
        "sv_presence_rate": sv_presence_records / total_hap if total_hap else 0.0,
        "n_compound_het_lof_samplegenes": compound_het_samples,
        "super_populations_covered": sorted(
            {manifest.get(r["sample"], {}).get("super_population", "UNKNOWN")
             for r in per_sample}
        ),
    }
    (args.out / "tally.json").write_text(json.dumps(tally, indent=2))
    print(json.dumps(tally, indent=2), file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
