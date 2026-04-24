#!/usr/bin/env python3
"""
Fetch the 1000G NYGC high-coverage phased WGS VCFs for the ACMG SF v3.2
gene regions and write a single merged, regioned, bgzipped-tabixed VCF.

Scale-aware: instead of downloading the full 500+ GB 3202-sample release,
uses `bcftools view` remote-region fetch (tabix over HTTP) to pull only
the ACMG gene intervals from the NYGC release, per-chromosome, and merges
them into one on-disk VCF.

Inputs (public):
  - 1000G NYGC GRCh38 phased release (Byrska-Bishop 2022):
      http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
  - GENCODE v45 GFF3 (for gene-to-coordinate resolution):
      https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gff3.gz

Outputs (in --out-dir):
  - acmg_regions.bed         — per-gene intervals with 1 kb buffer
  - acmg.phased.vcf.gz[.tbi] — subset VCF across all ACMG regions
  - samples.txt              — 1000G sample list (2504 unrelated + optional relatives)

Usage:
  uv run python3 benchmarks/1000g_acmg/fetch.py \\
      --genes benchmarks/1000g_acmg/acmg_sf_v3_2.tsv \\
      --gff /path/to/gencode.v45.annotation.gff3.gz \\
      --out-dir data/1000g_acmg/ \\
      --threads 4

  # To just plan the fetch without running:
  uv run python3 benchmarks/1000g_acmg/fetch.py --dry-run ...
"""
from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
from pathlib import Path


# NYGC 2022 3202-sample phased release. One per chromosome.
NYGC_BASE = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
)


def _gene_intervals(gff_path: Path, gene_symbols: set[str], buffer_bp: int) -> list[tuple[str, int, int, str]]:
    """Scan a GFF3 for gene features whose Name or gene_name attribute is in
    `gene_symbols`. Return [(chrom, start_1based, end_1based, symbol)].
    Buffer is added symmetrically; start is clamped at 1.
    """
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    found: dict[str, tuple[str, int, int]] = {}
    with opener(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs: dict[str, str] = {}
            for kv in parts[8].rstrip(";").split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k.strip()] = v.strip()
            name = attrs.get("gene_name") or attrs.get("Name") or ""
            if name not in gene_symbols:
                continue
            if name in found:
                # Some genes have multiple entries (e.g. biotype duplicates);
                # take the widest span.
                old = found[name]
                found[name] = (old[0], min(old[1], int(parts[3])),
                               max(old[2], int(parts[4])))
            else:
                found[name] = (parts[0], int(parts[3]), int(parts[4]))
    out = []
    for sym, (chrom, s, e) in found.items():
        out.append((chrom, max(1, s - buffer_bp), e + buffer_bp, sym))
    out.sort(key=lambda r: (r[0], r[1]))
    return out


def _merge_contig_intervals(regions: list[tuple[str, int, int, str]]) -> list[tuple[str, int, int]]:
    """Collapse overlapping/adjacent intervals on the same contig."""
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, s, e, _ in regions:
        by_chrom.setdefault(chrom, []).append((s, e))
    out: list[tuple[str, int, int]] = []
    for chrom, ivs in by_chrom.items():
        ivs.sort()
        cur_s, cur_e = ivs[0]
        for s, e in ivs[1:]:
            if s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                out.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        out.append((chrom, cur_s, cur_e))
    return out


def _write_bed(regions: list[tuple[str, int, int, str]], bed_path: Path) -> None:
    with bed_path.open("w") as f:
        for chrom, s, e, sym in regions:
            f.write(f"{chrom}\t{s - 1}\t{e}\t{sym}\n")


def _chrom_vcf_url(chrom: str) -> str:
    # NYGC files are named e.g. CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz
    return (
        f"{NYGC_BASE}"
        f"CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"
    )


def _run(cmd: list[str], dry: bool) -> None:
    print("+", " ".join(cmd), file=sys.stderr)
    if dry:
        return
    subprocess.run(cmd, check=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", type=Path, required=True,
                    help="TSV with hgnc_symbol column (ACMG list)")
    ap.add_argument("--gff", type=Path, required=True,
                    help="GENCODE v45 annotation GFF3 (gzipped or plain)")
    ap.add_argument("--out-dir", type=Path, required=True,
                    help="Directory for downloaded + merged outputs")
    ap.add_argument("--buffer-bp", type=int, default=1000,
                    help="Flank to pad each gene interval (default 1 kb)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--nygc-base", default=NYGC_BASE,
                    help="Override NYGC base URL if your mirror differs")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print the plan; don't actually fetch")
    args = ap.parse_args(argv)

    # Read gene symbols from column `hgnc_symbol`.
    with args.genes.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "hgnc_symbol" not in reader.fieldnames:
            print(f"Expected `hgnc_symbol` column in {args.genes}", file=sys.stderr)
            return 1
        gene_symbols = {row["hgnc_symbol"].strip()
                        for row in reader if row.get("hgnc_symbol", "").strip()}
    print(f"Gene list: {len(gene_symbols)} symbols", file=sys.stderr)

    # Resolve gene → (chrom, start, end).
    print(f"Resolving coordinates from {args.gff} …", file=sys.stderr)
    intervals = _gene_intervals(args.gff, gene_symbols, args.buffer_bp)
    unresolved = gene_symbols - {sym for _, _, _, sym in intervals}
    if unresolved:
        print(f"WARN: unresolved gene symbols: {sorted(unresolved)}", file=sys.stderr)
    print(f"Resolved {len(intervals)} / {len(gene_symbols)} genes to intervals.",
          file=sys.stderr)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    bed = args.out_dir / "acmg_regions.bed"
    _write_bed(intervals, bed)
    print(f"Wrote {bed}", file=sys.stderr)

    merged = _merge_contig_intervals(intervals)
    chroms = sorted({c for c, _, _ in merged})
    print(f"Merged to {len(merged)} intervals over {len(chroms)} chromosomes.",
          file=sys.stderr)

    # For each chromosome, `bcftools view --regions-file` over the remote NYGC
    # VCF, write a per-chrom local subset, then concat at the end.
    per_chrom = []
    for chrom in chroms:
        url = _chrom_vcf_url(chrom).replace(NYGC_BASE, args.nygc_base)
        local = args.out_dir / f"acmg.{chrom}.vcf.gz"
        per_chrom.append(local)
        # `bcftools view -R <bed>` keeps only overlapping records.
        cmd = [
            "bcftools", "view",
            "-R", str(bed),
            "-O", "z",
            "--threads", str(args.threads),
            "-o", str(local),
            url,
        ]
        _run(cmd, args.dry_run)
        _run(["bcftools", "index", "--tbi", "--threads", str(args.threads), str(local)], args.dry_run)

    # Concat all per-chromosome subsets into one file.
    merged_vcf = args.out_dir / "acmg.phased.vcf.gz"
    if per_chrom:
        cmd = ["bcftools", "concat", "-O", "z", "--threads", str(args.threads),
               "-o", str(merged_vcf)] + [str(p) for p in per_chrom]
        _run(cmd, args.dry_run)
        _run(["bcftools", "index", "--tbi", "--threads", str(args.threads), str(merged_vcf)], args.dry_run)

    # Samples list (2504 unrelated + 698 related ≈ 3202 total)
    samples_out = args.out_dir / "samples.txt"
    if not args.dry_run and merged_vcf.exists():
        with samples_out.open("w") as f:
            subprocess.run(["bcftools", "query", "-l", str(merged_vcf)],
                           check=True, stdout=f)
        n_samples = sum(1 for _ in samples_out.open())
        print(f"Wrote {samples_out} ({n_samples} samples)", file=sys.stderr)

    print("Done.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
