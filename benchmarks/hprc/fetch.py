#!/usr/bin/env python3
"""
Fetch HPRC Release 2 diploid assembly FASTAs for a subset of samples.

The HPRC R2 assemblies live in an AWS Open Data S3 bucket (no auth required):
    s3://human-pangenomics/working/HPRC/<SAMPLE>/assemblies/<release>/<SAMPLE>.<hap>.fa.gz

Each hap is ~600 MB gzipped, ~3 GB unzipped. At 20 samples the total
download is ~24 GB (unzipped); at the full 200+ samples, ~250 GB.

This script:
  1. Reads a sample manifest (TSV from `sample_list.tsv`).
  2. Downloads hap1 + hap2 for each requested sample into --out-dir.
     Uses `aws s3 cp --no-sign-request` if available (fast), else curl.
  3. Unzips + faidx's the downloaded FASTAs.
  4. Prints a YAML snippet for `workflows/config.yaml`'s `haps:` block.

Usage:
  # Fetch just 5 samples for a scale-down proof-of-concept:
  uv run python3 benchmarks/hprc/fetch.py \\
      --manifest benchmarks/hprc/sample_list.tsv \\
      --samples HG00438 HG00621 HG00673 HG00733 HG00735 \\
      --out-dir data/hprc/assemblies/

  # All samples in the manifest:
  uv run python3 benchmarks/hprc/fetch.py \\
      --manifest benchmarks/hprc/sample_list.tsv \\
      --out-dir data/hprc/assemblies/

  # Refresh the sample manifest from the canonical HPRC index
  # (requires internet + jq-compatible python):
  uv run python3 benchmarks/hprc/fetch.py --refresh-manifest \\
      --out-manifest benchmarks/hprc/sample_list.tsv
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


CANONICAL_MANIFEST_URL = (
    "https://raw.githubusercontent.com/human-pangenomics/hpp_pangenome_resources/"
    "main/release_2/assembly_index.txt"
)


def _read_manifest(path: Path) -> list[dict]:
    rows = []
    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            # Support either tab or whitespace
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if parts[0] == "sample_id":  # header
                continue
            row = {
                "sample_id": parts[0],
                "hap1_url": parts[1],
                "hap2_url": parts[2],
                "population": parts[3] if len(parts) > 3 else "",
                "super_population": parts[4] if len(parts) > 4 else "",
            }
            rows.append(row)
    return rows


def _download(url: str, dest: Path, use_aws_cli: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if use_aws_cli and url.startswith("s3://"):
        subprocess.run(
            ["aws", "s3", "cp", "--no-sign-request", url, str(dest)],
            check=True,
        )
        return
    # Convert s3:// to the HTTPS public endpoint for curl fallback.
    if url.startswith("s3://"):
        url = url.replace("s3://", "https://s3.amazonaws.com/", 1)
    subprocess.run(["curl", "-sSL", "-o", str(dest), url], check=True)


def _unzip_and_index(gz: Path) -> Path:
    """Return the uncompressed, faidx'd FASTA path."""
    fa = gz.with_suffix("")  # strip `.gz`
    if not fa.exists():
        subprocess.run(["gunzip", "-k", str(gz)], check=True)
    subprocess.run(["samtools", "faidx", str(fa)], check=True)
    return fa


def _refresh_manifest(out_manifest: Path) -> int:
    """Download the canonical HPRC release_2 assembly index and rewrite
    `sample_list.tsv` in our four-column layout. Best effort — the upstream
    format changes occasionally; fall back to leaving the current list in place
    if the parse fails.
    """
    try:
        print(f"Fetching {CANONICAL_MANIFEST_URL} …", file=sys.stderr)
        with urllib.request.urlopen(CANONICAL_MANIFEST_URL, timeout=30) as r:
            txt = r.read().decode("utf-8")
    except Exception as exc:
        print(f"Failed to fetch canonical manifest: {exc}", file=sys.stderr)
        return 1
    # Parse the canonical format. As of Release 2 the columns are roughly
    # sample<tab>haplotype<tab>path; group by sample, collect hap1/hap2 URLs.
    by_sample: dict[str, dict] = {}
    for line in txt.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        sample, hap, url = parts[0], parts[1], parts[2]
        d = by_sample.setdefault(sample, {"sample_id": sample,
                                          "hap1_url": "", "hap2_url": "",
                                          "population": "", "super_population": ""})
        key = f"hap{hap[-1]}_url" if hap[-1] in ("1", "2") else None
        if key:
            d[key] = url
    # Keep only samples with both haps
    complete = [d for d in by_sample.values() if d["hap1_url"] and d["hap2_url"]]
    print(f"Resolved {len(complete)} samples with both haps", file=sys.stderr)
    with out_manifest.open("w") as f:
        f.write("sample_id\thap1_url\thap2_url\tpopulation\tsuper_population\n")
        for d in complete:
            f.write("\t".join([d["sample_id"], d["hap1_url"], d["hap2_url"],
                               d["population"], d["super_population"]]) + "\n")
    print(f"Wrote {out_manifest}", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=Path("benchmarks/hprc/sample_list.tsv"))
    ap.add_argument("--samples", nargs="+", help="sample IDs to fetch (default: all)")
    ap.add_argument("--out-dir", type=Path, default=Path("data/hprc/assemblies"))
    ap.add_argument("--refresh-manifest", action="store_true",
                    help="Download the current HPRC R2 index and rewrite the manifest")
    ap.add_argument("--out-manifest", type=Path, default=None,
                    help="Target manifest when --refresh-manifest is used")
    ap.add_argument("--print-config", action="store_true",
                    help="After fetching, print a `haps:` YAML block for config.yaml")
    args = ap.parse_args(argv)

    if args.refresh_manifest:
        return _refresh_manifest(args.out_manifest or args.manifest)

    manifest = _read_manifest(args.manifest)
    print(f"Manifest has {len(manifest)} samples", file=sys.stderr)

    if args.samples:
        wanted = set(args.samples)
        manifest = [r for r in manifest if r["sample_id"] in wanted]
        missing = wanted - {r["sample_id"] for r in manifest}
        if missing:
            print(f"WARN: sample(s) not in manifest: {sorted(missing)}", file=sys.stderr)

    use_aws = shutil.which("aws") is not None
    if not use_aws:
        print("aws CLI not found; falling back to curl (slower)", file=sys.stderr)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    config_haps: list[tuple[str, Path, Path]] = []
    for row in manifest:
        sample = row["sample_id"]
        h1_gz = args.out_dir / f"{sample}.hap1.fa.gz"
        h2_gz = args.out_dir / f"{sample}.hap2.fa.gz"
        for gz, url in ((h1_gz, row["hap1_url"]), (h2_gz, row["hap2_url"])):
            if gz.exists() and gz.stat().st_size > 1_000_000:
                print(f"  [skip] {gz} already present", file=sys.stderr)
                continue
            print(f"  {sample}: downloading {url} -> {gz}", file=sys.stderr)
            try:
                _download(url, gz, use_aws)
            except subprocess.CalledProcessError as exc:
                print(f"  [fail] {sample} {gz.name}: {exc}", file=sys.stderr)
                continue
        try:
            h1 = _unzip_and_index(h1_gz)
            h2 = _unzip_and_index(h2_gz)
            config_haps.append((sample, h1, h2))
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            print(f"  [fail] index {sample}: {exc}", file=sys.stderr)

    print(f"\nFetched {len(config_haps)} samples, both haplotypes each.",
          file=sys.stderr)

    if args.print_config and config_haps:
        print("# Paste into config.yaml under `haps:`", file=sys.stderr)
        print("haps:")
        for sample, h1, h2 in config_haps:
            print(f"  {sample}:")
            print(f"    hap1: {h1}")
            print(f"    hap2: {h2}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
