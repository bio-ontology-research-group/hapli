#!/usr/bin/env python3
"""
Fetch a MaveDB scoreset + its target protein sequence.

Uses the MaveDB v1 API:
  - /api/v1/score-sets/<urn>                — metadata + target-gene info
  - /api/v1/score-sets/<urn>/scores         — the scores CSV (one row per variant)

and UniProt for the reference protein sequence.

Writes two files into the output directory:
  <urn>.scores.csv     — raw MaveDB scores CSV
  <urn>.reference.fa   — UniProt reference protein FASTA

Usage:
  ./fetch.py <urn> [-o OUTDIR]
  ./fetch.py urn:mavedb:00000081-a-1 -o data/brca1
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


MAVEDB = "https://api.mavedb.org/api/v1"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"


def _get(url: str, *, retries: int = 3, timeout: int = 60) -> bytes:
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, urllib.error.HTTPError) as exc:
            last_err = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"GET {url} failed after {retries} retries: {last_err}")


def fetch_scoreset_metadata(urn: str) -> dict:
    raw = _get(f"{MAVEDB}/score-sets/{urn}")
    return json.loads(raw)


def fetch_scores_csv(urn: str) -> bytes:
    return _get(f"{MAVEDB}/score-sets/{urn}/scores")


def uniprot_fasta(accession: str) -> str:
    raw = _get(UNIPROT_FASTA.format(acc=accession)).decode()
    return raw


def pick_uniprot_from_metadata(meta: dict) -> str | None:
    """Look through target-gene annotations for a UniProt accession."""
    for tg in meta.get("targetGenes", []):
        for ext in tg.get("externalIdentifiers", []):
            ident = ext.get("identifier") or {}
            if ident.get("dbName", "").lower() == "uniprot":
                return ident.get("identifier")
    return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("urn", help="MaveDB scoreset URN, e.g. urn:mavedb:00000081-a-1")
    ap.add_argument("-o", "--out", type=Path, default=Path("."), help="output directory")
    ap.add_argument("--uniprot", default=None,
                    help="override UniProt accession (skip metadata discovery)")
    args = ap.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    safe_name = args.urn.replace(":", "_").replace("/", "_")

    print(f"Fetching metadata for {args.urn} …", file=sys.stderr)
    meta = fetch_scoreset_metadata(args.urn)
    (args.out / f"{safe_name}.meta.json").write_text(json.dumps(meta, indent=2))

    print(f"Fetching scores CSV …", file=sys.stderr)
    scores_csv = fetch_scores_csv(args.urn)
    scores_path = args.out / f"{safe_name}.scores.csv"
    scores_path.write_bytes(scores_csv)
    print(f"  wrote {scores_path} ({len(scores_csv)} bytes)", file=sys.stderr)

    uniprot_acc = args.uniprot or pick_uniprot_from_metadata(meta)
    if uniprot_acc:
        print(f"Fetching UniProt {uniprot_acc} …", file=sys.stderr)
        fa = uniprot_fasta(uniprot_acc)
        ref_path = args.out / f"{safe_name}.reference.fa"
        ref_path.write_text(fa)
        print(f"  wrote {ref_path}", file=sys.stderr)
    else:
        print("  no UniProt accession in metadata; pass --uniprot to fetch the protein",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
