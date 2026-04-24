#!/usr/bin/env python3
"""
End-to-end smoke test of the 1000G + ACMG-SF-v3.2 benchmark on a
synthetic multi-sample fixture.

Validates:
  1. fetch.py-equivalent output (reference + GFF + phased VCF) structure.
  2. Snakemake Mode A runs across (5 samples × 3 genes) without errors.
  3. `hapli aggregate` produces non-empty per-sample + per-gene TSVs.
  4. analyze.py consumes the aggregator output and emits a per-gene
     summary with at least one compound-het-LoF event and one
     SV-presence event, mirroring the paper's Axis 2 + Axis 1 signals.

Runs in under a minute. Exits 0 if every assertion holds; exits 1 with
a human-readable error otherwise.

Cohort design:
  * 5 samples total: S1, S2, S3, S4, S5 on a fixed-seed pedigree
  * 3 genes: GACMG1 (AR, compound-het-LoF carrier), GACMG2 (AD,
    intact), GACMG3 (AD, <DEL> on hap1 in one sample)
  * Expected signals:
      - compound_het_lof=True:  S1 on GACMG1
      - presence=deleted:        S3 on GACMG3 hap1
      - everything else intact.
"""
from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam


REPO = Path(__file__).resolve().parent.parent.parent


def _build_reference(out: Path) -> None:
    """3-gene synthetic reference on chr1 over 10 kb.

    Each gene's CDS is deliberately distinguishable so minimap2 (and
    thus Liftoff) can disambiguate when lifting a per-gene haplotype
    region back to the reference. Also uses a varied (not all-A) spacer
    so minimap2 can find seeds on flank regions.
    """
    import random
    rng = random.Random(42)
    spacer_bp = lambda n: "".join(rng.choice("ACGT") for _ in range(n))

    def _gene_cds(ala_codons: int, filler_codon: str) -> str:
        # ATG + filler*ala_codons + TAA
        return "ATG" + filler_codon * ala_codons + "TAA"

    # Each gene gets a distinguishing filler codon so the reference
    # CDS of GACMG{1,2,3} are not mutually substring-identical.
    cds1 = _gene_cds(98, "GCT")   # Ala-run
    cds2 = _gene_cds(98, "TTC")   # Phe-run (different 1st/2nd base)
    cds3 = _gene_cds(98, "CAA")   # Gln-run
    seq_parts = []
    seq_parts.append(spacer_bp(1000))           # 1..1000
    seq_parts.append(cds1)                      # 1001..1300  GACMG1
    seq_parts.append(spacer_bp(4000 - 1300))    # 1301..4000
    seq_parts.append(cds2)                      # 4001..4300  GACMG2
    seq_parts.append(spacer_bp(7000 - 4300))
    seq_parts.append(cds3)                      # 7001..7300  GACMG3
    seq_parts.append(spacer_bp(10000 - 7300))
    full = "".join(seq_parts)
    assert len(full) == 10000
    with out.open("w") as f:
        f.write(">chr1\n")
        for i in range(0, len(full), 60):
            f.write(full[i:i + 60] + "\n")
    pysam.faidx(str(out))


def _build_gff(out: Path) -> None:
    def _block(gene: str, start: int, end: int) -> list[str]:
        return [
            f"chr1\tsynth\tgene\t{start}\t{end}\t.\t+\t.\tID={gene};Name={gene};biotype=protein_coding",
            f"chr1\tsynth\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene}.t1;Parent={gene}",
            f"chr1\tsynth\texon\t{start}\t{end}\t.\t+\t.\tID={gene}.t1.e1;Parent={gene}.t1",
            f"chr1\tsynth\tCDS\t{start}\t{end}\t.\t+\t0\tID={gene}.t1.c1;Parent={gene}.t1",
        ]
    lines = ["##gff-version 3"]
    lines += _block("GACMG1", 1001, 1300)
    lines += _block("GACMG2", 4001, 4300)
    lines += _block("GACMG3", 7001, 7300)
    out.write_text("\n".join(lines) + "\n")


def _build_vcf(out_uncompressed: Path, out_gz: Path) -> None:
    """5-sample synthetic VCF.

    S1 on GACMG1: compound-het LoF.
        chr1:1010 GCT>TAA stop_gained on hap1  (1|0)
        chr1:1148 G>GA    frameshift   on hap2  (0|1)
    S2 on GACMG1: synonymous on hap2, no LoF.
        chr1:1009 T>C     synonymous   0|1
    S3 on GACMG3: symbolic <DEL> removing the whole gene on hap1.
        chr1:7000 A>.<DEL> SVTYPE=DEL;END=7300 1|0
    S4, S5: reference (no variants) for their unique positions.
    """
    header = [
        "##fileformat=VCFv4.2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
        "##ALT=<ID=DEL,Description=\"Deletion\">",
        "##contig=<ID=chr1,length=10000>",
    ]
    samples = ["S1", "S2", "S3", "S4", "S5"]
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
                  "\t".join(samples))
    rows = [
        # Sorted by (chrom, pos) for tabix. S2 synonymous (lowest pos in GACMG1).
        ("chr1", 1009, ".", "T", "C", ".", "PASS", ".", "GT",
         ["0|0", "0|1", "0|0", "0|0", "0|0"]),
        # S1 compound-het: stop-gain hap1 (pos 1010), frameshift hap2 (pos 1148)
        ("chr1", 1010, ".", "GCT", "TAA", ".", "PASS", ".", "GT",
         ["1|0", "0|0", "0|0", "0|0", "0|0"]),
        ("chr1", 1148, ".", "G", "GA", ".", "PASS", ".", "GT",
         ["0|1", "0|0", "0|0", "0|0", "0|0"]),
        # S3 whole-gene <DEL> on GACMG3, hap1
        ("chr1", 7000, ".", "A", "<DEL>", ".", "PASS", "SVTYPE=DEL;END=7300", "GT",
         ["0|0", "0|0", "1|0", "0|0", "0|0"]),
    ]
    with out_uncompressed.open("w") as f:
        for h in header:
            f.write(h + "\n")
        for row in rows:
            chrom, pos, vid, ref, alt, qual, filt, info, fmt, gts = row
            f.write("\t".join([chrom, str(pos), vid, ref, alt, qual, filt, info, fmt] + gts) + "\n")
    pysam.tabix_compress(str(out_uncompressed), str(out_gz), force=True)
    pysam.tabix_index(str(out_gz), preset="vcf", force=True)


def _run_snakemake(workdir: Path, repo: Path, cfg: Path) -> None:
    r = subprocess.run(
        ["uv", "run", "snakemake",
         "-s", str(repo / "workflows" / "Snakefile"),
         "--configfile", str(cfg),
         "--cores", "2",
         "--forceall"],
        cwd=repo, capture_output=True, text=True, timeout=600,
    )
    if r.returncode != 0:
        print("snakemake failed:", file=sys.stderr)
        print("--- stdout ---", file=sys.stderr); print(r.stdout, file=sys.stderr)
        print("--- stderr ---", file=sys.stderr); print(r.stderr, file=sys.stderr)
        raise RuntimeError("smoke-test snakemake run failed")


def main() -> int:
    if shutil.which("liftoff") is None:
        print("SKIP: liftoff not on PATH", file=sys.stderr)
        return 0

    with tempfile.TemporaryDirectory(prefix="hapli_1000g_smoke_") as td:
        td_p = Path(td)
        ref = td_p / "reference.fa"
        gff = td_p / "annotation.gff3"
        vcf_pre = td_p / "phased.vcf"
        vcf = td_p / "phased.vcf.gz"

        _build_reference(ref)
        _build_gff(gff)
        _build_vcf(vcf_pre, vcf)
        print(f"[smoke] fixture built in {td_p}", file=sys.stderr)

        # Synthetic ACMG gene list (our 3 fake genes)
        gene_list = td_p / "acmg.tsv"
        gene_list.write_text(
            "hgnc_symbol\tdisease_category\tinheritance\tsource\n"
            "GACMG1\tsynth_recessive\tAR\tsmoke\n"
            "GACMG2\tsynth_dominant\tAD\tsmoke\n"
            "GACMG3\tsynth_sv\tAD\tsmoke\n"
        )

        # Snakemake config
        out_dir = td_p / "run_out"
        cfg = td_p / "config.yaml"
        cfg.write_text(
            "mode: A\n"
            f"output_dir: {out_dir}\n"
            f"reference: {ref}\n"
            f"gff: {gff}\n"
            f"vcf: {vcf}\n"
            "samples:\n"
            "  - S1\n  - S2\n  - S3\n  - S4\n  - S5\n"
            "genes:\n"
            "  - GACMG1\n  - GACMG2\n  - GACMG3\n"
            "with_esm: false\n"
            "threads_per_gene: 1\n"
            "mem_mb_per_gene: 2000\n"
        )

        # Run the workflow
        _run_snakemake(td_p, REPO, cfg)

        # Consume the aggregate output
        per_sample = out_dir / "aggregate" / "per_sample.tsv"
        per_gene   = out_dir / "aggregate" / "per_gene.tsv"
        assert per_sample.exists(), f"missing {per_sample}"
        assert per_gene.exists(),   f"missing {per_gene}"

        with per_sample.open() as f:
            ps_rows = list(csv.DictReader(f, delimiter="\t"))
        with per_gene.open() as f:
            pg_rows = list(csv.DictReader(f, delimiter="\t"))

        print(f"[smoke] per_sample rows: {len(ps_rows)}", file=sys.stderr)
        print(f"[smoke] per_gene rows:   {len(pg_rows)}", file=sys.stderr)
        assert len(ps_rows) == 5 * 3, f"expected 15 (sample,gene) rows; got {len(ps_rows)}"
        assert len(pg_rows) == 3,     f"expected 3 per-gene rows; got {len(pg_rows)}"

        # Axis-2 check: S1 on GACMG1 must have compound_het_lof=True
        s1_g1 = next((r for r in ps_rows if r["sample"] == "S1" and r["gene"] == "GACMG1"), None)
        assert s1_g1 is not None, "missing S1/GACMG1 row"
        assert s1_g1["compound_het_lof"] in ("1", "True", "true"), (
            f"expected S1 compound_het_lof on GACMG1; got row={s1_g1}"
        )

        # Axis-1 check: S3 on GACMG3 must have presence=deleted on hap1
        s3_g3 = next((r for r in ps_rows if r["sample"] == "S3" and r["gene"] == "GACMG3"), None)
        assert s3_g3 is not None, "missing S3/GACMG3 row"
        assert s3_g3["presence_hap1"] == "deleted", (
            f"expected presence_hap1=deleted for S3 on GACMG3; got {s3_g3}"
        )

        # Run the analyzer
        analyze_out = td_p / "analyze_out"
        r = subprocess.run(
            ["uv", "run", "python3",
             str(REPO / "benchmarks" / "1000g_acmg" / "analyze.py"),
             "--per-sample", str(per_sample),
             "--per-gene",   str(per_gene),
             "--gene-list",  str(gene_list),
             "--out",        str(analyze_out)],
            check=True, capture_output=True, text=True,
        )
        tally_path = analyze_out / "tally.json"
        assert tally_path.exists(), "analyze.py didn't emit tally.json"
        tally = json.loads(tally_path.read_text())
        print(f"[smoke] tally: {tally}", file=sys.stderr)
        assert tally["n_genes_axis2_disagreement"] >= 1, (
            f"expected ≥1 axis-2 gene; got {tally}"
        )
        # axis-1 presence=deleted shows up in the aggregator's per-gene
        # n_deleted_alleles column, not its n_sv_presence_hap_records one;
        # so we check the per-gene TSV directly for GACMG3.
        g3 = next((r for r in pg_rows if r["gene"] == "GACMG3"), None)
        assert g3 is not None
        assert int(g3.get("n_deleted_alleles", "0")) >= 1, (
            f"expected ≥1 deleted allele for GACMG3; got row={g3}"
        )

        print("[smoke] all assertions passed ✓", file=sys.stderr)
        return 0


if __name__ == "__main__":
    sys.exit(main())
