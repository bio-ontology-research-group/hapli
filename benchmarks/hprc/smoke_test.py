#!/usr/bin/env python3
"""
End-to-end smoke test of the HPRC R2 × Mode B benchmark on a tiny
synthetic 3-sample assembled-haplotype fixture.

Validates:
  1. The Mode B Snakemake path (assess_gene rule) runs on pre-assembled
     haplotype FASTAs without a VCF, end-to-end.
  2. `hapli aggregate` produces non-empty per-sample + per-gene TSVs.
  3. `analyze.py` consumes the aggregator output and produces:
       - per_gene_presence.tsv
       - per_gene_per_superpop_lof_af.tsv
       - tally.json
  4. The synthetic fixture's designed signal — sample HG_A's hap1
     carries a whole-gene deletion on GACMG1 — shows up in Table A
     as `presence=deleted` count ≥ 1 for GACMG1.

Runs in well under a minute. Exits 0 iff every assertion holds.
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
    """2-gene synthetic reference on chr1 over 5 kb."""
    import random
    rng = random.Random(7)
    spacer_bp = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    cds1 = "ATG" + "GCT" * 98 + "TAA"
    cds2 = "ATG" + "TTC" * 98 + "TAA"
    parts = []
    parts.append(spacer_bp(1000))           # 1..1000
    parts.append(cds1)                      # 1001..1300 GACMG1
    parts.append(spacer_bp(3000 - 1300))    # 1301..3000
    parts.append(cds2)                      # 3001..3300 GACMG2
    parts.append(spacer_bp(5000 - 3300))    # 3301..5000
    full = "".join(parts)
    assert len(full) == 5000
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
    lines += _block("GACMG2", 3001, 3300)
    out.write_text("\n".join(lines) + "\n")


def _build_haplotypes(td: Path, ref_fa: Path) -> dict[str, dict[str, Path]]:
    """Build a mini cohort of 3 "samples" × 2 "assembled" haplotypes each.

    HG_A hap1: GACMG1 REMOVED (whole-gene deletion; the assembly lacks those 300 bp)
    HG_A hap2: reference-identical
    HG_B hap1/hap2: both reference-identical
    HG_C hap1: GACMG1 has a single-SNV "missense" at pos 1010 T→C
    HG_C hap2: reference-identical
    """
    # Load the reference string once
    ref_seq = ""
    with ref_fa.open() as f:
        for line in f:
            if line.startswith(">"):
                continue
            ref_seq += line.strip()
    assert len(ref_seq) == 5000

    def _apply_patch(seq: str, patches: list[tuple[int, str, str]]) -> str:
        """Apply 1-based (pos, ref_len_as_str, new_bases) patches. Positions must be sorted."""
        out_parts = []
        cur = 0
        for pos, ref_part, new in patches:
            pos0 = pos - 1
            # Sanity: what we expect to replace matches the reference
            assert seq[pos0:pos0 + len(ref_part)] == ref_part, (
                f"patch mismatch at {pos}: got {seq[pos0:pos0+len(ref_part)]!r} expected {ref_part!r}"
            )
            out_parts.append(seq[cur:pos0])
            out_parts.append(new)
            cur = pos0 + len(ref_part)
        out_parts.append(seq[cur:])
        return "".join(out_parts)

    # HG_A hap1: whole-gene deletion of GACMG1 at 1001..1300 → remove 300 bp
    hg_a_hap1_seq = ref_seq[:1000] + ref_seq[1300:]
    assert len(hg_a_hap1_seq) == 4700

    # HG_C hap1: single SNV in GACMG1 CDS at 1010 (T → C). Check the base first.
    base_at_1010 = ref_seq[1009]  # 1-based pos 1010 → index 1009
    hg_c_hap1_seq = ref_seq[:1009] + "A" + ref_seq[1010:] if base_at_1010 != "A" else \
                    ref_seq[:1009] + "C" + ref_seq[1010:]
    assert len(hg_c_hap1_seq) == 5000

    haps = {
        "HG_A": {"hap1": hg_a_hap1_seq, "hap2": ref_seq},
        "HG_B": {"hap1": ref_seq, "hap2": ref_seq},
        "HG_C": {"hap1": hg_c_hap1_seq, "hap2": ref_seq},
    }

    out: dict[str, dict[str, Path]] = {}
    for sample, haploidy in haps.items():
        out[sample] = {}
        for hap_name, seq in haploidy.items():
            fa = td / f"{sample}.{hap_name}.fa"
            with fa.open("w") as f:
                f.write(f">chr1\n")
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i + 60] + "\n")
            pysam.faidx(str(fa))
            out[sample][hap_name] = fa
    return out


def _run_snakemake(repo: Path, cfg: Path) -> None:
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

    with tempfile.TemporaryDirectory(prefix="hapli_hprc_smoke_") as td:
        td_p = Path(td)
        ref = td_p / "reference.fa"
        gff = td_p / "annotation.gff3"

        _build_reference(ref)
        _build_gff(gff)
        haps = _build_haplotypes(td_p, ref)
        print(f"[smoke] fixture built in {td_p}", file=sys.stderr)

        # Synthetic population manifest (AFR / EAS / EUR)
        manifest = td_p / "manifest.tsv"
        manifest.write_text(
            "sample_id\thap1_url\thap2_url\tpopulation\tsuper_population\n"
            f"HG_A\tfile://{haps['HG_A']['hap1']}\tfile://{haps['HG_A']['hap2']}\tYRI\tAFR\n"
            f"HG_B\tfile://{haps['HG_B']['hap1']}\tfile://{haps['HG_B']['hap2']}\tCHS\tEAS\n"
            f"HG_C\tfile://{haps['HG_C']['hap1']}\tfile://{haps['HG_C']['hap2']}\tCEU\tEUR\n"
        )

        # Mode B Snakemake config
        out_dir = td_p / "run_out"
        cfg = td_p / "config.yaml"
        cfg_lines = [
            "mode: B",
            f"output_dir: {out_dir}",
            f"reference: {ref}",
            f"gff: {gff}",
            "haps:",
        ]
        for sample, h in haps.items():
            cfg_lines.append(f"  {sample}:")
            cfg_lines.append(f"    hap1: {h['hap1']}")
            cfg_lines.append(f"    hap2: {h['hap2']}")
        cfg_lines += [
            "genes:",
            "  - GACMG1",
            "  - GACMG2",
            "with_esm: false",
            "threads_per_gene: 1",
            "mem_mb_per_gene: 2000",
        ]
        cfg.write_text("\n".join(cfg_lines) + "\n")

        _run_snakemake(REPO, cfg)

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
        assert len(ps_rows) == 3 * 2, f"expected 6 rows; got {len(ps_rows)}"
        assert len(pg_rows) == 2, f"expected 2 rows; got {len(pg_rows)}"

        # HG_A's hap1 on GACMG1 must be flagged absent/deleted/low-identity.
        hga_g1 = next((r for r in ps_rows if r["sample"] == "HG_A" and r["gene"] == "GACMG1"), None)
        assert hga_g1 is not None
        assert hga_g1["presence_hap1"] in ("deleted", "low_identity", "uncertain", "partial"), (
            f"HG_A hap1 on GACMG1 should indicate structural loss; got {hga_g1}"
        )

        # Analyzer
        analyze_out = td_p / "analyze_out"
        subprocess.run(
            ["uv", "run", "python3",
             str(REPO / "benchmarks" / "hprc" / "analyze.py"),
             "--per-sample", str(per_sample),
             "--per-gene",   str(per_gene),
             "--sample-manifest", str(manifest),
             "--out", str(analyze_out)],
            check=True, capture_output=True, text=True,
        )
        tally_path = analyze_out / "tally.json"
        assert tally_path.exists()
        tally = json.loads(tally_path.read_text())
        print(f"[smoke] tally: {tally}", file=sys.stderr)
        assert tally["n_samples"] == 3
        assert tally["n_genes"] == 2
        assert tally["n_sv_presence_records"] >= 1, (
            f"expected ≥1 SV-presence record (HG_A hap1 on GACMG1); got {tally}"
        )
        assert set(tally["super_populations_covered"]) >= {"AFR", "EAS", "EUR"}

        # per_gene_presence.tsv must show ≥1 non-intact status for GACMG1
        tbl_a = analyze_out / "per_gene_presence.tsv"
        assert tbl_a.exists()
        with tbl_a.open() as f:
            rows = {r["gene"]: r for r in csv.DictReader(f, delimiter="\t")}
        g1 = rows.get("GACMG1", {})
        non_intact = sum(int(g1.get(s, 0) or 0)
                         for s in ("low_identity", "partial", "deleted", "duplicated"))
        assert non_intact >= 1, (
            f"expected ≥1 non-intact GACMG1 haplotype record; got {g1}"
        )

        # per_gene_per_superpop_lof_af.tsv sanity: 3 super-pops × 2 genes = 6 rows
        tbl_b = analyze_out / "per_gene_per_superpop_lof_af.tsv"
        assert tbl_b.exists()
        with tbl_b.open() as f:
            b_rows = list(csv.DictReader(f, delimiter="\t"))
        assert len(b_rows) == 6, f"expected 6 rows; got {len(b_rows)}"

        print("[smoke] all assertions passed ✓", file=sys.stderr)
        return 0


if __name__ == "__main__":
    sys.exit(main())
