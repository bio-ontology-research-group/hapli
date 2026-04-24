# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Hapli is a haplotype-level gene-function toolkit. Rather than analyzing variants in isolation, it materialises each of a diploid sample's two haplotype sequences for every gene (via `bcftools consensus` in Mode A, or pre-assembled FASTA in Mode B), lifts the reference GFF3 onto each with Liftoff, extracts + translates per-transcript proteins, diffs them against reference, and — for haplotypes with ≥2 coding substitutions — computes an **additive-vs-joint ESM2 epistasis residual**. Population-scale aggregation yields per-gene LoF allele frequencies + compound-het-LoF counts.

Managed with `uv`; Python >= 3.10.16.

## Commands

```bash
# Install deps (+ ML extras for --with-esm)
uv sync --extra ml

# Mode A: phased VCF → per-gene analysis JSON
uv run hapli analyze --gene GENE --sample S1 \
    --vcf variants.vcf.gz --reference ref.fa --gff annotations.gff3 \
    --output-dir results/ \
    [--with-esm --esm-model esm2_t6_8M_UR50D] \
    [--alphamissense-table am.tsv.gz --clinvar-vcf clinvar.vcf.gz] \
    [--gnomad-constraint gnomad.tsv --clingen-dosage clingen.tsv]

# Mode B: two pre-assembled haplotype FASTAs (HPRC, hifiasm, Verkko)
uv run hapli assess --gene GENE --sample HG002 \
    --hap1 HG002.hap1.fa --hap2 HG002.hap2.fa \
    --reference ref.fa --gff annotations.gff3 --output-dir results/

# Interactive TUI + LLM interpretation (needs OPENROUTER_API_KEY)
uv run hapli explore results/<sample>_<gene>_analysis.json
uv run hapli interpret results/<sample>_<gene>_analysis.json

# Population aggregation across many analysis JSONs
uv run hapli aggregate 'results/*_*_analysis.json' \
    --per-sample per_sample.tsv --per-gene per_gene.tsv

# Batch via Snakemake
uv run snakemake -s workflows/Snakefile --configfile workflows/config.yaml --cores <N>

# Tests
uv run pytest                     # testpaths=["tests"] in pyproject.toml; 171 pass
uv run pytest tests/test_x.py::test_name   # single test

# Regression-tested walkthrough (exits 0 if pipeline output matches recorded)
uvx showboat verify docs/walkthrough.md
```

External runtime dependencies on `$PATH`: `bcftools`, `tabix`, `samtools`, `liftoff`, `minimap2`. The `liftoff` CLI should be installed via `uv tool install git+https://github.com/agshumate/Liftoff` (the PyPI wheel pins pysam==0.16 which conflicts with hapli's pysam>=0.23). VCFs must be bgzipped + tabix-indexed.

## Architecture

The pipeline is orchestrated by `hapli/workflow/pipeline.py::HapliPipeline.run_gene_analysis`, which glues together four stages. Understanding this flow is the fastest way to get productive:

1. **GFF loading** (`hapli/core/io.py::GFFProcessor`) — does NOT use `gffutils.create_db`. Instead it does a targeted single-pass scan of the raw GFF3 keyed on the `--gene` argument, collecting only the gene and its descendants (via `Parent=` attributes and a 100 kb co-location window on the same seqid). Assumes the GFF is sorted and grouped by chromosome; it `break`s once it passes the gene region. A lightweight `Feature` class in the same file mimics the `gffutils.Feature` API surface the rest of the code uses.

2. **Haplotype construction** (`hapli/genotype/builder.py::HaplotypeBuilder`) — for the gene region (with a 1 kb buffer added by the pipeline), fetches reference bases via `pysam.FastaFile` and variants via `pysam.VariantFile`, then walks variants in order applying the sample's allele at `hap_idx` 0 or 1. Returns `{'hap1': seq, 'hap2': seq}`. **Coordinate convention is load-bearing**: region strings are 1-based inclusive (`chr:start-end`), FASTA fetch converts to 0-based half-open, and `_apply_variants` tracks both `current_idx` (0-based in the region slice) and `current_genomic_pos` (1-based) in parallel. Overlapping / out-of-order variants after a prior insertion/deletion are skipped with a warning. Unphased VCFs are treated as if `GT` order were hap1/hap2 (noted as unsafe in comments).

3. **Transcript extraction + alignment** (`hapli/alignment/aligner.py::SequenceAligner`) — the pipeline builds each mRNA's spliced cDNA by fetching exon slices on the + strand in genomic order and reverse-complementing the concatenation if `strand == '-'`. **Do not** use `SequenceExtractor.get_sequence` for per-exon splicing — it reverse-complements each exon individually, which is wrong when concatenating; the pipeline deliberately bypasses it via `seq_extractor.fasta.fetch` (see the inline comment in `pipeline.py`). The aligner shells out to `minimap2 -a -x splice --cs`, parses SAM through `pysam`, and the pipeline records identity = 1 − NM/len, CIGAR, and an `is_perfect` flag per haplotype.

4. **Output** — writes `<sample>_<gene>_haplotypes.fa` and `<sample>_<gene>_analysis.json` into `--output-dir`. The JSON shape (`gene`, `sample`, `region`, `transcripts[].alignments.{hap1,hap2}`) is the contract consumed by both `hapli/cli/tui.py` (Textual UI) and `hapli/interpretation/llm.py` (OpenRouter prompt builder). Changing that schema requires updating both consumers.

CLI entry point is `main.py` → `hapli/cli/main.py` (Typer app with `analyze` / `explore` / `interpret` / `version` subcommands). The Snakefile in `workflows/` is a thin wrapper that shells out to `uv run main.py analyze` per (sample, gene) pair from `config.yaml`.

## Repo state notes

- `git status` shows a lot of deletions and renames — the project was recently refactored from a flatter layout (single-file `cli.py`, `tui.py`, etc.) into the current `hapli/{cli,core,genotype,alignment,workflow,interpretation}/` package layout (see commit `791cdeb`). The old files have been moved under `obsolete/`. Treat `obsolete/` as historical reference, not live code.
- There is a `GEMINI.md` describing an older version of the architecture (with a `hapli/variation/` module and a `generate` / `align` CLI split). The current CLI is a single `analyze` command — prefer this file over `GEMINI.md` when they disagree.
- The repo root contains many large intermediate artifacts (`*.paf`, `*.gfa`, `*.fa`, `*.json`) from prior experimental runs; these are covered by `.gitignore` but present on disk. `hrnr.gfa` in particular is ~78 GB — don't `cat`/`grep` it.
