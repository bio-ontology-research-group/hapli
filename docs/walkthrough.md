# Hapli current-pipeline walkthrough (pre-refactor baseline)

*2026-04-23T08:17:39Z by Showboat 0.6.1*
<!-- showboat-id: 452615dc-c3f1-4fa0-8074-dedb19e66882 -->

This document is the **pre-refactor baseline** for the hapli pipeline. It captures every externally-visible behavior of the current `main.py analyze` command on the tiny test fixture in `data/test/`. Every code block is executed by showboat and its output captured; `showboat verify docs/walkthrough.md` must pass before and after every subsequent change during the Phase 1+ refactor. If this document drifts, either a bug has been fixed (re-record intentionally) or a regression has been introduced (fix).

## Environment

First, record the tool versions the baseline depends on, so 'verify' on a different machine can diagnose version drift rather than blaming hapli.

```bash
uv --version && python3 --version && minimap2 --version && bcftools --version | head -2
```

```output
uv 0.9.5
Python 3.13.5
2.27-r1193
bcftools 1.21
Using htslib 1.21
```

## Test fixtures

`data/test/` contains a tiny synthetic fixture emitted by `scripts/generate_comprehensive_test_data.py`: one chr1 + one chr2 reference, one gene annotation, a 3-variant phased VCF. This is the smallest case the end-to-end `analyze` command accepts.

```bash
ls -1 data/test/ | grep -v '_db$' && echo '---' && head -20 data/test/annotation.gff3
```

```output
annotation.gff3
reference.fa
reference.fa.fai
variants.vcf.gz
variants.vcf.gz.tbi
---
##gff-version 3
chr1	Sim	gene	10606	16344	.	+	.	ID=GENE00001;Name=Gene_1;biotype=protein_coding
chr1	Sim	mRNA	10655	14044	.	+	.	ID=GENE00001.t1;Parent=GENE00001
chr1	Sim	exon	10655	11059	.	+	.	ID=GENE00001.t1.exon1;Parent=GENE00001.t1
chr1	Sim	CDS	10655	11059	.	+	0	ID=GENE00001.t1.cds1;Parent=GENE00001.t1
chr1	Sim	exon	11210	11668	.	+	.	ID=GENE00001.t1.exon2;Parent=GENE00001.t1
chr1	Sim	CDS	11210	11668	.	+	0	ID=GENE00001.t1.cds2;Parent=GENE00001.t1
chr1	Sim	exon	12041	12489	.	+	.	ID=GENE00001.t1.exon3;Parent=GENE00001.t1
chr1	Sim	CDS	12041	12489	.	+	0	ID=GENE00001.t1.cds3;Parent=GENE00001.t1
chr1	Sim	exon	12811	13098	.	+	.	ID=GENE00001.t1.exon4;Parent=GENE00001.t1
chr1	Sim	CDS	12811	13098	.	+	0	ID=GENE00001.t1.cds4;Parent=GENE00001.t1
chr1	Sim	exon	13245	13639	.	+	.	ID=GENE00001.t1.exon5;Parent=GENE00001.t1
chr1	Sim	CDS	13245	13639	.	+	0	ID=GENE00001.t1.cds5;Parent=GENE00001.t1
chr1	Sim	exon	13776	14044	.	+	.	ID=GENE00001.t1.exon6;Parent=GENE00001.t1
chr1	Sim	CDS	13776	14044	.	+	0	ID=GENE00001.t1.cds6;Parent=GENE00001.t1
chr1	Sim	mRNA	10620	13081	.	+	.	ID=GENE00001.t2;Parent=GENE00001
chr1	Sim	exon	10620	11030	.	+	.	ID=GENE00001.t2.exon1;Parent=GENE00001.t2
chr1	Sim	CDS	10620	11030	.	+	0	ID=GENE00001.t2.cds1;Parent=GENE00001.t2
chr1	Sim	exon	11457	11865	.	+	.	ID=GENE00001.t2.exon2;Parent=GENE00001.t2
chr1	Sim	CDS	11457	11865	.	+	0	ID=GENE00001.t2.cds2;Parent=GENE00001.t2
```

```bash
bcftools view data/test/variants.vcf.gz 2>/dev/null | grep -v '^##bcftools_view' | tail -20
```

```output
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20251225
##reference=file:///home/leechuck/Public/software/hapli/data/test/reference.fa
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=50000>
##contig=<ID=chr2,length=10000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_001
chr1	10857	a1053fd6-ee78-4c7e-8832-04e5c3162e0d	T	A	99	PASS	SVTYPE=SNV	GT	1|0
chr1	11134	7c71bca7-302e-4913-a3a5-11208abb9d27	GGAGT	G	99	PASS	SVTYPE=Deletion	GT	0|1
chr1	13676	c7faf884-0453-4970-b46c-990882f565f4	CGTGCAAGCCCAAAACCTTTAGTACCGGGTTCTATGGTCAGAATGATGCATTCCCATCTGAAGTACTCATCCCACTATATTCTCTCTTGACTTACCACTGACACTTCTGTCCACCGGCCCCACTAGCGAATGGTCTCACGCAATAAACGACCTGGCGTAATGTGTATGATTGCTTTGAGTACCCTACGCATTATTGTAAAAATCTTAGTTTGGCCCTTTGGCAAACATCAGTGAATAATCGGAAGTCGTTTGCTCAGGAGGTGTCATTCCTGGCTTGCCGGTCGATCATGAATTCTCAAGAGACTAAATGGCCCAACCCTGCCAGATGGTTCTATCCGCCTGCGCAGATATCCATCACCGAAAAGTGCGTAATCCGGTAAACTGACCTGCCAGATAATTACGAATTGAACTACGAAAAATTTAAGAGCGTTCATTGACTGCTACTGAGGTCATATATCTTTATCCGGC	C	99	PASS	SVTYPE=Deletion	GT	1|1
```

## CLI help

The `hapli` CLI is a Typer app with four subcommands: `analyze`, `explore`, `interpret`, `version`.

```bash
uv run main.py --help 2>&1
```

```output
                                                                                
 Usage: main.py [OPTIONS] COMMAND [ARGS]...                                     
                                                                                
 Hapli: Genotype-centric variant analysis.                                      
                                                                                
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --install-completion          Install completion for the current shell.      │
│ --show-completion             Show completion for the current shell, to copy │
│                               it or customize the installation.              │
│ --help                        Show this message and exit.                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ───────────────────────────────────────────────────────────────────╮
│ analyze     Analyze a specific gene for a sample, generating haplotypes and  │
│             checking feature preservation.                                   │
│ assess      Mode B: assess a gene on pre-assembled haplotype FASTAs (HPRC,   │
│             hifiasm, Verkko).                                                │
│ explore     Explore alignment results interactively.                         │
│ interpret   Generate LLM-based functional interpretation.                    │
│ aggregate   Aggregate many `analyze` / `assess` outputs into                 │
│             population-level TSVs.                                           │
│ version     Show version.                                                    │
╰──────────────────────────────────────────────────────────────────────────────╯

```

```bash
uv run main.py analyze --help 2>&1
```

```output
                                                                                
 Usage: main.py analyze [OPTIONS]                                               
                                                                                
 Analyze a specific gene for a sample, generating haplotypes and checking       
 feature preservation.                                                          
                                                                                
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --gene                                   TEXT  Target gene name/ID        │
│                                                   [required]                 │
│ *  --sample                                 TEXT  Sample name in VCF         │
│                                                   [required]                 │
│ *  --vcf                                    PATH  Path to phased VCF         │
│                                                   [required]                 │
│ *  --reference                              PATH  Path to reference FASTA    │
│                                                   [required]                 │
│ *  --gff                                    PATH  Path to GFF3 annotation    │
│                                                   [required]                 │
│ *  --output-dir                             PATH  Directory for results      │
│                                                   [required]                 │
│    --with-esm                                     Compute the ESM2 epistasis │
│                                                   residual (requires the     │
│                                                   `ml` extra).               │
│    --esm-model                              TEXT  fair-esm checkpoint name   │
│                                                   (ignored unless            │
│                                                   --with-esm).               │
│                                                   [default:                  │
│                                                   esm2_t6_8M_UR50D]          │
│    --gnomad-constraint                      PATH  gnomAD constraint TSV      │
│                                                   (canonical per-gene pLI /  │
│                                                   MisZ / o/e).               │
│    --clingen-dosage                         PATH  ClinGen dosage sensitivity │
│                                                   TSV.                       │
│    --alphamissense-table                    PATH  bgzipped+tabix-indexed     │
│                                                   AlphaMissense TSV          │
│                                                   (per-variant lookup).      │
│    --clinvar-vcf                            PATH  ClinVar VCF                │
│                                                   (bgzipped+tabix-indexed) — │
│                                                   annotates consequences     │
│                                                   with clnsig.               │
│    --verbose                --no-verbose          [default: no-verbose]      │
│    --help                                         Show this message and      │
│                                                   exit.                      │
╰──────────────────────────────────────────────────────────────────────────────╯

```

## End-to-end run on the test fixture

Run `analyze` on the test fixture. This is the critical behavior under regression test: the JSON output schema, the protein-FASTA side artifact, and the logged events are all contractual surfaces for the TUI and LLM consumers.

```bash
rm -rf results/walkthrough && uv run main.py analyze --gene GENE00001 --sample sample_001 --vcf data/test/variants.vcf.gz --reference data/test/reference.fa --gff data/test/annotation.gff3 --output-dir results/walkthrough 2>&1 | sed -E 's/\[[0-9-]+ [0-9:]+\]/[ts]/'
```

```output
[ts] INFO: Starting analysis for gene GENE00001, sample sample_001
[ts] INFO: Searching for gene 'GENE00001'...
[ts] INFO: Found gene GENE00001 at chr1:10606-16344
[ts] INFO: Materialising haplotypes for chr1:9606-17344 sample=sample_001 via bcftools consensus
[ts] INFO: Liftoff hap1: gene GENE00001 status=low_identity
[ts] INFO: Liftoff hap2: gene GENE00001 status=low_identity
[ts] INFO: csq: 0 consequence records
[ts] INFO: hap1 GENE00001.t2: protein identity=0.9980 subs=1 pts=23 fs_region=no
[ts] INFO: hap1 GENE00001.t1: protein identity=0.8795 subs=1 pts=10 fs_region=no
[ts] INFO: hap2 GENE00001.t2: protein identity=1.0000 subs=0 pts=23 fs_region=no
[ts] INFO: hap2 GENE00001.t1: protein identity=0.8808 subs=0 pts=10 fs_region=no
[ts] INFO: GENE00001 diploid: hap1=0.000 hap2=0.000 min=0.000 max=0.000 compound_het_lof=True
[ts] INFO: Analysis saved to results/walkthrough/sample_001_GENE00001_analysis.json
```

```bash
ls -1 results/walkthrough/ | grep -vE 'liftoff_tmp|_polished|\.mmi$' && echo '---' && uv run python -c "
import json, sys
d = json.load(open('results/walkthrough/sample_001_GENE00001_analysis.json'))
# Emit v1 fields + a trimmed evidence block (per-haplotype presence status only)
# for a stable, reviewable snapshot. Full evidence is in the JSON on disk.
out = {k: d[k] for k in ('schema_version', 'gene', 'sample', 'region', 'transcripts')}
out['evidence'] = {
    'presence': {
        h: {k: d['evidence']['presence'][h][k]
            for k in ('status', 'source', 'coverage', 'sequence_identity', 'valid_orfs')}
        for h in ('hap1', 'hap2')
    },
    'n_consequence': len(d['evidence']['consequence']),
}
json.dump(out, sys.stdout, indent=2)
print()
"
```

```output
sample_001_GENE00001_analysis.json
sample_001_GENE00001.csq.vcf.gz
sample_001_GENE00001.csq.vcf.gz.tbi
sample_001_GENE00001_hap1.fa
sample_001_GENE00001_hap1.fa.fai
sample_001_GENE00001_hap1.lifted.gff3
sample_001_GENE00001_hap1.pep.fa
sample_001_GENE00001_hap1.unmapped.txt
sample_001_GENE00001_hap2.fa
sample_001_GENE00001_hap2.fa.fai
sample_001_GENE00001_hap2.lifted.gff3
sample_001_GENE00001_hap2.pep.fa
sample_001_GENE00001_hap2.unmapped.txt
sample_001_GENE00001_haplotypes.fa
sample_001_GENE00001_ref.pep.fa
---
{
  "schema_version": "2.0",
  "gene": "GENE00001",
  "sample": "sample_001",
  "region": "chr1:9606-17344",
  "transcripts": [
    {
      "id": "GENE00001.t2",
      "alignments": {
        "hap1": {
          "identity": 0.9888451443569554,
          "nm": 17,
          "cigar": "410M423N1M3D409M146N3D406M2D6M84N148M280N1M1I1M3I2M1D136M",
          "is_perfect": false
        },
        "hap2": {
          "identity": 0.989501312335958,
          "nm": 16,
          "cigar": "410M419N1M3D409M146N3D406M2D6M84N148M280N1M1I1M3I2M1D136M",
          "is_perfect": false
        }
      }
    },
    {
      "id": "GENE00001.t1",
      "alignments": {
        "hap1": {
          "identity": 0.9951434878587196,
          "nm": 11,
          "cigar": "403M147N3D461M368N4D448M2D1M319N289M146N395M268S",
          "is_perfect": false
        },
        "hap2": {
          "identity": 0.9955849889624724,
          "nm": 10,
          "cigar": "403M143N3D461M368N4D448M2D1M319N289M146N395M268S",
          "is_perfect": false
        }
      }
    }
  ],
  "evidence": {
    "presence": {
      "hap1": {
        "status": "low_identity",
        "source": "liftoff",
        "coverage": 0.929,
        "sequence_identity": 0.928,
        "valid_orfs": 0
      },
      "hap2": {
        "status": "low_identity",
        "source": "liftoff",
        "coverage": 0.929,
        "sequence_identity": 0.929,
        "valid_orfs": 0
      }
    },
    "n_consequence": 0
  }
}
```

```bash
head -4 results/walkthrough/sample_001_GENE00001_haplotypes.fa | cut -c1-80 && echo '...' && awk '/^>hap1/{h=1; next} /^>hap2/{exit} h{n+=length($0)} END{print "hap1 length:", n}' results/walkthrough/sample_001_GENE00001_haplotypes.fa
```

```output
>hap1
AGAAATGGTCAACAACTAAGACCTCCTGTGACGGCCTGAATGATAACATTGCTACTCATAGCGAGTTTCAGTATTCCTTG
>hap2
AGAAATGGTCAACAACTAAGACCTCCTGTGACGGCCTGAATGATAACATTGCTACTCATAGCGAGTTTCAGTATTCCTTG
...
hap1 length: 7272
```

### Baseline schema

The JSON above is **schema v2**. The v1 top-level keys (`gene`, `sample`, `region`, `transcripts[]`) are preserved byte-for-byte so the existing TUI (`hapli/cli/tui.py`) and LLM consumer (`hapli/interpretation/llm.py`) keep working. Phase 1 adds:

- `schema_version`: `"2.0"` marker.
- `evidence.presence.{hap1,hap2}`: Liftoff-derived per-haplotype presence call (status + coverage + sequence_identity + valid_orfs). On this fixture both haplotypes lift at ~93% coverage / ~93% identity (the homozygous 465 bp deletion of the last exon drags the numbers down), with `valid_orfs=0` because the truncated gene has no valid ORF — hence `status="low_identity"`.
- `evidence.consequence`: list of per-variant per-haplotype HGVS consequences from `bcftools csq -p a --force`. Empty here because the fixture GFF has CDS phases that csq can't fully reconcile even with `--force` (phase=0 throughout, whereas some cumulative CDS lengths aren't ÷3). Real-world GFFs (Ensembl, GENCODE) don't hit this.
- `evidence.splice`, `evidence.lof`, `evidence.missense_agg`, `evidence.protein`, `evidence.epistasis`, `evidence.diploid`: placeholders for Phase 2+ (SpliceAI, LOFTEE, AlphaMissense, ESM2, epistasis residual, constraint-annotated diploid report).

The companion `<sample>_<gene>_haplotypes.fa` contains both haplotype sequences for the gene region (reference coordinates ± 1 kb buffer). On this fixture hap1 is 7272 bp and hap2 is 7268 bp — both carry the 465 bp homozygous deletion at chr1:13676; hap2 additionally carries the 4 bp heterozygous deletion at chr1:11134; hap1 carries the SNV at chr1:10857. Phase 1 also writes per-haplotype single-record FASTAs (`*_hap1.fa`, `*_hap2.fa`) named as `>chr1` so Liftoff can consume them as target assemblies.

**Phase 0 bug fix recorded here**: the pre-refactor hand-rolled `_apply_variants` indexed `pysam.VariantRecordSample.alleles` (genotype-ordered) with `allele_indices` (VCF-ordered), silently swallowing the SNV on hap1 when phase was `1|0`. The SNV then appeared in *neither* haplotype, undercounting edit distance by 1 on every transcript alignment touching that locus. The new `bcftools consensus` path correctly applies it, visible as NM going from 10→11 on `GENE00001.t1/hap1` and 16→17 on `GENE00001.t2/hap1`.

## Other subcommands

### `version`

```bash
uv run main.py version 2>&1
```

```output
hapli 0.3.0
```

### `interpret`

Without an `OPENROUTER_API_KEY` set, `interpret` prints a graceful fallback instead of calling the LLM. Captured here to pin the no-API-key behavior.

```bash
unset OPENROUTER_API_KEY; uv run main.py interpret results/walkthrough/sample_001_GENE00001_analysis.json 2>&1
```

```output
No OpenRouter API key provided. Skipping LLM interpretation.
--- Interpretation Report ---
Interpretation skipped (No API Key).
-----------------------------
```

### `explore`

`uv run main.py explore results/walkthrough/sample_001_GENE00001_analysis.json` opens a Textual TUI for interactive inspection. **Not captured here** because showboat's exec blocks are non-interactive — running a TUI in the baseline would hang. Manual smoke-test instruction: launch the command, verify the transcript tree and per-haplotype alignment panels render, press q to quit. The Phase 1+ refactor must keep the TUI able to consume whatever JSON `analyze` emits.

## Verification

Re-run `uvx showboat verify docs/walkthrough.md` after any refactor. A non-zero exit signals that either:
1. A bug was intentionally fixed (re-record with `showboat exec` on the changed block), or
2. A regression was introduced (fix the regression, do not re-record).

The schema v1 contract (top-level `gene`, `sample`, `region`, `transcripts[].alignments.{hap1,hap2}.{identity,nm,cigar,is_perfect}`) is load-bearing for the TUI and LLM. Any Phase 1 addition must nest under `evidence` and leave v1 fields untouched.
