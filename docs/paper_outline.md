# Hapli — paper outline (working draft)

> Working consolidation of Phase 0–5 work into a paperable structure. Numbers
> in **bold** are headline results already produced; numbers in *italic* are
> targets for runs not yet executed (see [pending] markers).

## Title (working)

**hapli: genotype-level gene-function assessment. Two aggregation axes
(variants → haplotype, haplotypes → genotype) that per-variant tools
cannot express.**

## Abstract

Clinical variant interpretation tools (VEP / LOFTEE, BCFtools/csq,
AlphaMissense, SpliceAI) emit one record per variant. Collapsing those
records into an actionable call for a diploid individual requires two
distinct aggregation steps that no per-variant tool performs itself:

1. **Variants within a haplotype → whole-gene function on that haplotype.**
   Per-variant tools produce *N* consequences; a clinician needs "is this
   gene preserved on allele X?". The common practice is a hand-rolled
   worst-consequence-wins rule (LOFTEE HC) or max-pathogenicity-over-
   variants (AlphaMissense aggregation).
2. **Two haplotypes → one diploid genotype.** Whether a gene is functional
   in the individual depends on both alleles: dominant ("one hit is enough")
   vs. recessive / compound-heterozygous ("need to hit both"). Per-variant
   tools emit nothing at this layer; it is implemented as downstream custom
   logic in every clinical pipeline.

We present **hapli**, which makes both aggregations first-class outputs.
For each gene on each haplotype, hapli materialises the mutated protein
sequence (via `bcftools consensus` + Liftoff + translation), diffs it
against the reference, and reports a set of per-haplotype scalars that
are qualitatively richer than any per-variant aggregation can reach:
identity, premature-stop position, a frameshift-region tag (including
a restored-frame flag), and a Liftoff-derived presence class (intact /
low_identity / partial / deleted / duplicated). At the genotype level,
the diploid report exposes `{hap1_score, hap2_score, min, max,
compound_het_lof}` plus the paired presence state.

On 12 hand-designed motivating cases (`data/paper_cases/`), hapli flags
five distinct disagreements with bcftools/csq that are inexpressible at
the per-variant layer:
**(i)** compound-het LoF (stop-gain × frameshift in trans),
**(ii)** frameshift+rescue (+1/−1 indel pair, restored-frame flag),
**(iii)** whole-gene deletion via symbolic `<DEL>`,
**(iv)** inversion bisecting a gene body (low_identity presence),
**(v)** tandem duplication (duplicated presence).
An optional ESM2 full-sequence score, offered as a possible supplement,
fires with large positive residuals on designed epistatic haplotypes
(+57 on frameshift-rescue, +12 on compound-missense-pocket) but does not
recover a generic epistasis signal on bulk pairwise DMS of a 56-aa
protein (464 408-double GB1 scan, AUROC essentially flat; reported as
a null result, §3.2). Hapli's contribution is at the aggregation layer,
not at the per-variant scoring layer.

Hapli runs end-to-end in 2–4 h per whole-genome sample on 16 cores + 1 GPU.

## 1. Introduction

Clinical variant interpretation for a diploid individual requires *two*
aggregation steps that existing tools do not perform themselves:

```
        per-variant records (VEP, LOFTEE, AlphaMissense, SpliceAI, bcftools/csq)
                                │
                   [Axis 1]  variants → haplotype gene-function
                                │
                 per-haplotype records
                                │
                   [Axis 2]  haplotypes → diploid genotype call
                                │
                      actionable call
```

### Axis 1 — variants within a haplotype

Per-variant tools emit *N* consequence records per gene per haplotype.
A clinician needs one call per haplotype: "is the protein preserved on
this allele?" The state of practice is an ad-hoc rule applied downstream:

- LOFTEE's worst-consequence cascade (HC LoF > LC LoF > missense > …).
- AlphaMissense max or count-above-threshold over variants in the gene.
- BCFtools/csq concatenates HGVS consequences along a phased haplotype
  but stops there — it does not emit a scalar "gene functional" call.

None of these are wrong, but they are lossy collapses of variant-level
records. They cannot express three classes of signal that are only
visible from the whole mutated protein:

1. **Frameshift+rescue** — +1 / −1 indel pair on the same haplotype.
   Each variant is individually frame-disrupting; jointly the frame is
   restored after a small substitution window and the protein is full-
   length. Danecek 2017 quantifies the phenomenon on 1000G (501 / 5019
   compound-call flips when csq switches on phasing). No per-variant
   aggregation rule recovers this: both variants are still "LoF by
   cascade".
2. **Gene-level structural variation** — symbolic `<DEL>`, `<INV>`,
   `<DUP:TANDEM>` events remove, bisect, or duplicate entire exons or
   genes. csq is silent on these records by design. A per-haplotype
   annotation lift (Liftoff) maps the reference GFF onto the materialised
   haplotype and classifies presence (intact / low_identity / deleted /
   duplicated). This information is *structurally* absent from per-
   variant outputs.
3. **Compound missense with joint structural effect** — two individually
   benign substitutions that together disrupt a binding pocket. Per-
   variant tools emit two "benign missense" records; clinical rules
   taking the max pathogenicity say "benign × benign = benign". Joint
   effects are a real signal (Olson 2014; PRESCOTT 2025) but small
   enough that detecting them on bulk pairwise DMS is hard (§3.2).

### Axis 2 — two haplotypes into a genotype

Whether a gene is functional in the individual depends on both alleles.
Dominant disorders require one hit; recessive / compound-heterozygous
disorders require a hit on each. No per-variant tool emits this
aggregation: it is implemented as downstream custom logic in every
clinical pipeline, on top of per-variant records. Hapli emits the
diploid report — per-haplotype score, min / max, and the specific
compound-het-LoF flag — as a first-class JSON field. We deliberately do
not commit to a single dominance model; clinical consumers apply their
own rule on top.

### What hapli adds

Hapli makes both aggregations first-class outputs, driven by the mutated
protein sequence on each haplotype. Per-haplotype scalars (identity,
premature-stop position, frameshift-region tag, Liftoff presence class)
handle axis 1; a diploid aggregator (`hap1_score`, `hap2_score`,
`compound_het_lof`) handles axis 2. An optional ESM2 full-sequence score
is offered on axis 1 as a supplementary epistasis signal; §3.2 reports
an honest null on that signal in bulk DMS.

Functional readouts from MAVE / DMS provide a non-circular ground truth
for evaluating per-variant scoring (MaveDB v2 has 200+ scoresets) —
distinct from the computationally-derived labels (ClinVar, gnomAD
constraint) on which most per-variant predictors were trained. For the
aggregation axes where hapli's contribution sits, the natural evaluation
is not an AUROC against DMS but a qualitative demonstration on motivating
cases; §3.3 does that for 12 curated cases.

## 2. Methods

### 2.1 Pipeline architecture

Two operating modes:

| Mode | Input | Use case |
|---|---|---|
| A (`hapli analyze`) | reference + phased VCF + GFF3 + sample | clinical short-read pipelines |
| B (`hapli assess`)  | reference + 2 pre-assembled haplotype FASTAs + GFF3 | HPRC, hifiasm, Verkko |

Per (sample, gene) the pipeline is organised by the two aggregation axes
of §1. Stages 1–6 each contribute an Axis 1 signal (variants → whole-gene
function on one haplotype); stage 7 performs the Axis 2 collapse (two
haplotypes → genotype).

1. **Haplotype materialisation** *(prerequisite)* — `bcftools consensus
   -H 1/2 -s SAMPLE` for Mode A; pass-through for Mode B. Symbolic
   structural-variant ALTs that bcftools (through v1.21) refuses —
   `<INV>` and `<DUP:TANDEM>` — are resolved upstream of bcftools by hapli's
   own preprocessor (`hapli/external/consensus.py::_resolve_symbolic_svs`):
   each symbolic record is rewritten into an explicit REF/ALT sequence pair
   (REF = reference slice over POS..END; ALT = reverse-complement for `<INV>`,
   REF+REF for `<DUP*>`). `<BND>`, `<CNV>`, and symbolic `<INS>` without
   sequence are rejected at the preflight stage.
2. **Liftoff** *(Axis 1 signal — presence)* — transfer the reference GFF3
   onto each haplotype FASTA; yields per-gene `presence` call
   (`intact / partial / low_identity / duplicated / deleted`). This is
   where gene-level SVs surface.
3. **bcftools csq** *(per-variant records; Mode A only)* — per-variant
   per-haplotype HGVS consequences. hapli consumes these as input to the
   Axis 1 aggregation below, not as final output. Includes optional
   ClinVar `clnsig` annotation per consequence (`--clinvar-vcf`).
4. **Protein extraction + diff** *(Axis 1 signals — identity / PTC /
   frameshift-rescue)* — slice each transcript's CDS from the lifted
   haplotype GFF3, translate, align against the reference protein
   (`Bio.Align.PairwiseAligner`), emit `identity`, `premature_stop_at`,
   `frameshift_region.restored` flag. The mutated protein FASTA is written
   as a first-class artifact.
5. **AlphaMissense aggregation** *(Axis 1 signal, opt-in)* — per-
   (transcript, hap) max / mean / sum AM scores over the coding variants
   actually applied.
6. **ESM2 full-sequence score** *(Axis 1 supplement, opt-in)* — for any
   (transcript, hap) with ≥ 2 equal-length coding substitutions:
   `S_joint = pLL(hap_protein) − pLL(ref_protein)` and
   `S_additive = Σ log p(alt|ref) − log p(ref|ref)`, `residual = S_joint − S_additive`.
   Calibration / null on bulk DMS, fires on designed cases (§3.3 / §3.4).
7. **Diploid aggregator** *(Axis 2 collapse)* — takes per-haplotype Axis 1
   outputs and emits `{hap1_score, hap2_score, min, max, compound_het_lof,
   constraints: {pLI, MisZ, ClinGen}}`. Deliberately does **not** collapse
   to a single diploid call; clinical consumers apply their own dominance
   rule on top. `compound_het_lof` is the canonical Axis 2 signal.

Output JSON schema v2 in `hapli/core/schema.py`. Population aggregation
across samples via `hapli aggregate` → per-sample TSV + per-gene allele-frequency TSV.

### 2.2 Datasets

- **Curated motivating cases** (`scripts/generate_paper_cases.py`) — 12
  hand-designed cases covering the full LoF/SV matrix:
  01 synonymous pair, 02 single benign missense, 03 single stop-gain,
  04 compound-het LoF, **05 frameshift+rescue** *(hero case)*,
  06 compound missense (pocket), 07 splice donor loss + cryptic gain,
  08 symbolic `<DEL>` whole-gene, 09 upstream DEL coord-shift,
  10 `<INV>` gene body, 11 `<DUP:TANDEM>` small gene, 12 in-frame 3 bp del.
  Each case emits `{reference.fa, annotation.gff3, phased.vcf.gz, expected.json}`
  and is pinned by a dedicated pytest assertion.
  `benchmarks/paper_cases/run_summary.py` runs all 12 through the pipeline
  and emits a per-case TSV that becomes paper Table 1 — showing per-haplotype
  presence (`intact` / `deleted` / `duplicated` / `low_identity` / `uncertain`),
  diploid score, compound-het-LoF flag, per-variant consequence categories,
  and `frameshift_restored` where applicable.
- **MAVE benchmark** (`benchmarks/mave/`):
  - **Single-mutant**: BRCA1 RING/BRCT depletion (Starita 2015,
    `urn:mavedb:00000081-a-1`, 1055 missense variants).
  - **Double-mutant**: GB1 IgG-FC binding pairwise scan (Olson 2014,
    `urn:mavedb:00000105-a-1`, 315 759 double-mutant haplotypes).
- **Head-to-head with bcftools/csq** (`benchmarks/gnomad_flips/`): curated
  paper cases run through both tools, agreement matrix.
- **Synthetic SV-haplotype shape coverage** (`benchmarks/sv/`): 12 cases
  covering single missense, double missense, nonsense, frameshift,
  in-frame del, range del, large N-terminal del (whole-gene-loss proxy),
  compound LoF+missense.

### 2.3 Models + tools

- ESM2 (Lin 2023) — `esm2_t6_8M_UR50D` for development; `esm2_t33_650M_UR50D`
  for publication numbers (run on NVIDIA RTX 4090).
- AlphaMissense (Cheng 2023) — precomputed hg38 lookup table.
- bcftools 1.21 (consensus, csq), Liftoff 1.6.3 (gene lift-over, uses minimap2),
  pysam, Bio.Align.

## 3. Results

### 3.1 Head-to-head vs bcftools/csq on motivating cases (primary)

All 12 hand-designed cases in `data/paper_cases/` run through both bcftools/
csq (at the per-variant layer) and hapli (at both aggregation axes); reproduce
with `uv run python3 benchmarks/gnomad_flips/run_head_to_head.py --with-esm`.

| case                         | csq.hap1                  | csq.hap2    | hapli h1/h2   | hapli-only signals                                      |
|---                           |---                        |---          |---            |---                                                      |
| 01 synonymous_pair           | synonymous                | synonymous  | 1.00 / 1.00   | —                                                       |
| 02 single_benign_missense    | missense                  | —           | 0.99 / 1.00   | —                                                       |
| 03 single_stop_gain          | stop_gained               | —           | 0.20 / 1.00   | —                                                       |
| 04 compound_het_lof          | stop_gained               | frameshift  | 0.00 / 0.49   | **compound_het_lof = True**  *(Axis 2)*                 |
| 05 frameshift_rescue         | inframe_altering,missense | —           | 0.85 / 1.00   | **frameshift_region.restored = True**  *(Axis 1)*       |
| 06 compound_missense_pocket  | missense                  | —           | 0.96 / 1.00   | ESM residual +12 *(Axis 1 supplement, see §3.3)*        |
| 07 splice_donor_cryptic      | —                         | —           | 1.00 / 1.00   | —                                                       |
| 08 symbolic_del_whole_gene   | —                         | —           | 0.00 / 1.00   | **presence = deleted**  *(Axis 1)*                      |
| 09 upstream_del_shift        | —                         | —           | 1.00 / 1.00   | —                                                       |
| 10 inversion_body            | —                         | —           | 0.89 / 1.00   | **presence = low_identity**  *(Axis 1)*                 |
| 11 tandem_dup_small_gene     | —                         | —           | 1.00 / 1.00   | **presence = duplicated**  *(Axis 1)*                   |
| 12 inframe_domain_indel      | inframe_deletion          | —           | 0.99 / 1.00   | —                                                       |

Summary (`benchmarks/gnomad_flips/results/head_to_head.json`):
`hapli_compound_het_lof_calls=1, hapli_epistasis_flagged_calls=2,
csq_lof_haplotype_records=3`.

Five disagreement axes, each corresponding to an aggregation step csq
cannot perform:

- **Axis 2 aggregation (haplotype → genotype)**: case 04
  (compound_het_lof).
- **Axis 1 aggregation (variant → haplotype gene function)**: case 05
  (frameshift+rescue flagged via protein-diff), case 08 (gene-level DEL
  via Liftoff presence), case 10 (INV → low_identity), case 11 (DUP →
  duplicated). All four require scoring the *whole mutated gene*, which
  no per-variant record encodes.

Case 06 (compound-missense-pocket) flags the ESM residual as an optional
supplement; the generic version of that signal is evaluated as a null
result on bulk DMS in §3.3 below.

### 3.2 Per-variant scoring sanity check (BRCA1 RING)

ESM2 8M per-variant log-odds vs. Starita 2015 depletion labels:
**AUROC = 0.69** at LoF threshold "depleted in ≥ 2/4 replicates".

This is a baseline sanity check that ESM-1v / wildtype-marginals scoring
is behaving as published (AlphaMissense is in the 0.7 ballpark on the
same task). It does **not** evaluate either aggregation axis — it is the
per-variant layer, and the number is included only so downstream aggregation
results are not read as compensating for a broken per-variant scorer.

### 3.3 Optional ESM2 full-sequence score: bulk DMS is a null

The ESM residual — full-sequence pseudo-log-likelihood delta minus the
sum of per-variant log-odds — is offered as an **optional supplement**
on Axis 1 when a haplotype carries ≥ 2 coding substitutions of equal
length. We evaluate it on the GB1 pairwise DMS scan (Olson 2014,
MaveDB `urn:mavedb:00000105-a-1`, **n = 464,408 double mutants**,
ESM2 650M on RTX 4090, ~4.5 h wall across job 2948 on
`unimatrix01:node005`; reproduce via `benchmarks/mave/postprocess_650m.sh`).

| Predictor | AUROC vs `mave_score < 0.1` | Spearman vs MAVE | Spearman vs E_exp |
|---|---|---|---|
| sum of per-variant log-odds | 0.659 | −0.30 | −0.019 |
| full-sequence PLL delta | 0.648 | −0.27 | −0.008 |
| residual (PLL − sum) | 0.524 | −0.021 | +0.001 |

E_exp is the Olson-style experimentally-derived epistasis target,
`s(v1,v2) − ŝ(v1) − ŝ(v2) + s̄` with ŝ(v) estimated by marginal
averaging over all observed pairs containing v.

**Honest null.** On a 56-aa protein with exhaustive pairwise DMS coverage:
- The full-sequence PLL delta is not measurably better than the sum of
  per-variant log-odds (ΔAUROC ≈ −0.01).
- The residual is essentially uncorrelated with experimental epistasis
  (|ρ| ≈ 0.001).

**Caveats that matter.** The "sum of per-variant log-odds" comparator
here is not a clinical baseline — no Mendelian / tumor tool sums per-
variant pathogenicity scores within a gene. It is instead the model
output that mechanically isolates the "joint vs additive" axis for
scrutiny. On that axis, on this dataset, ESM2 650M has no detectable
advantage. On the short (56-aa) GB1 protein the two scores are
essentially redundant; longer proteins and structure-aware scorers
(SaProt) are the natural follow-ons. On *designed* epistatic haplotypes
the residual does fire reliably (§3.4).

*Historical note.* An earlier preliminary run on the first 2000 rows of
the GB1 CSV with ESM2 8M reported residual AUROC 0.76 and residual |ρ|
vs E_exp = 0.44 — impressive numbers that turned out to be a sampling
artifact of the head-of-file slice (`S_additive` |ρ| vs E_exp dropped
from 0.21 → 0.019 when extended to the full set, a 10× fall). Retained
as a transparency note on the fragility of small-sample MAVE validation.

### 3.4 Controlled paper-case residuals (where the residual does fire)

Run every paper-case bundle with `--with-esm`; record the residual per
haplotype. Expected signals designed into each case are recovered:

| case                         | n_vars | sum_per_variant | full_sequence_PLL | residual  | flagged |
|---                           |---     |---              |---                |---        |---      |
| 05 frameshift_rescue         | 6      | -58.48          | -1.43             | **+57.05** | ✓       |
| 06 compound_missense_pocket  | 2      | -21.91          | -9.70             | **+12.20** | ✓       |
| 07 splice_donor_cryptic      | 0      | 0               | 0                 | 0         | —       |
| 10 inversion_body            | 0      | 0               | 0                 | 0         | —       |

Reproduce: `uv run python3 benchmarks/paper_cases/run_epistasis.py`
(ESM2 8M on CPU, < 10 s / case). Cases 07 and 10 correctly emit n=0 —
the residual is defined for ≥ 2 coding substitutions of equal length,
which splice-only and SV cases by construction do not satisfy.

These designed cases are the illustrative limit: when two variants
structurally buffer each other (case 05 restored frame) or structurally
synergise (case 06 pocket clash), the residual is large in the expected
direction. On bulk averaged pairwise DMS over all 464K doubles, that
signal is apparently buried in per-variant noise (§3.3).

### 3.5 Variant-shape coverage (SV haplotype benchmark)

12 synthetic cases × 45 assertions, all PASS:
- Single missense, double missense, nonsense (early + late stop-gain),
  frameshift (with + without Ter offset), single-residue del, range del,
  large N-terminal del (SV whole-gene-loss proxy), compound LoF+missense,
  out-of-range/skipped, residual-is-None for length-changing haps.

### 3.6 Computational performance

- 8M model on CPU: ~32 ms / haplotype on a 221-aa protein.
- 650M model on GPU (RTX 4090): *32 ms / haplotype* (same wall time;
  the larger model is GPU-batched while the small model is CPU-bound).
- Per whole-genome sample, end-to-end: **2–4 h** at 16 cores + 1 GPU.

## 4. Discussion

### 4.1 What hapli adds over prior art

Hapli's contribution is at the **aggregation layer**, not at the
per-variant scoring layer. Two specific axes:

**Axis 1 — variants within a haplotype → whole-gene function.**

Prior art at this layer:
- LOFTEE worst-consequence: picks the worst HGVS cascade element.
- AlphaMissense MAX: take the most-pathogenic per-missense score.
- bcftools/csq: concatenates phased consequences; does not collapse.

Hapli replaces the hand-rolled collapse with a *protein-first* pass:
materialise the haplotype FASTA (bcftools consensus + symbolic INV/DUP
resolution), lift the GFF3 onto each haplotype (Liftoff), extract and
translate each transcript's CDS, then compute scalars of the whole
mutated protein:

- `identity` vs reference protein (pairwise alignment).
- `premature_stop_at` with N-terminal-LoF scaling.
- `frameshift_region.restored` — detects +1/−1 indel pairs that
  restore frame (case 05). No per-variant aggregation rule reaches this,
  because both variants individually still cascade to frameshift LoF.
- `presence` (from Liftoff) — `intact / low_identity / partial /
  deleted / duplicated`. Surfaces SV-level signal (cases 08, 10, 11)
  that is silent in any per-variant output.
- Optional ESM2 full-sequence score for compound missense; fires
  reliably on designed epistatic cases (case 06, residual +12), null
  on bulk pairwise DMS of a short protein (§3.3).

**Axis 2 — two haplotypes → one genotype.**

Prior art: none at the tool layer. Dominance rules are implemented as
downstream scripts on top of per-variant records; compound-het-LoF
detection in particular requires joining two haplotype-of-origin records,
which every clinical pipeline writes itself.

Hapli emits `{hap1_score, hap2_score, min, max, compound_het_lof,
constraints}` as a first-class JSON field. Deliberately, no single
dominance call is committed to — clinical consumers still apply their
own rule. The headline axis-2 signal is the `compound_het_lof` flag
(case 04).

**Framing note.** The §3.3 ESM2-joint-vs-sum benchmark is not a clinical
comparator (no Mendelian / tumor tool sums per-variant pathogenicity
within a gene — summation is specifically a PRS / burden-test idiom).
We keep it in the paper only to report an honest null on a mechanistic
decomposition of the per-haplotype score; the two aggregation axes above
are where the paper actually lives.

### 4.2 Limitations

- Mode A's `bcftools consensus` cannot materialise BND records (translocations);
  Mode B (pre-assembled haplotypes) is the answer for those. `<CNV>`
  (allele-count rather than allele-state) is similarly out of scope in Mode A;
  symbolic `<INS>` without an explicit ALT sequence is rejected at the
  preflight stage and must be pre-resolved upstream (pangenome graph /
  long-read-supported caller).
- Hapli's `_resolve_symbolic_svs` preprocessor rewrites `<INV>` and
  `<DUP:TANDEM>` into explicit REF/ALT sequences so bcftools consensus
  (≤ 1.21, which refuses these symbolic ALTs) can apply them. This is a
  linear encoding of each event in isolation; it does not capture interactions
  between overlapping SVs, nor alternative resolutions of ambiguous
  breakpoint pairs. True graph-aware SV reconciliation (pangenome graphs via
  vg / Minigraph-Cactus) is the correct long-term answer and is deferred to v2.
- ESM2 pseudo-log-likelihood is well-behaved on substitutions but undefined
  on length-changing haplotypes; protein-diff signals (premature_stop_at,
  identity, frameshift_region.restored) carry information for those cases
  instead. The residual is, by construction, computed only when ≥ 2 coding
  substitutions share a transcript on the same haplotype — splice-altering
  (case 07) and SV-bisecting (case 10) cases correctly emit n=0 residual
  but may still be epistatic in a biological sense not captured here.
- The "joint beats additive" claim inherits any bias of the protein-language
  model used as the judge. ESM2 is trained on evolutionary co-occurrence;
  co-evolved residue pairs get structural-compensation bonuses even if the
  function is disrupted. The residual therefore over-predicts rescue on
  co-evolved pairs and under-predicts synergy on novel interactions. A SaProt
  (structure-aware) comparison is a planned ablation.
- The diploid score's LoF threshold (0.5) is a heuristic chosen to reproduce
  a "below-this-is-LoF" aggregate; clinical application requires per-gene
  calibration against MAVE or ClinVar-curated variants, ideally gene-by-gene.
- Liftoff is the authoritative source of per-haplotype gene-annotation lift;
  a mis-lift (e.g., gene incorrectly mapped across an inversion boundary)
  propagates to protein extraction and the residual. The planned miniprot
  protein-level cross-check (stub only in v1) is the intended safeguard.
- E_exp marginal-averaging assumes the GB1 pairwise scan covers each single
  position uniformly (the synonymous-self-pair pad is a small approximation).
- The paper's core epistasis claim is currently shown on one MAVE dataset (GB1).
  PSD95 PDZ3 (the natural second pairwise scan) was not fetchable from the
  MaveDB API for this submission and is documented as a follow-on.

### 4.3 Future directions

- **HPRC scale demo** [pending]: Snakefile is ready (`mode: B` config); needs
  a 200-sample assembly fetch + ~24h on the cluster.
- **gnomAD compound-het flip reproduction** [pending]: full Danecek 2017
  protocol on 1000G Phase 3 (~80 GB download) — Danecek reported that
  phasing-aware consequence calling flipped **501 / 5019** compound-variant
  predictions on 1000G; the hapli replication would add the aggregated
  diploid `compound_het_lof` tag and the epistasis residual as additional
  axes of flipping on top of csq's baseline. Runbook in
  `benchmarks/gnomad_flips/README.md`.
- **Saudi Human Genome Program**: population-specific compound-variant
  flip catalog; access-restricted, planned as discussion-section vignette.
- **SpliceAI / LOFTEE / TOGA wrappers**: stubs in `hapli/external/`;
  full integration is a self-contained sprint.

## 5. Reproducibility

Everything in `hapli` is one `uv sync && uv tool install git+https://github.com/agshumate/Liftoff` away from runnable. The benchmarks are
each a `./run.sh <urn>` or `bash benchmarks/.../run.sh`. Tests:
`uv run pytest` (171 passed at submission). Coverage highlights: 12 end-to-end
paper-case tests, resolver pins for `<INV>` / `<DUP:TANDEM>`, presence-
identity-capping regressions for low-identity/partial Liftoff lifts, Mode A
and Mode B Snakefile end-to-end smoke tests (`analyze → aggregate` and
`assess → aggregate`), LLM-prompt regressions (compound-het-LoF flag,
epistasis residual, ClinVar clnsig), TUI schema-v2 construction check,
direct unit tests for the minimap2 aligner wrapper and GFFProcessor.
`testpaths = ["tests"]` in `pyproject.toml` so `uv run pytest` without args
DTRT (excludes obsolete/). Walkthrough:
`uvx showboat verify docs/walkthrough.md` (executable end-to-end demo
on the data/test/ fixture, regression-tested).

## Headline numbers for the abstract / cover page

**Primary results (the paper's actual claims):**

| | value |
|---|---|
| **Axis-2 signal — compound-het-LoF calls hapli emits that csq cannot** | case 04 (1/1 motivating case) |
| **Axis-1 signal — frameshift-rescue calls hapli emits that csq cannot** | case 05 (`frameshift_region.restored=True`) |
| **Axis-1 signal — SV presence calls hapli emits that csq cannot** | cases 08 (`deleted`), 10 (`low_identity`), 11 (`duplicated`) |
| **Total disagreement axes with bcftools/csq on 12 motivating cases** | 5 (1× Axis-2, 4× Axis-1) |
| **Per-variant ESM2 scoring sanity check (BRCA1 RING 8M)** | AUROC 0.69 |
| **Compute per whole-genome sample** | 2–4 h (16-core + 1 GPU) |
| **Tests passing** | 172 |
| **SV-haplotype shape coverage** | 12 cases, 45/45 assertions |

**Honest-null and illustrative numbers (§3.3, §3.4):**

| | value |
|---|---|
| Full GB1 pairwise — sum-of-per-variant AUROC (n=464,408) | 0.659 |
| Full GB1 pairwise — full-sequence PLL AUROC (n=464,408) | 0.648 (flat) |
| Full GB1 pairwise — residual \|ρ\| vs E_exp | 0.001 (null) |
| Designed case 05 (frameshift rescue) residual | +57.05, flagged |
| Designed case 06 (compound missense pocket) residual | +12.20, flagged |
| GB1 doubles scored at submission | 464,408 (ESM2 650M on RTX 4090, ~4.5 h wall) |
