# hapli — Nature Communications "Structural genetic variants" submission plan

> Target collection: <https://www.nature.com/collections/djjhfigfge>
> Deadline: **31 August 2026** (17 weeks from 2026-05-01)
> Repo state: paper outline at `docs/paper_outline.md`; HPRC R2 (20×80) and 1000G ACMG (100×80) demos already complete. Pipeline architecture: `hapli/workflow/pipeline.py`.

---

## 0. Framing

The current paper outline (`docs/paper_outline.md`) is built around **two aggregation axes** (variants → haplotype, haplotypes → genotype). For the *Structural genetic variants* collection, the value proposition we sell on is **Axis 1's SV-presence signal**: hapli emits per-haplotype `presence ∈ {intact, partial, low_identity, deleted, duplicated}` calls that no per-variant tool produces. SVs are exactly the variant class where the per-variant aggregation rules underlying LOFTEE/AlphaMissense fail (csq is silent on symbolic records).

We will **re-foreground** the paper around SV-driven gene-function loss without removing the genotype-level (Axis 2) story. Concretely:

- Title pivots from "two aggregation axes" → "**hapli: SV-aware per-haplotype gene-function calls from short-read VCFs and long-read assemblies**".
- The frameshift+rescue + compound-het-LoF cases are kept as Axis-1 / Axis-2 demonstrators but moved into a unified "what hapli sees that csq cannot" subsection.
- ESM2 epistasis + AlphaMissense aggregator are *kept as supplementary material only* — they are the next paper.

**Hard constraint**: 17 weeks must include experiments + writing + internal review + revision pad. Budget: ~10 weeks experiments, 4 weeks writing, 2 weeks review/revision pad, 1 week submission logistics.

---

## 1. Tiered experiment list (locked)

### Tier 1 (must-do; cuts trigger paper non-submission)

| # | Experiment | Status | Owner | Wall-time est. |
|---|---|---|---|---|
| T1.1 | **HPRC R2 SV-focused vignette** — re-slice existing HPRC R2 outputs against the HPRC R2 dipcall SV truth set, intersect with `presence ∈ {deleted, low_identity, duplicated}` calls, build per-SV-class confusion matrix and a curated "gallery of 10 hapli-flagged gene losses" figure | New analysis on **existing** Mode B outputs + one new tool wrapper | — | 2 cluster days + 1 wk human curation |
| T1.2 | **Cancer long-read cohort** (user data, in hand) — Mode B (if assemblies available) or Mode A on tumor variant calls; emit per-gene tumor-vs-normal differential `presence` + `score`; tumor-suppressor-panel focus (TP53, RB1, PTEN, BRCA1/2, ATM, CDKN2A, NF1, APC, VHL, STK11) | New mode + new wrappers + new cohort | — | 5 cluster days |
| T1.3 | **AnnotSV + SvAnna head-to-head** — same VCF in, compare per-gene SV impact calls. Show the genotype-aware aggregation hapli adds on top | New wrappers + new benchmark | — | 3 days |
| T1.4 | **Synthetic SV-shape benchmark, expanded** — already 12 cases / 45 assertions; add 8 more covering tandem-segmental dup, exonic dup-tandem, intronic SV, partial-CDS deletion, in-frame deletion of >1 exon, gene-fusion-via-deletion, paired-allele compound SV+SNV | Existing scaffold (`benchmarks/sv/`) | — | 2 days |
| T1.5 | **HPRC R2 scaling**: bump 20 → 60 samples, full ACMG SF v3.2 panel, retain LoF-AF-vs-gnomAD calibration story | Existing pipeline; cluster only | — | ~3 cluster days (wall) |

### Tier 2 (strong addition — go/no-go at week 4)

| # | Experiment | Trigger to drop |
|---|---|---|
| T2.1 | **Rare-disease short-read cohort** (user data, in hand) — proband-only Mode A; positive-control: does `compound_het_lof=True` recover known recessive diagnoses where available? Run as discovery analysis if no truth labels | Drop if cohort lacks any genes-of-interest annotation by week 6 |
| T2.2 | **gnomAD v4 LoF-AF calibration figure** — across both 1000G and HPRC, plot hapli per-gene LoF-AF vs gnomAD pLoF-AF, separately for SNVs/indels and SV-driven losses | Drop if HPRC scaling (T1.5) doesn't finish |

### Tier 3 (nice-to-have — only if Tier 1+2 are wrapped by week 11)

- T3.1: Rare-disease **long-read** acquisition + analysis. Cut from this paper unless sequencing is already scheduled and turn-around < 6 weeks. Default position: **drop**.
- T3.2: 1000G Phase 3 expansion 100 → 500 samples. Pure scaling; only if a marginal cluster window is free.

### Cut from this submission (do not negotiate)

- ESM2 epistasis residual on cohorts (next paper).
- AlphaMissense aggregator on cohorts (next paper).
- Pangenome-graph workflows.
- BND / translocation reconciliation.
- UDN dbGaP-restricted runs.

---

## 2. Code / infrastructure milestones (independent of experiments)

These ship to `main` regardless of which Tier 2/3 experiments survive. All paths follow existing conventions in `hapli/external/` (one wrapper per external tool, `*NotAvailable` exception, dataclass result, optional logger).

### 2.1 `hapli/external/annotsv.py` — AnnotSV wrapper

```python
class AnnotSVNotAvailable(RuntimeError): ...

@dataclass
class AnnotSVResult:
    annotated_tsv: Path
    per_sv_calls: list[AnnotSVCall]   # SV_chrom, SV_start, SV_end, SV_type, gene, ACMG_class

def run_annotsv(
    vcf_path: Path, *, genome_build: str = "GRCh38",
    annotation_dir: Path | None = None,
    out_tsv: Path | None = None,
    annotsv_path: str = "AnnotSV",
    logger: logging.Logger | None = None,
) -> AnnotSVResult: ...
```

Mirrors `hapli/external/csq.py::run_csq` shape exactly. AnnotSV install via conda (`bioconda::annotsv`); a `*NotAvailable` raise lets pipeline degrade gracefully like `LiftoffNotAvailable`. Test fixture: `tests/test_annotsv.py` with a tiny SV VCF + recorded AnnotSV TSV stub (don't shell out in CI).

### 2.2 `hapli/external/svanna.py` — SvAnna wrapper

Same shape. Java tool — wrap via `subprocess` exactly like Liftoff. SvAnna emits per-SV pathogenicity + per-gene impact; we want the per-gene impact join. Less mature ecosystem than AnnotSV; keep wrapper but clearly mark as a *secondary* comparator if integration is brittle.

### 2.3 `hapli/external/dipcall_sv.py` — HPRC R2 truth-set parser

Not a tool wrapper but a parser for the HPRC R2 dipcall SV BED/VCF outputs. Parses the per-sample SV truth set so T1.1's confusion matrix has a denominator. Tiny module; conventions match `hapli/external/clinvar.py::ClinVarLookup`.

### 2.4 Cancer-mode flag — pipeline + CLI

Modify `hapli/workflow/pipeline.py::HapliPipeline.__init__` to accept a `tumor_sample: str | None` and `normal_sample: str | None`. When both set: run the per-(sample, gene) pipeline twice on the same VCF (once with `--samples tumor`, once with `--samples normal`), emit a third differential record under `evidence.tumor_normal`:

```python
class TumorNormalDiff(BaseModel):
    tumor_presence: dict[str, PresenceCall]      # hap1, hap2 — tumor calls
    normal_presence: dict[str, PresenceCall]
    delta_score: float | None                    # tumor.min_score - normal.min_score
    new_lof_alleles: int                         # LoF in tumor not in normal
    loh_event: bool                              # one allele lost in tumor only
```

Goes into `hapli/core/schema.py::GeneEvidence` as `tumor_normal: TumorNormalDiff | None`. New CLI subcommand:

```bash
hapli analyze-cancer --gene TP53 --tumor T01 --normal N01 \
    --vcf somatic.vcf.gz --reference ref.fa --gff anno.gff3 \
    --output-dir results/
```

Exposed in `hapli/cli/main.py` alongside `analyze`/`assess`. Internal implementation: run `analyze` twice with different `--sample` and merge the two `AnalysisResult` JSONs into one. Reuses the existing whole-haplotype consensus + Liftoff cache via `_file_lock` already in `consensus.py`.

**Phased / unphased VCF handling for somatic data**: somatic VCFs are typically unphased. The current `HaplotypeBuilder` treats unphased GT as if it were `hap1/hap2` ordered (commented as "unsafe"). For somatic VAFs <0.5 we cannot ascribe variants to specific haplotypes anyway, so the pragmatic policy is:

- Treat tumor sample as **diploid-bag**: both ALT alleles applied to one synthetic "tumor consensus" haplotype, normal applied to the other.
- Emit a `phasing: "tumor_unphased"` flag on the `TumorNormalDiff`.
- Score with min(hap_tumor, hap_normal) — preserves the LoF detection without overclaiming compound-het.

Document this clearly in the methods. It is not a phasing solution; it is a defensible policy choice for the SV story (where the lesion is allele-state regardless of phase).

### 2.5 SV-evidence reporting in JSON — link presence calls to the SV records that caused them

Currently `evidence.presence[hap1]` says `status=deleted` but doesn't tell you *which* `<DEL>` record at which coordinates was responsible. For the SV paper this is load-bearing.

Add to `hapli/core/schema.py::PresenceCall`:

```python
class CausingSVRecord(BaseModel):
    chrom: str
    pos: int                          # 1-based
    end: int
    sv_type: Literal["DEL", "INV", "DUP", "INS", "BND", "CNV", "complex"]
    sv_id: str | None                 # VCF ID column
    length: int                       # END - POS
    overlap_fraction: float           # of the gene's CDS overlapped

class PresenceCall(BaseModel):
    ...                                # existing fields
    causing_svs: list[CausingSVRecord] = Field(default_factory=list)
```

Population in `hapli/workflow/pipeline.py::_run_liftoff`: after computing `PresenceCall`, for any non-`intact` status, scan the input VCF over the gene region and attach SV records that overlap the CDS on the same haplotype. New helper in `hapli/external/consensus.py::find_overlapping_svs(vcf_path, sample, hap_idx, region) -> list[CausingSVRecord]`. Mode B: walk Liftoff's intermediate PAF instead and derive the gap regions.

This single change is what makes the SV story machine-readable in the paper TSV and is what SvAnna/AnnotSV will be compared against gene-by-gene.

### 2.6 Snakefile changes for cancer cohort

Add a `mode: "cancer"` branch to `workflows/Snakefile` mirroring the existing Mode A/B split:

- `config["pairs"]` is a `{tumor_id: normal_id}` mapping.
- `analyze_cancer_gene` rule shells out to `hapli analyze-cancer ...`.
- Aggregator `hapli aggregate` already groks any `*_<gene>_analysis.json`; extend it to recognise `evidence.tumor_normal` and emit a `per_pair_per_gene.tsv` with LoH and new-LoF columns.

Touch points:
- `workflows/Snakefile` — new mode branch + new rule.
- `workflows/config.yaml` — example cancer block.
- `hapli/workflow/aggregate.py::SampleGeneRow` — add `tumor_loh: bool`, `tumor_new_lof: int`.
- `hapli/cli/main.py::aggregate` — already glob-driven; no change.

### 2.7 Testing

- `tests/test_annotsv.py` — parser unit tests with frozen TSV fixtures.
- `tests/test_svanna.py` — same.
- `tests/test_cancer_mode.py` — synthetic tumor + normal VCF, asserts LoH detection on a known-engineered single-allele DEL on tumor.
- `tests/test_sv_evidence.py` — paper case 08 (`<DEL>` whole-gene) emits exactly one `CausingSVRecord` with the right coordinates.
- Integration: extend `tests/test_snakefile_e2e.py` with a cancer-mode end-to-end on the synthetic fixture.

Target: 195+ tests passing by week 6 (currently 179).

---

## 3. Datasets — concrete URLs + versions

### Already in repo / on cluster

- **HPRC R2 phased assemblies**: `s3://human-pangenomics/working/HPRC_R2/` — sample manifest in `benchmarks/hprc/sample_list.tsv` (20 samples, scalable to 200+).
- **1000G NYGC phased VCFs**: `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/` — already region-subsetted in `data/1000g_acmg/`.
- **GRCh38 + GENCODE v45 GFF3** — already on disk per `benchmarks/1000g_acmg/README.md`.
- **gnomAD v4 constraint** — `https://gnomad.broadinstitute.org/downloads` → `gnomad.v4.0.constraint_metrics.tsv`.
- **ClinGen dosage** — `ftp.clinicalgenome.org/.../ClinGen_gene_curation_list_GRCh38.tsv`.

### Required for SV truth (T1.1)

- **HPRC R2 dipcall SV calls** per sample: `s3://human-pangenomics/working/HPRC_R2/<sample>/assemblies/year1_freeze_assembly_v2/dipcall/` (per-sample VCFs of the assembly-vs-GRCh38 SV catalog). For each of the 20 (→60) samples: `<sample>.dip.vcf.gz` + `<sample>.dip.bed`.

### Required for AnnotSV/SvAnna (T1.3)

- **AnnotSV annotation bundle** v3.4: `https://www.lbgi.fr/AnnotSV/Documentation/AnnotSV_3.4.tar.gz` (~3 GB). Install via conda `bioconda::annotsv` and run once to populate `$ANNOTSV/share/AnnotSV/Annotations_Human/`.
- **SvAnna v1.0.4** jar: `https://github.com/TheJacksonLaboratory/SvAnn/releases`. JDK 17 on cluster.
- **SvAnna data files**: download via `svanna download --build hg38`.

### TBD — user must confirm by **end of week 1** (hard gate)

- **Cancer LR data** (T1.2):
  - Sequencing platform (PacBio HiFi vs ONT vs both)?
  - Sample count?
  - Tumor/normal pairing? Matched normal blood/buccal?
  - Cancer type(s)?
  - Pre-existing variant calls (from where? PEPPER-DeepVariant? Sniffles2? PBSV?)
  - Pre-existing assemblies (hifiasm? Verkko?) → if yes, **prefer Mode B**
  - Ethics / data-use status — IRB clearance for redistribution of derived call statistics? Cluster storage path (IBEX scratch?)
- **Rare-disease SR data** (T2.1):
  - Sample count, family structure (probands only? trios?)
  - Phenotype labels available (HPO terms?)
  - Known molecular diagnoses (positive controls)?
  - Phasing status (statistically phased or unphased?)

### Will not acquire for this submission

- Long-read rare-disease (T3.1) — cut unless user reports week-1 that sequencing is already scheduled with turnaround < 6 weeks.

---

## 4. Compute budget

IBEX 16-core × 128 GB nodes, `pi-hohndor` account, `batch` partition. Reference numbers from completed runs:

- Mode A (1000G): **21 h wall** for 100 samples × 80 ACMG genes, 16 cores, one batch node.
- Mode B (HPRC): **~17 min wall** for 20 samples × 80 genes once Liftoff is cached, plus **~30 min Liftoff array**.

| Experiment | Samples | Genes | Mode | Est. wall | Cluster days |
|---|---|---|---|---|---|
| T1.1 (re-analysis only) | 20 | 80 | n/a (post-hoc) | <2h on a single node | 0 |
| T1.5 (HPRC scale-up to 60) | 60 | 80 | B | ~3h Liftoff + ~1h analysis | 1 |
| T1.2 (cancer LR, est. 20 pairs) | 20 pairs (40 samples) | full ACMG + tumor-suppressor panel ~120 genes | A or B (TBD) | scaling Mode A: 2 × 21h × 0.4 ≈ 17h on 4 nodes; Mode B: 2h | 2-5 (depending on mode) |
| T1.3 (AnnotSV/SvAnna head-to-head) | 100 (1000G existing) | ACMG | n/a | minutes | 0 |
| T1.4 (synthetic SV expansion) | n/a | n/a | n/a | minutes | 0 |
| T2.1 (rare-disease SR, est. 50 probands) | 50 | full ACMG + RD panel ~300 genes | A | ~10h on 4 nodes (parallel) | 2 |

**Total cluster footprint**: ~10 cluster-days across the 4 months. Comfortably fits one rolling allocation. Storage: budget ~2 TB scratch (HPRC R2 60 samples × 30 GB for assemblies + cancer cohort).

---

## 5. Week-by-week schedule

> Today (assumed start of execution): **Mon 2026-05-04 (week 1)**.
> Deadline: **Mon 2026-08-31 (week 17)**.

### Week 1 (May 4–10): Lock requirements + scaffolding

- **Hard decision**: user confirms cancer LR cohort metadata (platform, count, pairing, ethics, paths). **No further experimental work proceeds without this.** If cohort is not Mode-B-ready and Mode A turnaround would exceed week 8, decide between (a) Mode A on small cohort (b) cohort + assembly via hifiasm (c) cut T1.2.
- Submit AnnotSV + SvAnna installation MR (`hapli/external/annotsv.py`, `hapli/external/svanna.py`) — wrappers + tests, no live tool integration yet.
- Install AnnotSV + SvAnna binaries on IBEX, write smoke tests.
- Branch + scaffold `hapli/core/schema.py::CausingSVRecord` + `find_overlapping_svs` helper in `consensus.py`.
- **Writing**: Begin updating `docs/paper_outline.md` toward the SV-foreground framing. Title + abstract first draft.

### Week 2 (May 11–17): SV evidence wiring + HPRC scaling kickoff

- Land §2.5 (causing-SV evidence) on main; backfill into the 12 paper cases — case 08, 10, 11 should now emit explicit `CausingSVRecord`.
- Pytest regressions for SV-evidence schema additions.
- Kick off HPRC R2 sample-list expansion (20 → 60). Run `benchmarks/hprc/fetch.py --refresh-manifest`.
- Begin downloading 40 new HPRC R2 assemblies (~50 GB).
- T1.4: synthetic SV-shape benchmark expansion — design 8 new cases in `scripts/generate_paper_cases.py` extending the registry.

### Week 3 (May 18–24): T1.1 vignette + HPRC scale run

- Write the HPRC dipcall SV truth-set parser (`hapli/external/dipcall_sv.py`).
- Run HPRC scale-up on cluster (T1.5, ~1 cluster day, ~3 days end-to-end with download + Liftoff cache build).
- T1.1 starts as soon as HPRC outputs land: confusion matrix `presence ∈ {deleted, low_identity, duplicated}` × `dipcall SV ∈ {DEL, INV, DUP, INS}` per gene; gallery curation begins.
- AnnotSV smoke test on 1000G subset.

### Week 4 (May 25–31): Cancer mode + T2 go/no-go

- Land §2.4 (cancer mode + `analyze-cancer` CLI) + §2.6 (Snakefile branch) on main.
- Tests `tests/test_cancer_mode.py` (synthetic LoH detection).
- Smoke-run cancer pipeline on 1 tumor/normal pair from user data.
- **Hard decision (week 4 close)**: T2.1 (rare-disease SR) go/no-go. Criteria:
  - Cohort metadata + access confirmed by user.
  - At least 1 known molecular diagnosis available as positive control, or willingness to publish as a discovery analysis without ground truth.
  - If criteria not met: cut T2.1. T2.2 (gnomAD calibration) becomes the only Tier 2 deliverable.

### Week 5 (June 1–7): T1.1 figure lock + cancer cohort kick-off

- T1.1 confusion matrix + curated 10-gene gallery → paper Figure 4.
- T1.1 README + reproducibility script lands at `benchmarks/hprc/sv_vignette.py`.
- Cancer cohort: full T1.2 run on cluster (5 cluster-days budget).
- AnnotSV head-to-head (T1.3): start design — pick comparison gene × sample matrix.

### Week 6 (June 8–14): T1.3 head-to-head + T1.2 first cut

- T1.3: run AnnotSV on the 1000G ACMG VCF subset; build agreement matrix (per-gene per-sample SV impact) — hapli vs AnnotSV vs SvAnna.
- T1.2 first cut: per-pair tumor-vs-normal `presence` differential TSV; sanity-check that known tumor-suppressor LoH events surface.
- **Writing milestone**: Methods section first draft (sections 2.1–2.4 in `docs/paper_outline.md`).
- Internal mini-review of methods draft with collaborators.

### Week 7 (June 15–21): T2 execution

- If T2.1 alive: rare-disease SR cohort run (~2 cluster days).
- If T2.2 alive: build the LoF-AF-vs-gnomAD figure across HPRC + 1000G.
- Post-hoc analysis layer for cancer cohort: which tumor-suppressor genes lost on which tumor types? Pull per-cancer-type aggregate.

### Week 8 (June 22–28): All experiments closed; tables freeze

- Final rerun of all benchmarks for paper-final numbers (Snakemake `--rerun-incomplete`).
- **Lock**: Tables 2-5 (per-gene scaling, head-to-head, cancer LoH summary, rare-disease findings).
- Begin assembling supplementary tables.

### Week 9 (June 29–July 5): Figures lock

- All paper figures rendered from final tables — matplotlib scripts living next to each `benchmarks/*/analyze.py`.
- Style pass (consistent palette, font sizes, panel labels).
- **Lock**: Figures 1-5 + supplementary figures.

### Week 10 (July 6–12): Full first draft

- Complete first draft of all sections: Abstract, Intro, Methods, Results, Discussion.
- Reproducibility checklist completed (every numerical claim traceable to a `benchmarks/*/analyze.py` invocation; commit hashes pinned).
- Code: tag `v1.0.0-rc1`.

### Week 11 (July 13–19): Internal review #1

- Distribute draft to internal reviewers (PI + 2-3 collaborators).
- Address feedback in real-time.

### Week 12 (July 20–26): Revision pass + supplementary lock

- Address week-11 feedback.
- Lock supplementary methods, supplementary figures, supplementary tables.
- Generate static reproducibility archive (Zenodo DOI request).

### Week 13 (July 27–Aug 2): Internal review #2 + final bug fix pad

- Final internal review — sign-off from PI + co-authors.
- Final bug-fix pad: known unknowns, edge-case rerun budget.
- Run full pytest suite + showboat walkthrough one last time on a clean checkout.

### Week 14 (Aug 3–9): External review (optional)

- If a friendly external reviewer is available (e.g., HPRC consortium contact), send for outside-eyes pass.

### Week 15 (Aug 10–16): Final revisions

- Incorporate external feedback if any.
- Final figure tweaks.
- Submission system formatting (Nature Comm. requires LaTeX or Word; figures at journal-spec resolution).

### Week 16 (Aug 17–23): Submission package assembly

- Cover letter draft.
- Author contributions matrix.
- Data availability statement (paths to all datasets, code Zenodo DOI, GitHub release).
- COI statements.
- Suggested reviewers list (5 names with rationale).

### Week 17 (Aug 24–31): Submission

- 24th: Final read-through.
- 25-27th: Last-minute fixes (revision pad).
- 28th: Submit.
- 29-31st: **Buffer** for system errors / requested format changes.

---

## 6. Risk register

| # | Risk | Probability | Impact | Mitigation |
|---|---|---|---|---|
| R1 | Cancer LR data ethics/access not cleared by week 1 | Medium | Critical — kills T1.2 | Hard week-1 gate; if not cleared, swap T1.2 with T2.1 (rare-disease SR) as the third "real-data" Tier-1 |
| R2 | HPRC R2 dipcall SV truth set is incomplete or inconsistent across the 60 samples | Low-medium | Medium — weakens T1.1 | Fall back to 20-sample subset where dipcall is consistent; report as known limitation |
| R3 | AnnotSV install / runtime issues on cluster | Medium | Low — could degrade T1.3 | Wrapper degrades to `*NotAvailable` like Liftoff; SvAnna covers comparator ground if AnnotSV blocks |
| R4 | Cancer somatic VCFs are unphased + low-VAF subclonal — `compound_het_lof` is meaningless | High | Medium | Document explicitly in methods (§2.4 above); report `tumor_loh` and `tumor_new_lof` instead, which don't require phasing |
| R5 | Liftoff misses gene boundaries in cancer-rearranged genomes (high SV burden); spurious `low_identity` | Medium | Medium | Cross-check with the AnnotSV BED on the same VCF; manually curate the top 10 disagreements |
| R6 | IBEX scheduled maintenance + queue contention eats > 1 cluster week | Medium | Medium | The 21h-on-1-node budget plus our 10-cluster-day plan has 4× slack. Keep Snakemake `--rerun-incomplete` discipline |
| R7 | One reviewer / co-author unreachable during week 11 review | Medium | Low | Schedule reviews with 1-week notice; have a backup external reviewer in week 14 slot |
| R8 | Schema change to add `CausingSVRecord` breaks existing TUI/LLM consumers | Low | Low | `extra="allow"` on every pydantic model already; field is additive. Tested via `tests/test_llm_prompt.py` |
| R9 | Whole-haplotype consensus disk usage (3.3 GB × 2 × 60 samples = ~400 GB scratch) | Medium | Low | Cluster scratch budget allows; clean up after Liftoff via `rm -f *.whole.fa` post-run hook |
| R10 | Reviewers demand a per-variant tool comparison (LOFTEE, SnpEff) we haven't run | Medium | High during revision | Pre-build LOFTEE wrapper stub now (week 5) so revision turnaround is tractable. Already noted in `docs/paper_outline.md` §4.3 as a planned wrapper |

---

## 7. Cut-down plan: minimum-viable submission

**Hierarchy of cuts** if month-3 execution exposes a Tier-1 blocker:

### Level 1 cut (can lose without dropping submission)

If T1.2 (cancer LR) cannot complete by week 8: substitute the cancer real-data demonstration with a **methods-only cancer vignette** — 1-2 pairs as a proof-of-concept rather than a cohort. The paper still claims SV+haplotype-level LoH detection with the full pipeline; it just doesn't claim cohort statistics. Discussion-section vignette.

### Level 2 cut (lose Tier 2 entirely)

If both T2.1 and T2.2 collapse: paper is HPRC vignette + cancer demo + AnnotSV head-to-head + synthetic. This is the **minimum viable submission**.

### Submission floor (below which we do not submit)

The paper must contain at least:

1. The 12 motivating cases + head-to-head with bcftools/csq (existing — `benchmarks/paper_cases/`).
2. HPRC R2 SV vignette (T1.1) with at least 20 samples of curated SV-driven LoF events with `CausingSVRecord` evidence chain — this is the SV story's anchor.
3. Synthetic SV benchmark with ≥ 12 cases / ≥ 45 assertions (existing).
4. **Either** AnnotSV/SvAnna head-to-head (T1.3) **or** a cancer-cohort demonstration (T1.2).

If neither (3) nor (4) are achievable by week 11, submission shifts to a different (more methods-only) venue (eg Bioinformatics, GenomeBio) rather than this collection. **Do not submit a Tier-1-incomplete paper to a structural-variants collection.**

---

## 8. Hard decision points (calendar)

| Date | Owner | Decision |
|---|---|---|
| End of week 1 (May 10) | User | Cancer LR cohort: type, count, ethics, paths |
| End of week 4 (May 31) | User + lead | T2.1 rare-disease SR: keep or drop |
| End of week 5 (June 7) | Lead | T1.1 vignette: dipcall truth set workable on full 60? If not, lock to 20 |
| End of week 6 (June 14) | Lead | T1.2 cancer scope: tumor-suppressor panel only, or full ACMG? |
| End of week 8 (June 28) | Lead | Tables freeze — no new experiments after this point |
| End of week 9 (July 5) | Lead | Figures freeze |
| End of week 10 (July 12) | All | First draft complete, distributed |
| End of week 13 (Aug 2) | PI | Submission go/no-go |

---

## 9. Critical files for implementation

These are the files a follow-on engineer needs to open first:

- `/home/leechuck/Public/software/hapli/hapli/workflow/pipeline.py` — the heart of the pipeline; cancer-mode wiring, `CausingSVRecord` population, and Mode B branch all touch here.
- `/home/leechuck/Public/software/hapli/hapli/external/consensus.py` — where the new `find_overlapping_svs` helper lives and where the cancer somatic-VCF policy is enforced.
- `/home/leechuck/Public/software/hapli/hapli/core/schema.py` — every new evidence class (`CausingSVRecord`, `TumorNormalDiff`) is added here; downstream consumers (TUI, LLM, aggregator) read defensively.
- `/home/leechuck/Public/software/hapli/hapli/workflow/aggregate.py` — population aggregator; needs new `tumor_loh` / `tumor_new_lof` columns + extended per-gene SV-impact tally.
- `/home/leechuck/Public/software/hapli/workflows/Snakefile` — Mode A/B/cancer fan-out; the new `analyze_cancer_gene` rule lives here.
- `/home/leechuck/Public/software/hapli/benchmarks/hprc/README.md` and `/home/leechuck/Public/software/hapli/benchmarks/1000g_acmg/README.md` — runbook patterns to copy for the cancer + SV-vignette benchmarks.
- `/home/leechuck/Public/software/hapli/docs/paper_outline.md` — kept as the working paper; SV foregrounding edits land here through weeks 1-10.
