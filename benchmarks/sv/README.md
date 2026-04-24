# Phase 5-D — Structural-variant coverage benchmark

MaveDB scoresets are amino-acid DMS data — they do not contain SVs. This
benchmark validates that hapli's haplotype-level scoring correctly classifies
SV-shaped haplotypes (whole-gene deletions, large in-frame deletions,
frameshift truncations at varying positions) as LoF, end-to-end through the
unified scorer.

The SVs here are expressed in HGVS protein notation (`p.Met1_Val50del`,
`p.Arg100Ter`, `p.Leu300GlyfsTer10`) — which is exactly the shape the
upstream VCF→protein pipeline would emit when bcftools consensus has been
applied to a symbolic `<DEL>` or `<INV>` in the VCF. The pipeline converts
symbolic variants to protein-level consequences at the Liftoff step;
this benchmark exercises the downstream scoring in isolation.

## What the benchmark asserts

For each synthetic SV-bearing haplotype:

1. `category ∈ {lof, indel}` — the shape router correctly classifies the event.
2. `length_changed = True` — downstream consumers know the protein diverged structurally.
3. `premature_stop_at` is populated for stop-gain / frameshift cases.
4. `hapli_score_heuristic ≤ 0.5` — a haplotype with a truncating lesion
   scores low enough to trigger the Phase 4 compound-het LoF detector.
5. `residual is None` — ESM residual is not defined for length-changing
   haplotypes; the shape signals carry the information instead.

## Layout

```
cases.py      — Synthetic SV-bearing haplotype definitions
run.py        — Score each case and report PASS/FAIL per assertion
```

## Quickstart

```bash
./run.py --out results/
```

Expected: all cases PASS on the tiny 50-residue synthetic protein fixture.
```
