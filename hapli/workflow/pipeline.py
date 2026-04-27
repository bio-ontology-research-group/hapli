"""
Hapli analysis pipeline — Phase 1 integration.

Stages per (sample, gene):

  1. Load gene + exon structure from the reference GFF3 via GFFProcessor.
  2. Materialise per-haplotype sequences over the gene region ± buffer via
     `bcftools consensus` (hapli/external/consensus.py).
  3. Per-transcript reference-cDNA-vs-haplotype splice alignment (retained
     for schema v1 compatibility — TUI + LLM consumers depend on it).
  4. [Phase 1 new] Per-haplotype gene-annotation lift-over with Liftoff →
     evidence.presence (intact / partial / deleted / duplicated / …).
  5. [Phase 1 new] Per-variant haplotype-aware consequence calling with
     `bcftools csq -p a` → evidence.consequence.
  6. Assemble the schema-v2 AnalysisResult and write JSON.

Phases 2+ will add per-haplotype protein extraction (using the lifted GFF3
from step 4), AlphaMissense / SpliceAI / LOFTEE / ESM2 scoring, and the
epistasis-residual computation, all under `evidence`.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import pysam

from ..alignment.aligner import SequenceAligner
from ..core.io import GFFProcessor, SequenceExtractor
from ..core.protein import (
    ProteinRecord,
    extract_proteins_for_gene,
    write_protein_fasta,
)
from ..core.protein_diff import diff_proteins
from ..core.schema import (
    AnalysisResult,
    ConsequenceCall,
    GeneEvidence,
    PresenceCall,
)
from ..external.alphamissense import AlphaMissenseLookup
from ..external.clinvar import ClinVarLookup
from ..external.consensus import consensus_region
from ..external.constraint import ConstraintLookup
from ..external.csq import BcftoolsNotAvailable, run_csq
from ..external.liftoff import LiftoffNotAvailable, run_liftoff
from ..interpretation.diploid import build_diploid_report


BUFFER_BP = 1000


class HapliPipeline:
    def __init__(
        self,
        ref_path: Path,
        vcf_path: Path | None,
        gff_path: Path,
        output_dir: Path,
        *,
        with_esm: bool = False,
        esm_model: str = "esm2_t6_8M_UR50D",
        gnomad_constraint: Path | None = None,
        clingen_dosage: Path | None = None,
        alphamissense_table: Path | None = None,
        clinvar_vcf: Path | None = None,
    ):
        """Mode A (VCF): pass `vcf_path`. Mode B (pre-assembled): pass None and
        supply hap FASTAs to `run_gene_assess`."""
        self.ref_path = Path(ref_path)
        self.vcf_path = Path(vcf_path) if vcf_path is not None else None
        self.gff_path = Path(gff_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.with_esm = with_esm
        self.esm_model = esm_model
        self.gnomad_constraint = gnomad_constraint
        self.clingen_dosage = clingen_dosage
        self.alphamissense_table = alphamissense_table
        self.clinvar_vcf = clinvar_vcf
        self.logger = logging.getLogger(__name__)
        self._esm_scorer = None                      # lazy-initialised
        self._constraint_lookup: ConstraintLookup | None = None
        self._alphamissense: AlphaMissenseLookup | None = None
        self._clinvar: ClinVarLookup | None = None

    # ─────────────────────────────────────────────────────────────
    # Entry point
    # ─────────────────────────────────────────────────────────────
    def run_gene_analysis(self, gene_name: str, sample_name: str) -> dict[str, Any]:
        self.logger.info("Starting analysis for gene %s, sample %s", gene_name, sample_name)

        gff, gene_feature = self._locate_gene(gene_name)
        region_str = (
            f"{gene_feature.seqid}:"
            f"{max(1, gene_feature.start - BUFFER_BP)}-"
            f"{gene_feature.end + BUFFER_BP}"
        )

        haplotypes = self._build_haplotypes(region_str, sample_name)
        self._write_haplotypes(haplotypes, sample_name, gene_name, gene_feature.seqid)

        transcripts_v1 = self._per_transcript_alignment(gff, gene_feature, haplotypes)

        evidence = GeneEvidence(gene=gene_name)
        lifted_gffs = self._run_liftoff(
            sample_name=sample_name,
            gene_name=gene_name,
            chrom=gene_feature.seqid,
            evidence=evidence,
        )
        self._run_csq(
            sample_name=sample_name,
            gene_name=gene_name,
            evidence=evidence,
        )
        self._extract_and_diff_proteins(
            sample_name=sample_name,
            gene_name=gene_name,
            lifted_gffs=lifted_gffs,
            evidence=evidence,
        )
        self._aggregate_alphamissense(evidence=evidence)
        self._annotate_clinvar(evidence=evidence)
        if self.with_esm:
            self._compute_epistasis(evidence=evidence)
        self._build_diploid(gene_name=gene_name, evidence=evidence)

        result = AnalysisResult(
            gene=gene_name,
            sample=sample_name,
            region=region_str,
            transcripts=transcripts_v1,
            evidence=evidence,
        )
        return self._write_result(result, sample_name, gene_name)

    # ─────────────────────────────────────────────────────────────
    # Mode B entry: pre-assembled haplotype FASTAs (HPRC / hifiasm / Verkko)
    # ─────────────────────────────────────────────────────────────
    def run_gene_assess(
        self,
        gene_name: str,
        sample_name: str,
        hap1_fasta: Path,
        hap2_fasta: Path,
    ) -> dict[str, Any]:
        """Run the downstream pipeline on pre-assembled haplotype FASTAs.

        No bcftools consensus, no bcftools csq — variant-equivalent information
        is recovered from the Liftoff-transferred GFF3 and the protein diff
        against the reference. Otherwise identical to Mode A end-to-end.
        """
        self.logger.info(
            "Starting Mode B assessment for gene %s, sample %s "
            "(hap1=%s, hap2=%s)",
            gene_name, sample_name, hap1_fasta, hap2_fasta,
        )
        gff, gene_feature = self._locate_gene(gene_name)
        region_str = (
            f"{gene_feature.seqid}:"
            f"{max(1, gene_feature.start - BUFFER_BP)}-"
            f"{gene_feature.end + BUFFER_BP}"
        )

        # Register the supplied haplotype FASTAs in the canonical per-hap layout.
        # Liftoff consumes these directly; we don't rewrite them.
        self._per_hap_fastas = {
            "hap1": Path(hap1_fasta),
            "hap2": Path(hap2_fasta),
        }
        for hap_name, p in self._per_hap_fastas.items():
            if not p.exists():
                raise FileNotFoundError(f"{hap_name} FASTA not found: {p}")
            fai = p.with_suffix(p.suffix + ".fai")
            if not fai.exists():
                self.logger.info("Indexing %s (fai missing)", p)
                pysam.faidx(str(p))

        # Skip v1 transcript alignment — it needs the raw hap sequences in
        # memory; for Mode B we rely on the Liftoff + protein-diff signals.
        transcripts_v1: list[dict[str, Any]] = []

        evidence = GeneEvidence(gene=gene_name)
        lifted_gffs = self._run_liftoff(
            sample_name=sample_name,
            gene_name=gene_name,
            chrom=gene_feature.seqid,
            evidence=evidence,
        )
        # Mode B skips bcftools csq: there is no VCF to annotate.
        self._extract_and_diff_proteins(
            sample_name=sample_name,
            gene_name=gene_name,
            lifted_gffs=lifted_gffs,
            evidence=evidence,
        )
        # AlphaMissense aggregation needs evidence.consequence which Mode B
        # doesn't populate; the aggregator no-ops gracefully on empty input.
        self._aggregate_alphamissense(evidence=evidence)
        self._annotate_clinvar(evidence=evidence)
        if self.with_esm:
            self._compute_epistasis(evidence=evidence)
        self._build_diploid(gene_name=gene_name, evidence=evidence)

        result = AnalysisResult(
            gene=gene_name,
            sample=sample_name,
            region=region_str,
            transcripts=transcripts_v1,
            evidence=evidence,
        )
        return self._write_result(result, sample_name, gene_name)

    # ─────────────────────────────────────────────────────────────
    # Stages
    # ─────────────────────────────────────────────────────────────
    def _locate_gene(self, gene_name: str):
        """Load the GFF and return (processor, gene_feature). GFFProcessor.__init__
        already logs the gene's coordinates, so we don't duplicate that here."""
        gff = GFFProcessor(self.gff_path, target_gene=gene_name)
        for f in gff.features_by_id.values():
            if f.featuretype == "gene":
                return gff, f
        raise ValueError(f"Gene {gene_name!r} not found in GFF.")

    def _build_haplotypes(self, region_str: str, sample_name: str) -> dict[str, str]:
        self.logger.info(
            "Materialising haplotypes for %s sample=%s via bcftools consensus",
            region_str, sample_name,
        )
        return consensus_region(
            reference_fasta=self.ref_path,
            vcf_path=self.vcf_path,
            sample=sample_name,
            region=region_str,
        )

    def _write_haplotypes(
        self,
        haplotypes: dict[str, str],
        sample_name: str,
        gene_name: str,
        chrom: str,
    ) -> None:
        """Write the combined (v1-style) FASTA *and* per-haplotype single-record
        FASTAs with the reference chromosome name as the header, ready for
        Liftoff to consume as target assemblies.
        """
        combined = self.output_dir / f"{sample_name}_{gene_name}_haplotypes.fa"
        with combined.open("w") as f:
            f.write(f">hap1\n{haplotypes['hap1']}\n")
            f.write(f">hap2\n{haplotypes['hap2']}\n")

        self._per_hap_fastas: dict[str, Path] = {}
        for hap_name in ("hap1", "hap2"):
            p = self.output_dir / f"{sample_name}_{gene_name}_{hap_name}.fa"
            with p.open("w") as f:
                f.write(f">{chrom}\n{haplotypes[hap_name]}\n")
            pysam.faidx(str(p))
            self._per_hap_fastas[hap_name] = p

    def _per_transcript_alignment(
        self,
        gff: GFFProcessor,
        gene_feature,
        haplotypes: dict[str, str],
    ) -> list[dict[str, Any]]:
        """Schema-v1 transcript-vs-haplotype splice alignment.

        Retained verbatim from the pre-Phase-1 pipeline so the TUI and LLM
        consumers keep working and the showboat walkthrough keeps verifying.
        Phase 2 will supersede this with protein-level extraction from the
        Liftoff-produced per-haplotype GFF3.
        """
        seq_extractor = SequenceExtractor(self.ref_path)
        aligner = SequenceAligner()

        transcripts: list[dict[str, Any]] = []
        for transcript in gff.get_children(gene_feature.id):
            # GENCODE/Ensembl use "transcript"; legacy / other GFFs use "mRNA".
            if transcript.featuretype not in ("mRNA", "transcript"):
                continue
            exons = [e for e in gff.get_children(transcript.id) if e.featuretype == "exon"]
            exons.sort(key=lambda e: e.start)

            # Correct reverse-strand splicing: extract each exon on + strand, concat, RC if -.
            raw_parts = [
                seq_extractor.fasta.fetch(e.seqid, e.start - 1, e.end)
                for e in exons
            ]
            full_tx_seq = "".join(raw_parts)
            if transcript.strand == "-":
                full_tx_seq = seq_extractor._reverse_complement(full_tx_seq)

            t_res: dict[str, Any] = {"id": transcript.id, "alignments": {}}
            for hap_name, hap_seq in haplotypes.items():
                alignments = aligner.align(full_tx_seq, transcript.id, hap_seq, hap_name)
                if not alignments:
                    t_res["alignments"][hap_name] = "Unmapped"
                    continue
                best = alignments[0]
                nm = best.get_tag("NM") if best.has_tag("NM") else -1
                length = len(full_tx_seq)
                identity = 1.0 - (nm / length) if length > 0 else 0
                t_res["alignments"][hap_name] = {
                    "identity": identity,
                    "nm": nm,
                    "cigar": best.cigarstring,
                    "is_perfect": nm == 0,
                }
            transcripts.append(t_res)
        return transcripts

    # ─────────────────────────────────────────────────────────────
    # Phase 1 new stages: Liftoff + csq
    # ─────────────────────────────────────────────────────────────
    def _run_liftoff(
        self,
        sample_name: str,
        gene_name: str,
        chrom: str,
        evidence: GeneEvidence,
    ) -> dict[str, Path]:
        """Lift the reference GFF3 onto each haplotype. Populate evidence.presence.
        Returns {hap_name: path_to_lifted_gff} so downstream stages can slice
        proteins out of the haplotype-native coordinates.
        """
        from ..external.liftoff import parse_lifted_gff
        lifted: dict[str, Path] = {}
        try:
            for hap_name, hap_fa in self._per_hap_fastas.items():
                # Output paths are per (sample, hap), NOT per (sample, gene, hap),
                # so the lifted GFF is reused across all genes for this sample.
                # First gene call pays the Liftoff cost; remaining genes for the
                # same sample short-circuit and just re-parse the existing GFF.
                out_gff = self.output_dir / f"{sample_name}_{hap_name}.lifted.gff3"
                polished_gff = out_gff.parent / (out_gff.name + "_polished")
                unmapped = self.output_dir / f"{sample_name}_{hap_name}.unmapped.txt"
                tmp_dir = self.output_dir / f"{sample_name}_{hap_name}.liftoff_tmp"

                presence = None
                # Prefer the polished GFF when present — it carries Liftoff's
                # CDS-validation pass that re-classifies marginal lifts as
                # truly partial vs. intact. Falling back to the un-polished
                # GFF would inflate "partial" calls and tank scores spuriously.
                cache_candidates = [p for p in (polished_gff, out_gff) if p.exists() and p.stat().st_size > 0]
                for cache_path in cache_candidates:
                    candidate = parse_lifted_gff(cache_path)
                    if gene_name in candidate:
                        presence = candidate
                        self.logger.info(
                            "Liftoff %s: reusing cached lifted GFF (%s)",
                            hap_name, cache_path.name,
                        )
                        break
                if presence is None:
                    try:
                        result = run_liftoff(
                            haplotype_fasta=hap_fa,
                            reference_fasta=self.ref_path,
                            reference_gff=self.gff_path,
                            out_gff=out_gff,
                            unmapped_txt=unmapped,
                            intermediate_dir=tmp_dir,
                            polish=True,
                            logger=self.logger,
                        )
                    except RuntimeError as exc:
                        self.logger.warning("Liftoff failed on %s: %s", hap_name, exc)
                        evidence.presence[hap_name] = PresenceCall(status="not_run", source="liftoff")
                        continue
                    presence = result.presence

                call = presence.get(gene_name)
                if call is None:
                    call = PresenceCall(status="uncertain", source="liftoff")
                evidence.presence[hap_name] = call
                lifted[hap_name] = out_gff
                self.logger.info("Liftoff %s: gene %s status=%s", hap_name, gene_name, call.status)
        except LiftoffNotAvailable as exc:
            self.logger.warning("Skipping Liftoff: %s", exc)
            for hap_name in ("hap1", "hap2"):
                evidence.presence[hap_name] = PresenceCall(status="not_run", source="liftoff")
        return lifted

    def _run_csq(self, sample_name: str, gene_name: str, evidence: GeneEvidence) -> None:
        """Run bcftools csq and pull consequences relevant to this gene's transcripts."""
        out_vcf = self.output_dir / f"{sample_name}_{gene_name}.csq.vcf.gz"
        try:
            result = run_csq(
                vcf_path=self.vcf_path,
                reference_fasta=self.ref_path,
                gff3=self.gff_path,
                sample=sample_name,
                out_vcf=out_vcf,
                logger=self.logger,
            )
        except BcftoolsNotAvailable as exc:                 # pragma: no cover
            self.logger.warning("Skipping bcftools csq: %s", exc)
            return
        except RuntimeError as exc:
            self.logger.warning("bcftools csq failed: %s", exc)
            return

        # Consequences are per-variant; we retain all — downstream consumers
        # (TUI / LLM / diploid report) filter by transcript themselves.
        evidence.consequence.extend(result.consequences)
        self.logger.info("csq: %d consequence records", len(result.consequences))

    def _extract_and_diff_proteins(
        self,
        sample_name: str,
        gene_name: str,
        lifted_gffs: dict[str, Path],
        evidence: GeneEvidence,
    ) -> None:
        """Phase 2: translate each haplotype's copy of the gene and diff it
        against the reference protein. Populates evidence.protein and writes
        per-haplotype + reference protein FASTAs as first-class artifacts.
        """
        # Reference proteins (one per transcript) — extracted once from the
        # reference FASTA + original GFF3.
        try:
            ref_proteins = extract_proteins_for_gene(
                self.ref_path, self.gff_path, gene_name, logger=self.logger
            )
        except Exception as exc:
            self.logger.warning("Could not extract reference proteins for %s: %s", gene_name, exc)
            return
        if not ref_proteins:
            self.logger.info("No reference protein records for %s — skipping protein stage", gene_name)
            return

        ref_fa_out = self.output_dir / f"{sample_name}_{gene_name}_ref.pep.fa"
        write_protein_fasta(ref_proteins.values(), ref_fa_out, header_prefix="ref.")

        # Per-haplotype proteins, only on haplotypes where Liftoff returned a GFF.
        # (deleted haplotypes skip extraction — they're flagged in evidence.presence).
        for hap_name in ("hap1", "hap2"):
            presence = evidence.presence.get(hap_name)
            if presence is not None and presence.status == "deleted":
                self.logger.info(
                    "%s: gene %s deleted — skipping protein extraction", hap_name, gene_name
                )
                continue
            lifted_gff = lifted_gffs.get(hap_name)
            if lifted_gff is None or not lifted_gff.exists():
                self.logger.info(
                    "%s: no lifted GFF — skipping protein extraction", hap_name
                )
                continue
            hap_fa = self._per_hap_fastas[hap_name]
            try:
                hap_proteins = extract_proteins_for_gene(
                    hap_fa, lifted_gff, gene_name, logger=self.logger
                )
            except Exception as exc:
                self.logger.warning(
                    "%s: protein extraction failed: %s", hap_name, exc
                )
                continue
            hap_fa_out = self.output_dir / f"{sample_name}_{gene_name}_{hap_name}.pep.fa"
            write_protein_fasta(hap_proteins.values(), hap_fa_out, header_prefix=f"{hap_name}.")

            hap_idx = 1 if hap_name == "hap1" else 2
            for tx_id, hap_protein in hap_proteins.items():
                ref_protein = ref_proteins.get(tx_id)
                if ref_protein is None:
                    continue
                diff = diff_proteins(
                    transcript_id=tx_id,
                    haplotype=hap_idx,
                    ref_seq=ref_protein.protein,
                    hap_seq=hap_protein.protein,
                    logger=self.logger,
                )
                evidence.protein.append(diff)
                self.logger.info(
                    "%s %s: protein identity=%.4f subs=%d pts=%s fs_region=%s",
                    hap_name, tx_id, diff.identity, len(diff.substitutions),
                    diff.premature_stop_at,
                    "yes" if diff.frameshift_region else "no",
                )

    def _compute_epistasis(self, evidence: GeneEvidence) -> None:
        """Phase 3: for each haplotype protein-diff with ≥2 substitutions,
        compute the ESM2 additive-vs-joint residual. Populates evidence.epistasis.
        """
        # Lazy import so the core pipeline can run without torch installed.
        try:
            from ..interpretation.esm_scoring import ESM2Scorer, EsmNotAvailable
            from ..interpretation.epistasis import compute_residuals_for_diff
        except ImportError as exc:                   # pragma: no cover
            self.logger.warning("Epistasis dependencies missing: %s", exc)
            return

        if self._esm_scorer is None:
            try:
                self._esm_scorer = ESM2Scorer(model_name=self.esm_model)
            except EsmNotAvailable as exc:
                self.logger.warning(
                    "Skipping epistasis (ESM2 not available): %s", exc
                )
                return

        # Cross-reference each ProteinDiff to a (ref_protein, hap_protein) pair
        # by re-extracting from the artifact FASTAs written in Phase 2. This
        # decouples the epistasis step from the in-memory protein objects.
        ref_by_tx = self._read_protein_fasta(
            next(self.output_dir.glob("*_ref.pep.fa"), None)
        )
        hap_by_tx: dict[int, dict[str, str]] = {}
        for hap_idx, hap_name in ((1, "hap1"), (2, "hap2")):
            hap_fa = next(self.output_dir.glob(f"*_{hap_name}.pep.fa"), None)
            if hap_fa is not None:
                hap_by_tx[hap_idx] = self._read_protein_fasta(hap_fa)

        for diff in evidence.protein:
            hap_idx = diff.haplotype
            tx = diff.transcript
            ref_seq = ref_by_tx.get(tx)
            hap_seq = hap_by_tx.get(hap_idx, {}).get(tx)
            if ref_seq is None or hap_seq is None:
                continue
            residual = compute_residuals_for_diff(
                self._esm_scorer,
                transcript_id=tx,
                haplotype=hap_idx,
                ref_seq=ref_seq,
                hap_seq=hap_seq,
                logger=self.logger,
            )
            evidence.epistasis.append(residual)
            if residual.n_variants >= 2:
                self.logger.info(
                    "%s hap%d: epistasis n=%d S_add=%+.2f S_joint=%+.2f "
                    "residual=%+.2f flagged=%s",
                    tx, hap_idx, residual.n_variants,
                    residual.s_additive, residual.s_joint,
                    residual.residual, residual.flagged,
                )

    def _aggregate_alphamissense(self, evidence: GeneEvidence) -> None:
        """Look up an AlphaMissense score for every missense ConsequenceCall in
        evidence.consequence, then aggregate per (transcript, haplotype) into
        evidence.missense_agg. No-op if no AM table was provided or if the
        evidence has no missense records.
        """
        if self.alphamissense_table is None:
            return
        if not evidence.consequence:
            return
        if self._alphamissense is None:
            try:
                self._alphamissense = AlphaMissenseLookup(
                    self.alphamissense_table, logger=self.logger,
                )
            except FileNotFoundError as exc:
                self.logger.warning("AlphaMissense table missing: %s", exc)
                return
        # Bucket: (transcript, haplotype) → list of am scores.
        from collections import defaultdict
        from ..core.schema import MissenseAggregate
        bucket: dict[tuple[str, int], list[float]] = defaultdict(list)
        for c in evidence.consequence:
            if c.consequence != "missense":
                continue
            # AlphaMissense scores per-position SNVs only. Decompose MNVs
            # (ref="GCT" alt="TGG") into per-position substitutions, look each
            # up, and take the MAX (most-pathogenic-component aggregation).
            if len(c.ref) != len(c.alt) or not c.ref or not c.alt:
                # Indels and frameshifts don't have AM scores; skip.
                continue
            scores: list[float] = []
            for offset, (rb, ab) in enumerate(zip(c.ref, c.alt)):
                if rb == ab:
                    continue
                hit = self._alphamissense.get(c.chrom, c.pos + offset, rb, ab)
                if hit is not None:
                    scores.append(hit.score)
            if scores:
                # Per-MNV: keep the worst-component score for the aggregation.
                bucket[(c.transcript, c.haplotype)].append(max(scores))
        for (tx, hap), scores in bucket.items():
            evidence.missense_agg.append(MissenseAggregate(
                transcript=tx, haplotype=hap,
                n_missense=len(scores),
                max_am=max(scores),
                mean_am=sum(scores) / len(scores),
                sum_am=sum(scores),
                source="alphamissense_v1",
            ))
        if bucket:
            self.logger.info(
                "alphamissense: aggregated %d missense scores into %d (transcript,hap) buckets",
                sum(len(s) for s in bucket.values()), len(bucket),
            )

    def _annotate_clinvar(self, evidence: GeneEvidence) -> None:
        """Stamp each ConsequenceCall (or per-component SNV of an MNV) with the
        ClinVar `clnsig` label, when a ClinVar VCF is provided. Annotation is
        attached as an extra `clnsig` attribute on the ConsequenceCall instance
        — the pydantic schema is extra=allow, so it round-trips into the JSON.
        """
        if self.clinvar_vcf is None:
            return
        if not evidence.consequence:
            return
        if self._clinvar is None:
            try:
                self._clinvar = ClinVarLookup.from_vcf(self.clinvar_vcf, logger=self.logger)
            except FileNotFoundError as exc:
                self.logger.warning("ClinVar VCF missing: %s", exc)
                return
        n_annotated = 0
        for c in evidence.consequence:
            # Decompose MNVs into per-position SNVs (same pattern as AM lookup).
            if not c.ref or not c.alt or len(c.ref) != len(c.alt):
                continue
            best_clnsig: str | None = None
            best_priority = 999
            priority = {
                "pathogenic": 0, "likely_pathogenic": 1, "vus": 2,
                "conflicting": 3, "drug_response": 4, "association": 5,
                "risk_factor": 6, "likely_benign": 7, "benign": 8,
            }
            for offset, (rb, ab) in enumerate(zip(c.ref, c.alt)):
                if rb == ab:
                    continue
                hit = self._clinvar.get(c.chrom, c.pos + offset, rb, ab)
                if hit is None or hit.clnsig is None:
                    continue
                p = priority.get(hit.clnsig, 99)
                if p < best_priority:
                    best_priority = p
                    best_clnsig = hit.clnsig
            if best_clnsig is not None:
                c.clnsig = best_clnsig         # type: ignore[attr-defined]
                n_annotated += 1
        if n_annotated:
            self.logger.info(
                "clinvar: annotated %d / %d consequence records",
                n_annotated, len(evidence.consequence),
            )

    def _build_diploid(self, gene_name: str, evidence: GeneEvidence) -> None:
        """Phase 4: build the two-number diploid report + constraint priors."""
        if self._constraint_lookup is None:
            self._constraint_lookup = ConstraintLookup.load(
                gnomad_path=self.gnomad_constraint,
                clingen_path=self.clingen_dosage,
                logger=self.logger,
            )
        priors = self._constraint_lookup.get(gene_name)
        evidence.diploid = build_diploid_report(
            evidence, constraints=priors, logger=self.logger,
        )

    @staticmethod
    def _read_protein_fasta(path: Path | None) -> dict[str, str]:
        """Return {transcript_id: protein_sequence}. Header format written by
        `core.protein.write_protein_fasta` is `>prefix.transcript_id ...`
        — we key on the bare transcript ID."""
        out: dict[str, str] = {}
        if path is None or not path.exists():
            return out
        cur_id: str | None = None
        parts: list[str] = []
        with path.open() as f:
            for line in f:
                if line.startswith(">"):
                    if cur_id is not None:
                        out[cur_id] = "".join(parts).rstrip("*")
                    header = line[1:].strip().split()[0]
                    # header looks like 'ref.T1' or 'hap1.T1' → strip prefix
                    cur_id = header.split(".", 1)[-1] if "." in header else header
                    parts = []
                else:
                    parts.append(line.strip())
            if cur_id is not None:
                out[cur_id] = "".join(parts).rstrip("*")
        return out

    # ─────────────────────────────────────────────────────────────
    # Output
    # ─────────────────────────────────────────────────────────────
    def _write_result(
        self,
        result: AnalysisResult,
        sample_name: str,
        gene_name: str,
    ) -> dict[str, Any]:
        json_out = self.output_dir / f"{sample_name}_{gene_name}_analysis.json"

        # Dump to a dict and strip empty evidence sub-fields so the v1 shape
        # on the surface is identical to the pre-Phase-1 output for fixtures
        # that don't trigger any new stage (keeps showboat + TUI happy).
        dumped = result.model_dump()
        # schema_version sits outside the v1 keys; keep it at the top level.
        payload: dict[str, Any] = {
            "schema_version": dumped.pop("schema_version"),
            "gene": dumped["gene"],
            "sample": dumped["sample"],
            "region": dumped["region"],
            "transcripts": dumped["transcripts"],
            "evidence": dumped["evidence"],
        }

        with json_out.open("w") as f:
            json.dump(payload, f, indent=2)
        self.logger.info("Analysis saved to %s", json_out)
        return payload
