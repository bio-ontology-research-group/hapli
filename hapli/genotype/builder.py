"""
Haplotype construction from a phased VCF + reference.

As of the Phase 0 refactor, this module is a thin adapter over
`hapli.external.consensus.consensus_region` (samtools faidx |
bcftools consensus). The previous hand-rolled `_apply_variants`
was removed because it:

  1. Silently corrupted symbolic SV records by inserting the
     literal string "<DEL>" into the haplotype FASTA, and
  2. Misindexed pysam's `sample_rec.alleles` (genotype-ordered)
     with `allele_indices` (VCF-ordered), so SNVs on hap1 of a
     phased GT=1|0 record were not applied.

See `docs/walkthrough.md` for the updated end-to-end baseline.
"""

from __future__ import annotations

import logging
from pathlib import Path

from ..external.consensus import consensus_region


class HaplotypeBuilder:
    def __init__(self, reference_path: Path, vcf_path: Path):
        self.reference_path = Path(reference_path)
        self.vcf_path = Path(vcf_path)
        self.logger = logging.getLogger(__name__)

    def build_haplotypes(self, region: str, sample: str) -> dict[str, str]:
        """Return {'hap1': seq, 'hap2': seq} for the 1-based inclusive region."""
        self.logger.info("Materialising haplotypes for %s sample=%s via bcftools consensus", region, sample)
        return consensus_region(
            reference_fasta=self.reference_path,
            vcf_path=self.vcf_path,
            sample=sample,
            region=region,
        )
