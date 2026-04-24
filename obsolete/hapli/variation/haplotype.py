import pysam
import logging
from pathlib import Path
from typing import List, Tuple, Dict
import os

class HaplotypeGenerator:
    """Generates haplotype sequences from Reference + VCF."""
    
    def __init__(self, reference_path: Path, vcf_path: Path):
        self.reference_path = reference_path
        self.vcf_path = vcf_path
        self.ref_fasta = pysam.FastaFile(str(reference_path))
        self.vcf = pysam.VariantFile(str(vcf_path))

    def generate_haplotypes_for_region(self, chrom: str, start: int, end: int, sample: str, output_file: Path):
        """
        Generates two haplotype sequences for a specific region and sample.
        Writes them to output_file as FASTA.
        """
        logging.info(f"Generating haplotypes for {sample} in {chrom}:{start}-{end}")
        
        # Fetch reference sequence
        # Expand region slightly to ensure we catch variants at boundaries if needed, 
        # but for now let's stick to the requested region.
        ref_seq = self.ref_fasta.fetch(chrom, start, end)
        
        # Fetch variants
        variants = list(self.vcf.fetch(chrom, start, end))
        
        # We need to construct two sequences
        # This is a complex operation if there are overlapping variants or indels.
        # Simple approach: iterate variants and apply them.
        
        # Pysam variants have samples[sample]['GT'] -> (0, 1) or (1, 1) etc.
        
        hap1_seq = list(ref_seq)
        hap2_seq = list(ref_seq)
        
        # Offset tracking is crucial because Indels shift coordinates.
        # However, applying variants to a string directly changes indices for subsequent variants.
        # Better approach: Build chunks.
        
        # But wait, variants in VCF are sorted by position.
        # We can iterate from end to start to avoid index shifting!
        
        variants.sort(key=lambda x: x.pos, reverse=True)
        
        for record in variants:
            # VCF pos is 1-based.
            # Record.pos is 1-based start.
            # Record.ref is reference allele.
            # Record.alts is tuple of alt alleles.
            
            # Check overlap with our region
            # We fetched with start, end (0-based in pysam fetch? No, fetch uses region string or start/end)
            # pysam.FastaFile.fetch(start, end) uses 0-based.
            # pysam.VariantFile.fetch(start, end) uses 0-based.
            
            # Relative position in our ref_seq
            # global pos P -> relative pos p = P - start - 1 (because P is 1-based)
            # Wait, pysam record.pos is 1-based. 
            # ref_seq starts at `start` (0-based).
            # So index in ref_seq is `record.pos - 1 - start`.
            
            rel_pos = record.pos - 1 - start
            if rel_pos < 0 or rel_pos >= len(ref_seq):
                continue
                
            genotype = record.samples[sample]['GT']
            if not genotype: continue
            
            # genotype is e.g. (0, 1). 
            # 0 = ref, 1 = alt1.
            
            # Apply to Hap1
            self._apply_variant(hap1_seq, rel_pos, record.ref, record.alts, genotype[0])
            
            # Apply to Hap2
            # Handle haploid/diploid
            gt2 = genotype[1] if len(genotype) > 1 else genotype[0]
            self._apply_variant(hap2_seq, rel_pos, record.ref, record.alts, gt2)
            
        # Join
        h1 = "".join(hap1_seq)
        h2 = "".join(hap2_seq)
        
        with open(output_file, 'a') as f:
            f.write(f">{sample}_1_{chrom}:{start}-{end}\n{h1}\n")
            f.write(f">{sample}_2_{chrom}:{start}-{end}\n{h2}\n")

    def _apply_variant(self, seq_list: List[str], rel_pos: int, ref: str, alts: Tuple[str], allele_idx: int):
        """
        Apply variant to sequence list in place.
        Since we iterate backwards, we can replace chunks without affecting future indices (which are smaller).
        """
        if allele_idx == 0:
            return # Ref allele, no change
            
        if allele_idx is None:
             return # Missing call
             
        alt = alts[allele_idx - 1]
        
        # Check if ref matches (sanity check)
        # Note: ref in VCF might be longer than 1 base for Indels.
        # seq_list[rel_pos : rel_pos + len(ref)] should match ref.
        
        # We replace seq_list[rel_pos : rel_pos + len(ref)] with alt.
        # But seq_list is a list of characters.
        
        # Verify ref match roughly (optional but good)
        # current_ref = "".join(seq_list[rel_pos : rel_pos + len(ref)])
        # if current_ref != ref: 
        #    logging.warning(f"Ref mismatch at {rel_pos}: expected {ref}, got {current_ref}")
        
        # Substitution
        seq_list[rel_pos : rel_pos + len(ref)] = list(alt)
