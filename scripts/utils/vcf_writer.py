from pathlib import Path
from typing import List, Dict, Any
import datetime
from .variant_models import Variant, SNV, Insertion, Deletion

class VCFWriter:
    def __init__(self, output_path: Path, sample_names: List[str], reference_path: Path):
        self.output_path = output_path
        self.sample_names = sample_names
        self.reference_path = reference_path

    def write_header(self, f, contigs: List[Dict[str, Any]]):
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={datetime.date.today().strftime('%Y%m%d')}\n")
        f.write(f"##reference=file://{self.reference_path.absolute()}\n")
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        
        for contig in contigs:
            f.write(f"##contig=<ID={contig['id']},length={contig['length']}>\n")
            
        header_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + self.sample_names
        f.write("\t".join(header_cols) + "\n")

    def format_variant(self, variant: Variant, genotype: str = "0|1") -> str:
        # Convert internal Variant model to VCF fields
        chrom = variant.chromosome
        pos = variant.position
        ref = variant.ref_allele
        alt = variant.alt_allele
        
        # Handle indels for VCF normalization (VCF requires one preceding base for indels)
        # Note: The variant models in scripts/utils/variant_models.py might store exact changes.
        # Standard VCF requires the base BEFORE the event to be included in REF and ALT.
        # However, looking at the previous generate_variants.py, it seemed to handle some of this.
        # For this prototype, I will assume the Variant models provide the exact string to put in VCF.
        # If not, I might need to adjust based on the 'variant_models.py' logic I read earlier.
        
        # Checking variant_models.py again:
        # SNV: ref/alt are single bases. OK.
        # Insertion: ref=base, alt=base+insert. This IS VCF compliant (anchor base included).
        # Deletion: ref=sequence, alt=first_base. This IS VCF compliant (anchor base included).
        
        info = "."
        if hasattr(variant, 'variant_type'):
            info = f"SVTYPE={variant.variant_type}"

        return "\t".join([
            chrom,
            str(pos),
            variant.variant_id or ".",
            ref,
            alt,
            "99", # QUAL
            "PASS", # FILTER
            info,
            "GT",
            genotype
        ])

    def write_variants(self, variants_by_sample: Dict[str, List[Variant]]):
        # This is a bit complex because VCF is position-sorted, not sample-sorted.
        # I need to merge all variants from all samples, sort them, and then write lines.
        # For simplicity, assuming 1 sample or strictly distinct variants for now?
        # No, "generate_test_data" is usually for 1 sample in my plan.
        
        # Flatten all variants: key = (chrom, pos, ref, alt) -> {sample: genotype}
        # But wait, we might have different variants at same position.
        
        # Simplified: We only support 1 sample for the initial generation script for simplicity.
        if len(self.sample_names) > 1:
            raise NotImplementedError("Multi-sample VCF writing not yet implemented in this simple writer.")
            
        sample = self.sample_names[0]
        variants = variants_by_sample[sample]
        
        # Sort by chrom (lexicographical unfortunately unless we have an order) and pos
        variants.sort(key=lambda v: (v.chromosome, v.position))
        
        with open(self.output_path, 'w') as f:
            # We need contig info. I'll need to pass it in or infer it.
            # For now, I'll cheat and infer from variants or just write minimal header
            # Actually, I'll require contigs to be passed to write_variants or init.
            # Let's assume write_header was called already.
            pass 

    def write(self, contigs: List[Dict[str, Any]], variants: List[Variant], sample_name: str, phasing: Dict[str, str] = None):
        """
        Write a single-sample VCF.
        phasing: dict mapping variant_id to genotype string (e.g. "0|1" or "1|0")
        """
        with open(self.output_path, 'w') as f:
            self.write_header(f, contigs)
            
            # Sort variants
            variants.sort(key=lambda v: (v.chromosome, v.position))
            
            for v in variants:
                gt = "0|1" # Default
                if phasing and v.variant_id in phasing:
                    gt = phasing[v.variant_id]
                
                f.write(self.format_variant(v, gt) + "\n")
