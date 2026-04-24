#!/usr/bin/env python3
import sys
import shutil
from pathlib import Path
import random
import logging

# Add the script directory to path to import utils
sys.path.append(str(Path(__file__).parent))

from utils.genome_utils import get_chromosome_lengths, generate_chromosome
from utils.variant_models import Variant, SNV, Insertion, Deletion
from utils.vcf_writer import VCFWriter
from utils.gff_generator import GFFGenerator, Gene
from utils.variant_models import HGVSFormatter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Configure Logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")

def main():
    output_dir = Path("data/test")
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)
    
    # 1. Generate Reference
    logging.info("Generating Reference Genome...")
    chroms = []
    # Create 1 main chromosome and 1 small one
    seq1 = generate_chromosome("chr1", 50000, 0.45) # 50kb
    seq2 = generate_chromosome("chr2", 10000, 0.45) # 10kb
    
    chroms = [seq1, seq2]
    
    ref_path = output_dir / "reference.fa"
    with open(ref_path, "w") as f:
        SeqIO.write(chroms, f, "fasta")
        
    # 2. Generate Genes
    logging.info("Generating Gene Annotations...")
    gff_gen = GFFGenerator(output_dir / "annotation.gff3")
    genes = []
    
    # Place a gene on chr1
    gene1 = gff_gen.generate_random_gene("chr1", (5000, 40000), "GENE", 1)
    genes.append(gene1)
    
    gff_gen.write(genes)
    
    # 3. Generate Variants
    logging.info("Generating Variants...")
    variants = []
    
    # Target specific regions to ensure impact
    # Variant 1: SNV in the first exon of the first transcript of Gene 1
    t1 = gene1.transcripts[0]
    exon1 = t1.exons[0]
    pos_snv = (exon1.start + exon1.end) // 2
    
    ref_base_snv = str(seq1.seq[pos_snv-1]) # 0-based index
    alt_base_snv = "T" if ref_base_snv != "T" else "A"
    
    v1 = SNV("chr1", pos_snv, ref_base_snv, alt_base_snv, "SNV")
    variants.append(v1)
    
    # Variant 2: Deletion in Intron (if exists) or Exon 2
    if len(t1.exons) > 1:
        # Try to put it in the intron between exon 1 and 2
        intron_start = t1.exons[0].end + 10
        intron_end = t1.exons[1].start - 10
        if intron_end > intron_start:
            pos_del = (intron_start + intron_end) // 2
            v2 = Deletion("chr1", pos_del, "A"*5, "A", "Deletion", deleted_length=5) # Dummy ref/alt for model, fix below
            # Fix ref allele from actual sequence
            real_ref = str(seq1.seq[pos_del-1:pos_del-1+5])
            v2.ref_allele = real_ref
            v2.alt_allele = real_ref[0]
            variants.append(v2)
            
    # Variant 3: Structural Variant (Large Deletion) covering last exon
    last_exon = t1.exons[-1]
    pos_sv = last_exon.start - 100
    len_sv = (last_exon.end - last_exon.start) + 200
    if pos_sv > 0:
        # Get Ref
        real_ref = str(seq1.seq[pos_sv-1:pos_sv-1+len_sv])
        v3 = Deletion("chr1", pos_sv, real_ref, real_ref[0], "Deletion", deleted_length=len_sv)
        variants.append(v3)

    # 4. Write VCF
    logging.info("Writing VCF...")
    vcf_path = output_dir / "variants.vcf"
    contig_info = [{"id": c.id, "length": len(c.seq)} for c in chroms]
    
    writer = VCFWriter(vcf_path, ["sample_001"], ref_path)
    
    # Phasing: 
    # v1 (SNV) -> Haplotype 1 (1|0)
    # v2 (Intron Del) -> Haplotype 2 (0|1)
    # v3 (Exon Del) -> Homozygous (1|1)
    
    phasing = {
        v1.variant_id: "1|0",
        v2.variant_id: "0|1",
        v3.variant_id: "1|1"
    }
    
    writer.write(contig_info, variants, "sample_001", phasing)
    
    # Index VCF
    import pysam
    logging.info("Indexing VCF...")
    pysam.tabix_index(str(vcf_path), preset="vcf", force=True)
    
    logging.info("Test Data Generation Complete!")
    logging.info(f"Directory: {output_dir}")
    logging.info(f"Files: reference.fa, annotation.gff3, variants.vcf.gz, variants.vcf.gz.tbi")

if __name__ == "__main__":
    main()
