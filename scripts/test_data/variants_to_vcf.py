#!/usr/bin/env python3
"""Convert variant JSON files to VCF format."""

import argparse
import json
import glob
import os
from datetime import datetime
from pathlib import Path
from pyfaidx import Fasta

def parse_arguments():
    parser = argparse.ArgumentParser(description="Convert variant descriptions to VCF")
    parser.add_argument("--variants-dir", required=True, help="Directory containing variant files")
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--output", required=True, help="Output VCF file")
    return parser.parse_args()

def load_variant_json(json_file):
    """Load variant data from JSON file."""
    with open(json_file, 'r') as f:
        return json.load(f)

def variant_to_vcf_record(variant, ref_fasta):
    """Convert variant object to VCF record."""
    chrom = variant['chromosome']
    pos = variant['position']
    var_type = variant['variant_type']
    
    # Get reference allele
    ref_seq = str(ref_fasta[chrom][pos-1:pos])  # VCF is 1-based
    
    # Build VCF fields based on variant type
    if var_type == 'SNV':
        ref = variant['ref_allele']
        alt = variant['alt_allele']
        info = "."
        vcf_pos = pos
    
    elif var_type == 'Insertion':
        # For insertions, we need the base before
        ref = str(ref_fasta[chrom][pos-1:pos])
        alt = ref + variant['insert_sequence']
        info = f"SVTYPE=INS;SVLEN={len(variant['insert_sequence'])}"
        vcf_pos = pos
    
    elif var_type == 'Deletion':
        # Include one base before deletion
        del_len = variant['deleted_length']
        ref = str(ref_fasta[chrom][pos-1:pos+del_len])
        alt = ref[0]
        info = f"SVTYPE=DEL;SVLEN=-{del_len};END={pos+del_len}"
        vcf_pos = pos
    
    elif var_type == 'Inversion':
        ref = "N"  # Symbolic allele
        alt = "<INV>"
        end = variant['end_position']
        info = f"SVTYPE=INV;END={end};SVLEN={end-pos+1}"
        vcf_pos = pos
    
    elif var_type == 'Duplication':
        ref = "N"
        alt = "<DUP>"
        end = variant['end_position']
        copies = variant.get('copy_number', 2)
        info = f"SVTYPE=DUP;END={end};SVLEN={end-pos+1};CN={copies}"
        vcf_pos = pos
    
    else:
        # Default for unknown types
        ref = "N"
        alt = "<UNK>"
        info = f"SVTYPE={var_type}"
        vcf_pos = pos
    
    # Build VCF line: CHROM POS ID REF ALT QUAL FILTER INFO
    vcf_id = variant.get('variant_id', '.')
    return f"{chrom}\t{vcf_pos}\t{vcf_id}\t{ref}\t{alt}\t.\tPASS\t{info}"

def main():
    args = parse_arguments()
    
    # Load reference
    ref_fasta = Fasta(args.reference)
    
    # Write VCF header
    with open(args.output, 'w') as vcf:
        # Write header
        vcf.write("##fileformat=VCFv4.3\n")
        vcf.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        vcf.write(f"##reference={args.reference}\n")
        
        # INFO fields
        vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        vcf.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
        vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
        vcf.write('##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number of segment containing variant">\n')
        
        # ALT fields
        vcf.write('##ALT=<ID=DEL,Description="Deletion">\n')
        vcf.write('##ALT=<ID=INS,Description="Insertion">\n')
        vcf.write('##ALT=<ID=DUP,Description="Duplication">\n')
        vcf.write('##ALT=<ID=INV,Description="Inversion">\n')
        
        # Contig lines
        for chr_name in ref_fasta.keys():
            length = len(ref_fasta[chr_name])
            vcf.write(f"##contig=<ID={chr_name},length={length}>\n")
        
        # Column header
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        # Process variants
        variant_files = sorted(glob.glob(os.path.join(args.variants_dir, "*_variant.json")))
        
        for variant_file in variant_files:
            try:
                variant = load_variant_json(variant_file)
                record = variant_to_vcf_record(variant, ref_fasta)
                vcf.write(record + "\n")
            except Exception as e:
                print(f"Error processing {variant_file}: {e}")

if __name__ == "__main__":
    main()
