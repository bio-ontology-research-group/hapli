#!/usr/bin/env python3
"""
Convert multi-sample haplotype calls to diploid genotypes.

This script takes a multi-sample VCF with haplotype columns and converts
it to a multi-sample VCF with proper diploid genotypes for each sample.
"""

import sys
import re

def convert_haplotype_to_genotype(hap1_gt, hap2_gt):
    """Convert haplotype genotypes to diploid genotype"""
    # Handle missing data
    if hap1_gt == '.' or hap2_gt == '.':
        return './.'
    
    # Convert to integers
    try:
        h1 = int(hap1_gt)
        h2 = int(hap2_gt)
        return str(h1) + "|" + str(h2)  # Phased genotype
    except ValueError:
        return './.'

def get_reference_chromosome_name(reference_file):
    """Extract the first chromosome name from reference FASTA"""
    try:
        with open(reference_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Extract chromosome name (first word after >)
                    chrom_name = line[1:].split()[0]
                    return chrom_name
    except:
        pass
    return None

def main():
    if len(sys.argv) < 4:
        print("Usage: convert_multi_sample_haplotypes.py input_vcf output_vcf sample1 [sample2 ...]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    sample_names = sys.argv[3:]

    # Try to get reference file path from VCF header to get correct chromosome name
    reference_file = None
    reference_chrom = None
    
    # Read the input VCF and find haplotype columns
    with open(input_file, 'r') as f_in:
        header_line = None
        for line in f_in:
            if line.startswith('##reference=file://'):
                # Extract reference file path
                ref_path = line.strip().split('file://', 1)[1]
                # Remove any /data prefix that might be from Docker mount
                if ref_path.startswith('/data/'):
                    ref_path = ref_path[6:]
                reference_chrom = get_reference_chromosome_name(ref_path)
                print(f"Found reference file: {ref_path}", file=sys.stderr)
                print(f"Found reference chromosome: {reference_chrom}", file=sys.stderr)
            elif line.startswith('#CHROM'):
                header_line = line.strip()
                break
    
    if not header_line:
        print("Error: Could not find header line in VCF")
        sys.exit(1)
    
    # Parse header to find column positions
    header_parts = header_line.split('\t')
    sample_hap_positions = {}
    
    for sample in sample_names:
        hap1_col = f"{sample}_hap1"
        hap2_col = f"{sample}_hap2"
        
        hap1_pos = None
        hap2_pos = None
        
        for i, col in enumerate(header_parts):
            if col == hap1_col:
                hap1_pos = i
            elif col == hap2_col:
                hap2_pos = i
        
        sample_hap_positions[sample] = (hap1_pos, hap2_pos)
        print(f"Sample {sample}: hap1 at position {hap1_pos}, hap2 at position {hap2_pos}", file=sys.stderr)

    # Process the VCF
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('##'):
                # Copy header lines, but update contig if we have reference chromosome
                if line.startswith('##contig=<ID=reference') and reference_chrom:
                    # Replace reference with actual chromosome name
                    line = line.replace('ID=reference', f'ID={reference_chrom}')
                f_out.write(line)
            elif line.startswith('#CHROM'):
                # Create new header with sample names
                parts = line.strip().split('\t')
                new_header = parts[:9] + sample_names  # Keep first 9 columns, add sample names
                f_out.write('\t'.join(new_header) + '\n')
            else:
                # Process variant lines
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                
                # Update chromosome name if we have reference chromosome
                if reference_chrom and parts[0] == 'reference':
                    parts[0] = reference_chrom
                
                # Build new line with diploid genotypes
                new_parts = parts[:9]  # Keep first 9 columns (CHROM through FORMAT)
                
                # Check if any sample has a variant (not 0|0)
                has_variant = False
                sample_genotypes = []
                
                for sample in sample_names:
                    hap1_pos, hap2_pos = sample_hap_positions[sample]
                    
                    # Get haplotype genotypes
                    hap1_gt = parts[hap1_pos] if hap1_pos and hap1_pos < len(parts) else '.'
                    hap2_gt = parts[hap2_pos] if hap2_pos and hap2_pos < len(parts) else '.'
                    
                    # Convert to diploid genotype
                    diploid_gt = convert_haplotype_to_genotype(hap1_gt, hap2_gt)
                    sample_genotypes.append(diploid_gt)
                    
                    # Check if this sample has a variant
                    if diploid_gt != '0|0' and diploid_gt != './.':
                        has_variant = True
                
                # Only output variants where at least one sample has a non-reference genotype
                if has_variant:
                    new_parts.extend(sample_genotypes)
                    f_out.write('\t'.join(new_parts) + '\n')

if __name__ == "__main__":
    main()
