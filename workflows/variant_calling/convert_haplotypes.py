#!/usr/bin/env python3
"""
Convert haplotype calls to diploid genotypes.

This script takes a multi-sample VCF with haplotype columns and converts
it to a single-sample VCF with proper diploid genotypes.
"""

import sys
import os

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

def get_all_reference_chromosome_names(reference_file):
    """Extract all chromosome names from reference FASTA"""
    chromosomes = []
    try:
        # Handle both absolute and relative paths
        if not os.path.exists(reference_file):
            # Try without /data prefix
            alt_path = reference_file.replace('/data/', '')
            if os.path.exists(alt_path):
                reference_file = alt_path
            else:
                print(f"Warning: Reference file not found: {reference_file}", file=sys.stderr)
                return []
                
        with open(reference_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Extract chromosome name (first word after >)
                    chrom_name = line[1:].split()[0]
                    chromosomes.append(chrom_name)
                    print(f"Found chromosome: {chrom_name}", file=sys.stderr)
    except Exception as e:
        print(f"Error reading reference file {reference_file}: {e}", file=sys.stderr)
    return chromosomes

def main():
    if len(sys.argv) != 6:
        print("Usage: convert_haplotypes.py input_vcf output_vcf sample_name hap1_pos hap2_pos")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    sample_name = sys.argv[3]
    hap1_pos = int(sys.argv[4]) if sys.argv[4] and sys.argv[4] != 'None' else None
    hap2_pos = int(sys.argv[5]) if sys.argv[5] and sys.argv[5] != 'None' else None

    print(f"Converting haplotypes for sample: {sample_name}", file=sys.stderr)
    print(f"Hap1 position: {hap1_pos}, Hap2 position: {hap2_pos}", file=sys.stderr)

    # Try to get reference file path from VCF header to get correct chromosome names
    reference_chroms = []
    
    with open(input_file, 'r') as f_in:
        for line in f_in:
            if line.startswith('##reference=file://'):
                # Extract reference file path
                ref_path = line.strip().split('file://', 1)[1]
                # Remove any /data prefix that might be from Docker mount
                if ref_path.startswith('/data/'):
                    ref_path = ref_path[6:]
                reference_chroms = get_all_reference_chromosome_names(ref_path)
                print(f"Found reference file: {ref_path}", file=sys.stderr)
                print(f"Found reference chromosomes: {reference_chroms}", file=sys.stderr)
                break
            elif line.startswith('#CHROM'):
                break

    # If we couldn't get chromosomes from reference, use defaults
    if not reference_chroms:
        reference_chroms = ["chr1", "chr2", "chr3"]  # Add more as needed
        print(f"Using default chromosome names: {reference_chroms}", file=sys.stderr)

    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('##'):
                # Copy header lines, but update contig if we have reference chromosomes
                if line.startswith('##contig=<ID=reference') and reference_chroms:
                    # Replace reference with actual chromosome names - add all chromosomes
                    for chrom in reference_chroms:
                        new_line = line.replace('ID=reference', f'ID={chrom}')
                        f_out.write(new_line)
                    continue  # Skip the original line
                f_out.write(line)
            elif line.startswith('#CHROM'):
                # Modify column header
                parts = line.strip().split('\t')
                new_header = parts[:9] + [sample_name]  # Keep first 9 columns, replace samples
                f_out.write('\t'.join(new_header) + '\n')
            else:
                # Process variant lines
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                
                # Update chromosome name if we have reference chromosomes
                # Map "reference" to the appropriate chromosome based on position or other logic
                if reference_chroms and parts[0] == 'reference':
                    # For now, just use the first chromosome - this might need more sophisticated logic
                    parts[0] = reference_chroms[0]  # This is a simplification
                    
                # Get haplotype genotypes
                hap1_gt = parts[hap1_pos - 1] if hap1_pos and hap1_pos <= len(parts) else '.'
                hap2_gt = parts[hap2_pos - 1] if hap2_pos and hap2_pos <= len(parts) else '.'
                
                # Convert to diploid genotype
                diploid_gt = convert_haplotype_to_genotype(hap1_gt, hap2_gt)
                
                # Skip variants where both haplotypes are reference (0|0)
                if diploid_gt == '0|0':
                    continue
                    
                # Build new line
                new_parts = parts[:9] + [diploid_gt]
                f_out.write('\t'.join(new_parts) + '\n')

    print(f"Conversion completed for sample: {sample_name}", file=sys.stderr)

if __name__ == "__main__":
    main()
