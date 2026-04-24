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

def get_chromosome_lengths(reference_file):
    """Get chromosome lengths from reference FASTA"""
    chrom_lengths = {}
    current_chrom = None
    current_length = 0
    
    try:
        if not os.path.exists(reference_file):
            alt_path = reference_file.replace('/data/', '')
            if os.path.exists(alt_path):
                reference_file = alt_path
            else:
                return {}
                
        with open(reference_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous chromosome
                    if current_chrom:
                        chrom_lengths[current_chrom] = current_length
                    # Start new chromosome
                    current_chrom = line[1:].split()[0]
                    current_length = 0
                else:
                    current_length += len(line)
            
            # Save last chromosome
            if current_chrom:
                chrom_lengths[current_chrom] = current_length
                
    except Exception as e:
        print(f"Error reading reference file {reference_file}: {e}", file=sys.stderr)
    
    return chrom_lengths

def map_position_to_chromosome(position, chrom_lengths):
    """Map a position on concatenated reference to actual chromosome"""
    cumulative_pos = 0
    for chrom, length in chrom_lengths.items():
        if position <= cumulative_pos + length:
            # Position is in this chromosome
            chrom_pos = position - cumulative_pos
            return chrom, chrom_pos
        cumulative_pos += length
    
    # If we get here, position is beyond all chromosomes
    # Return the last chromosome
    if chrom_lengths:
        last_chrom = list(chrom_lengths.keys())[-1]
        return last_chrom, position - (cumulative_pos - chrom_lengths[last_chrom])
    
    return "chr1", position

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
    reference_file = None
    
    with open(input_file, 'r') as f_in:
        for line in f_in:
            if line.startswith('##reference=file://'):
                # Extract reference file path
                ref_path = line.strip().split('file://', 1)[1]
                # Remove any /data prefix that might be from Docker mount
                if ref_path.startswith('/data/'):
                    ref_path = ref_path[6:]
                reference_file = ref_path
                reference_chroms = get_all_reference_chromosome_names(ref_path)
                print(f"Found reference file: {ref_path}", file=sys.stderr)
                print(f"Found reference chromosomes: {reference_chroms}", file=sys.stderr)
                break
            elif line.startswith('#CHROM'):
                break

    # Get chromosome lengths for position mapping
    chrom_lengths = {}
    if reference_file:
        chrom_lengths = get_chromosome_lengths(reference_file)
        print(f"Chromosome lengths: {chrom_lengths}", file=sys.stderr)

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
                        if chrom in chrom_lengths:
                            new_line = f"##contig=<ID={chrom},length={chrom_lengths[chrom]}>\n"
                        else:
                            new_line = f"##contig=<ID={chrom}>\n"
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
                
                # Map position to actual chromosome if we have chromosome lengths
                original_chrom = parts[0]
                original_pos = int(parts[1])
                
                if original_chrom == 'reference' and chrom_lengths:
                    # Map position to actual chromosome
                    actual_chrom, actual_pos = map_position_to_chromosome(original_pos, chrom_lengths)
                    parts[0] = actual_chrom
                    parts[1] = str(actual_pos)
                    print(f"Mapped reference:{original_pos} -> {actual_chrom}:{actual_pos}", file=sys.stderr)
                elif reference_chroms and parts[0] == 'reference':
                    # Fallback: just use the first chromosome
                    parts[0] = reference_chroms[0]
                    
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
