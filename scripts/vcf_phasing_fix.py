import sys

def create_phased_vcf(input_file, output_file, sample_name, hap1_name, hap2_name):
    """
    Reads an input VCF, extracts specific haplotype columns, and writes a new
    VCF with a single phased genotype column for the given sample.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        hap1_col = -1
        hap2_col = -1
        ref_col = -1
        
        for line in infile:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    # Parse header to find column indices
                    parts = line.strip().split('\t')
                    for i, part in enumerate(parts[9:], 9):  # Sample columns start at index 9
                        if part == 'reference':
                            ref_col = i
                        elif part == hap1_name:
                            hap1_col = i
                        elif part == hap2_name:
                            hap2_col = i
                    
                    # Create new header with just our sample
                    new_header = parts[:9] + [sample_name]
                    outfile.write('\t'.join(new_header) + '\n')
                else:
                    outfile.write(line)
            else:
                # Process variant lines
                parts = line.strip().split('\t')
                if len(parts) >= 10:
                    # Extract genotypes for reference and haplotypes
                    # Default to "0" if column not found or empty
                    ref_gt = parts[ref_col] if ref_col != -1 and ref_col < len(parts) else "0"
                    hap1_gt = parts[hap1_col] if hap1_col != -1 and hap1_col < len(parts) else "0"
                    hap2_gt = parts[hap2_col] if hap2_col != -1 and hap2_col < len(parts) else "0"
                    
                    # Clean up genotype strings (remove . or empty)
                    if hap1_gt == "." or hap1_gt == "":
                        hap1_gt = "0"
                    if hap2_gt == "." or hap2_gt == "":
                        hap2_gt = "0"
                    
                    # Create phased genotype
                    phased_gt = hap1_gt + "|" + hap2_gt
                    
                    # Create output line with phased genotype
                    new_parts = parts[:9] + [phased_gt]
                    outfile.write('\t'.join(new_parts) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python vcf_phasing_fix.py <input_vcf> <output_vcf> <sample_name> <hap1_path_name> <hap2_path_name>", file=sys.stderr)
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    sample_name = sys.argv[3]
    hap1_path_name = sys.argv[4]
    hap2_path_name = sys.argv[5]
    
    create_phased_vcf(input_vcf, output_vcf, sample_name, hap1_path_name, hap2_path_name)

