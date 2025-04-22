#!/usr/bin/env python3
"""
Generate test data for feature analysis with various sequence modifications.

This script creates:
1. A reference FASTA file with sample sequences
2. A GFF3 file with example features on the reference
3. A GFA file with paths containing SNPs, indels, and structural variants
"""

import os
import random
import argparse
from typing import List, Dict, Tuple

# Ensure data directory exists
data_dir = os.path.dirname(os.path.abspath(__file__))
os.makedirs(data_dir, exist_ok=True)

def generate_random_sequence(length: int) -> str:
    """Generate a random DNA sequence of specified length."""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))

def create_reference_fasta(filename: str, seq_length: int = 10000) -> str:
    """Create a reference FASTA file with one sequence."""
    ref_seq = generate_random_sequence(seq_length)
    
    with open(filename, 'w') as f:
        f.write(f">reference\n")
        
        # Write in 80-character lines
        for i in range(0, len(ref_seq), 80):
            f.write(f"{ref_seq[i:i+80]}\n")
    
    return ref_seq

def create_gff_features(ref_length: int) -> List[Dict]:
    """Create a list of GFF features."""
    features = []
    
    # Create 3 genes of different sizes
    gene_coords = [
        (500, 2500),    # Gene 1
        (4000, 5500),   # Gene 2
        (7000, 9000)    # Gene 3
    ]
    
    for i, (start, end) in enumerate(gene_coords):
        gene_id = f"gene_{i+1}"
        gene_length = end - start
        strand = '+' if random.random() > 0.3 else '-'
        
        # Add gene feature
        features.append({
            'type': 'gene',
            'start': start,
            'end': end,
            'strand': strand,
            'id': gene_id,
            'attributes': {
                'ID': gene_id,
                'Name': f"GENE{i+1}"
            }
        })
        
        # Add mRNA feature
        mrna_id = f"mRNA_{i+1}"
        features.append({
            'type': 'mRNA',
            'start': start,
            'end': end,
            'strand': strand,
            'id': mrna_id,
            'attributes': {
                'ID': mrna_id,
                'Parent': gene_id,
                'Name': f"TRANSCRIPT{i+1}"
            }
        })
        
        # Add 3-7 exons per gene
        num_exons = random.randint(3, 7)
        exon_starts = sorted(random.sample(range(start, end-50), num_exons))
        
        for j, exon_start in enumerate(exon_starts):
            # Exon length between 50-300bp
            exon_length = random.randint(50, min(300, end - exon_start))
            exon_end = exon_start + exon_length
            
            exon_id = f"exon_{i+1}_{j+1}"
            features.append({
                'type': 'exon',
                'start': exon_start,
                'end': exon_end,
                'strand': strand,
                'id': exon_id,
                'attributes': {
                    'ID': exon_id,
                    'Parent': mrna_id
                }
            })
    
    return features

def write_gff3_file(filename: str, features: List[Dict]) -> None:
    """Write features to a GFF3 file."""
    with open(filename, 'w') as f:
        f.write("##gff-version 3\n")
        
        for feature in features:
            attr_str = ';'.join([f"{k}={v}" for k, v in feature['attributes'].items()])
            f.write(f"reference\t.\t{feature['type']}\t{feature['start']}\t{feature['end']}\t.\t"
                    f"{feature['strand']}\t.\t{attr_str}\n")

def create_variant_sequences(ref_seq: str, features: List[Dict]) -> Dict[str, str]:
    """Create variant sequences based on the reference."""
    variant_seqs = {}
    
    # Path 1: SNPs only (10 random SNPs)
    path1_seq = list(ref_seq)
    snp_positions = random.sample(range(len(ref_seq)), 10)
    for pos in snp_positions:
        original = path1_seq[pos]
        bases = [b for b in ['A', 'C', 'G', 'T'] if b != original]
        path1_seq[pos] = random.choice(bases)
    variant_seqs["path1_snps"] = ''.join(path1_seq)
    
    # Path 2: Small indels (3 small insertions, 3 small deletions)
    path2_seq = list(ref_seq)
    # Insertions
    for _ in range(3):
        pos = random.randint(100, len(ref_seq) - 100)
        insertion = generate_random_sequence(random.randint(5, 15))
        path2_seq[pos:pos] = insertion
    
    # Convert to string for deletions
    path2_seq = ''.join(path2_seq)
    # Deletions
    for _ in range(3):
        pos = random.randint(100, len(path2_seq) - 100)
        del_len = random.randint(5, 15)
        path2_seq = path2_seq[:pos] + path2_seq[pos+del_len:]
    
    variant_seqs["path2_indels"] = path2_seq
    
    # Path 3: Feature-affecting variants (targeted at gene, exon features)
    path3_seq = ref_seq
    
    # Find gene and exon features
    genes = [f for f in features if f['type'] == 'gene']
    exons = [f for f in features if f['type'] == 'exon']
    
    # Introduce deletion in an exon (first exon of gene 1)
    gene1_exons = [e for e in exons if e['attributes']['Parent'] == "mRNA_1"]
    if gene1_exons:
        exon = gene1_exons[0]
        exon_start, exon_end = exon['start'], exon['end']
        # Delete 30% of the exon
        del_start = exon_start + int((exon_end - exon_start) * 0.3)
        del_end = del_start + int((exon_end - exon_start) * 0.3)
        path3_seq = path3_seq[:del_start] + path3_seq[del_end:]
    
    # Introduce SNPs in gene 2
    gene2 = next((g for g in genes if g['id'] == 'gene_2'), None)
    if gene2:
        gene_start, gene_end = gene2['start'], gene2['end']
        # Add 20 SNPs in gene2
        path3_seq_list = list(path3_seq)
        snp_positions = random.sample(range(gene_start, gene_end), 20)
        for pos in snp_positions:
            if pos < len(path3_seq_list):
                original = path3_seq_list[pos]
                bases = [b for b in ['A', 'C', 'G', 'T'] if b != original]
                path3_seq_list[pos] = random.choice(bases)
        path3_seq = ''.join(path3_seq_list)
    
    # Introduce large structural change in gene 3
    gene3 = next((g for g in genes if g['id'] == 'gene_3'), None)
    if gene3:
        gene_start, gene_end = gene3['start'], gene3['end']
        # Invert a large segment
        seg_len = (gene_end - gene_start) // 2
        seg_start = gene_start + (gene_end - gene_start) // 4
        segment = path3_seq[seg_start:seg_start+seg_len]
        inverted = segment[::-1]  # Reverse the segment
        path3_seq = path3_seq[:seg_start] + inverted + path3_seq[seg_start+seg_len:]
    
    variant_seqs["path3_feature_mods"] = path3_seq
    
    return variant_seqs

def create_gfa_file(filename: str, ref_seq: str, variant_seqs: Dict[str, str]) -> None:
    """Create a GFA file with the reference and variant paths."""
    with open(filename, 'w') as f:
        # Header
        f.write("H\tVN:Z:1.0\n")
        
        # Segments (just use whole sequences as single segments for simplicity)
        f.write(f"S\tref\t{ref_seq}\n")
        
        for path_id, seq in variant_seqs.items():
            f.write(f"S\t{path_id}\t{seq}\n")
        
        # Paths
        f.write(f"P\treference\tref+\n")
        for path_id in variant_seqs:
            f.write(f"P\t{path_id}\t{path_id}+\n")

def main():
    """Generate test data files."""
    parser = argparse.ArgumentParser(description="Generate test data for feature analysis.")
    parser.add_argument('--output-dir', default=data_dir, help='Output directory for test files')
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("Generating test data with feature variants...")
    
    # Create reference sequence
    ref_fasta = os.path.join(args.output_dir, "test_reference.fasta")
    ref_seq = create_reference_fasta(ref_fasta)
    print(f"Created reference FASTA: {ref_fasta}")
    
    # Create annotation features
    features = create_gff_features(len(ref_seq))
    
    # Write GFF3 file
    gff_file = os.path.join(args.output_dir, "test_features.gff3")
    write_gff3_file(gff_file, features)
    print(f"Created GFF3 file: {gff_file}")
    
    # Create variant sequences
    variant_seqs = create_variant_sequences(ref_seq, features)
    
    # Create GFA file
    gfa_file = os.path.join(args.output_dir, "test_variants.gfa")
    create_gfa_file(gfa_file, ref_seq, variant_seqs)
    print(f"Created GFA file: {gfa_file}")
    
    print(f"All test data generated in {args.output_dir}")
    print("\nTest data includes:")
    print("- Path 1 (path1_snps): 10 random SNPs")
    print("- Path 2 (path2_indels): 3 small insertions and 3 small deletions")
    print("- Path 3 (path3_feature_mods): Exon deletion, gene with SNPs, and inverted gene segment")

if __name__ == "__main__":
    main()
