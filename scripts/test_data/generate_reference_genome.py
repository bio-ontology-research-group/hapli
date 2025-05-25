#!/usr/bin/env python3
"""
Generate a synthetic reference genome resembling hg38.

This script creates a synthetic reference genome with realistic characteristics
similar to the human genome (hg38). It generates chromosomes with:
- Weighted nucleotide selection to match target GC content
- Repetitive regions (5-10% of each chromosome)
- N regions representing unknown sequences (1-2% of length)
- CpG islands with elevated GC content
- Telomeric repeats at chromosome ends
- Centromeric regions with tandem repeats
- Realistic chromosome naming and sizing

Algorithm:
1. Calculate chromosome lengths using exponential decay
2. For each chromosome:
   - Generate base sequence with target GC content
   - Add CpG islands near chromosome start
   - Add telomeric repeats at chromosome ends
   - Add centromeric region in the middle
   - Add repetitive regions (simple tandem repeats)
   - Insert N regions randomly
   - Write to FASTA with proper headers

Parameters can be adjusted to create genomes of different sizes and characteristics
for testing purposes.
"""

import argparse
import logging
import random
import sys
from pathlib import Path
from typing import List, Tuple

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)


def setup_logging() -> None:
    """Configure logging for progress reporting."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def generate_chromosome_names(num_chromosomes: int) -> List[str]:
    """
    Generate realistic chromosome names.
    
    Args:
        num_chromosomes: Number of chromosomes to generate
        
    Returns:
        List of chromosome names (chr1, chr2, ..., chrX, chrY, chrM)
    """
    names = []
    
    # Regular autosomes
    regular_count = min(num_chromosomes, 22)
    for i in range(1, regular_count + 1):
        names.append(f"chr{i}")
    
    # Special chromosomes if we need more
    remaining = num_chromosomes - regular_count
    if remaining > 0:
        special = ["chrX", "chrY", "chrM"]
        for i in range(min(remaining, len(special))):
            names.append(special[i])
        remaining -= min(remaining, len(special))
    
    # Additional chromosomes if still needed
    for i in range(remaining):
        names.append(f"chr{22 + len(special) + i + 1}")
    
    return names


def calculate_chromosome_lengths(num_chromosomes: int, base_length: int, decay: float) -> List[int]:
    """
    Calculate chromosome lengths using exponential decay.
    
    Args:
        num_chromosomes: Number of chromosomes
        base_length: Length of first chromosome
        decay: Decay factor for subsequent chromosomes
        
    Returns:
        List of chromosome lengths
    """
    lengths = []
    current_length = base_length
    
    for i in range(num_chromosomes):
        lengths.append(int(current_length))
        current_length *= decay
    
    return lengths


def generate_weighted_sequence(length: int, gc_content: float) -> str:
    """
    Generate a random DNA sequence with target GC content.
    
    Args:
        length: Length of sequence to generate
        gc_content: Target GC content (0.0 to 1.0)
        
    Returns:
        Random DNA sequence string
    """
    # Calculate nucleotide probabilities
    gc_prob = gc_content / 2  # Split GC content between G and C
    at_prob = (1 - gc_content) / 2  # Split AT content between A and T
    
    nucleotides = ['A', 'T', 'G', 'C']
    weights = [at_prob, at_prob, gc_prob, gc_prob]
    
    sequence = random.choices(nucleotides, weights=weights, k=length)
    return ''.join(sequence)


def add_cpg_islands(sequence: str, chr_length: int) -> str:
    """
    Add CpG islands with elevated GC content near chromosome start.
    
    Args:
        sequence: Input DNA sequence
        chr_length: Total chromosome length
        
    Returns:
        Sequence with CpG islands added
    """
    seq_list = list(sequence)
    
    # Add 2-5 CpG islands per chromosome
    num_islands = random.randint(2, 5)
    
    # Place within first 10% of chromosome
    search_region = int(chr_length * 0.1)
    
    for _ in range(num_islands):
        # Island size: 500-2000 bp
        island_length = random.randint(500, 2000)
        
        # Make sure we don't exceed sequence bounds
        if search_region < island_length:
            continue
            
        # Random position within search region
        start_pos = random.randint(0, search_region - island_length)
        
        # Generate high GC content sequence (70-80%)
        high_gc_content = random.uniform(0.70, 0.80)
        island_seq = generate_weighted_sequence(island_length, high_gc_content)
        
        # Replace the region
        for i, nucleotide in enumerate(island_seq):
            if start_pos + i < len(seq_list):
                seq_list[start_pos + i] = nucleotide
    
    return ''.join(seq_list)


def add_telomeric_repeats(sequence: str) -> str:
    """
    Add telomeric repeats (TTAGGG) at chromosome ends.
    
    Args:
        sequence: Input DNA sequence
        
    Returns:
        Sequence with telomeric repeats added
    """
    seq_list = list(sequence)
    seq_length = len(sequence)
    
    # Number of repeats at each end: 50-200
    num_repeats_start = random.randint(50, 200)
    num_repeats_end = random.randint(50, 200)
    
    # Telomeric repeat motif
    telomere_motif = "TTAGGG"
    
    # Add to start of chromosome
    start_length = min(num_repeats_start * len(telomere_motif), seq_length // 4)
    start_seq = (telomere_motif * (start_length // len(telomere_motif) + 1))[:start_length]
    
    # Add some degenerate copies (5% chance of mutation per base)
    start_seq_list = list(start_seq)
    for i in range(len(start_seq_list)):
        if random.random() < 0.05:  # 5% mutation rate
            start_seq_list[i] = random.choice(['A', 'T', 'G', 'C'])
    start_seq = ''.join(start_seq_list)
    
    # Replace start region
    for i, nucleotide in enumerate(start_seq):
        if i < len(seq_list):
            seq_list[i] = nucleotide
    
    # Add to end of chromosome
    end_length = min(num_repeats_end * len(telomere_motif), seq_length // 4)
    end_seq = (telomere_motif * (end_length // len(telomere_motif) + 1))[:end_length]
    
    # Add some degenerate copies
    end_seq_list = list(end_seq)
    for i in range(len(end_seq_list)):
        if random.random() < 0.05:  # 5% mutation rate
            end_seq_list[i] = random.choice(['A', 'T', 'G', 'C'])
    end_seq = ''.join(end_seq_list)
    
    # Replace end region
    start_pos = seq_length - end_length
    for i, nucleotide in enumerate(end_seq):
        if start_pos + i < len(seq_list):
            seq_list[start_pos + i] = nucleotide
    
    return ''.join(seq_list)


def add_centromeric_region(sequence: str, chr_length: int) -> str:
    """
    Add centromeric region with tandem repeats and Ns.
    
    Args:
        sequence: Input DNA sequence
        chr_length: Total chromosome length
        
    Returns:
        Sequence with centromeric region added
    """
    seq_list = list(sequence)
    
    # Centromere size: 1-5% of chromosome length
    centromere_fraction = random.uniform(0.01, 0.05)
    centromere_length = int(chr_length * centromere_fraction)
    
    # Position: 40-60% of chromosome length
    center_position = random.uniform(0.40, 0.60)
    start_pos = int(chr_length * center_position - centromere_length // 2)
    start_pos = max(0, min(start_pos, chr_length - centromere_length))
    
    # Centromeric repeat motifs (alpha satellite-like)
    centromere_motifs = ["ATATACATAG", "CATTCCATTC", "AATCAACCC", "TTCCATTCCATTC"]
    
    # Fill centromere with repeats and Ns
    centromere_seq = []
    pos = 0
    while pos < centromere_length:
        if random.random() < 0.3:  # 30% chance of N region
            n_length = min(random.randint(5, 50), centromere_length - pos)
            centromere_seq.extend(['N'] * n_length)
            pos += n_length
        else:
            # Add tandem repeat
            motif = random.choice(centromere_motifs)
            repeat_length = min(random.randint(20, 200), centromere_length - pos)
            repeat_seq = (motif * (repeat_length // len(motif) + 1))[:repeat_length]
            centromere_seq.extend(list(repeat_seq))
            pos += repeat_length
    
    # Replace the centromeric region
    for i, nucleotide in enumerate(centromere_seq):
        if start_pos + i < len(seq_list):
            seq_list[start_pos + i] = nucleotide
    
    return ''.join(seq_list)


def add_repetitive_regions(sequence: str, repeat_fraction: float = 0.075) -> str:
    """
    Add repetitive regions to simulate tandem repeats.
    
    Args:
        sequence: Input DNA sequence
        repeat_fraction: Fraction of sequence to make repetitive
        
    Returns:
        Sequence with repetitive regions added
    """
    seq_list = list(sequence)
    seq_length = len(sequence)
    repeat_length = int(seq_length * repeat_fraction)
    
    # Common repeat motifs
    repeat_motifs = ["CA", "GT", "ATCG", "TAGA", "CACA", "GTGT"]
    
    # Add several repetitive regions
    num_regions = random.randint(3, 8)
    for _ in range(num_regions):
        if repeat_length <= 0:
            break
            
        # Choose random position and length for this repeat
        region_length = min(random.randint(50, 500), repeat_length)
        start_pos = random.randint(0, seq_length - region_length)
        
        # Choose a repeat motif and create the repetitive sequence
        motif = random.choice(repeat_motifs)
        repeat_seq = (motif * (region_length // len(motif) + 1))[:region_length]
        
        # Replace the region
        for i, nucleotide in enumerate(repeat_seq):
            if start_pos + i < seq_length:
                seq_list[start_pos + i] = nucleotide
        
        repeat_length -= region_length
    
    return ''.join(seq_list)


def add_n_regions(sequence: str, n_fraction: float = 0.015) -> str:
    """
    Add N regions to simulate unknown sequences.
    
    Args:
        sequence: Input DNA sequence
        n_fraction: Fraction of sequence to convert to Ns
        
    Returns:
        Sequence with N regions added
    """
    seq_list = list(sequence)
    seq_length = len(sequence)
    n_length = int(seq_length * n_fraction)
    
    # Add several N regions
    num_regions = random.randint(2, 6)
    for _ in range(num_regions):
        if n_length <= 0:
            break
            
        # Choose random position and length for this N region
        region_length = min(random.randint(10, 200), n_length)
        start_pos = random.randint(0, seq_length - region_length)
        
        # Replace with Ns
        for i in range(region_length):
            if start_pos + i < seq_length:
                seq_list[start_pos + i] = 'N'
        
        n_length -= region_length
    
    return ''.join(seq_list)


def generate_chromosome(name: str, length: int, gc_content: float, 
                       include_cpg_islands: bool = True, include_telomeres: bool = True, 
                       include_centromeres: bool = True) -> SeqRecord:
    """
    Generate a single chromosome with realistic features.
    
    Args:
        name: Chromosome name
        length: Chromosome length in bp
        gc_content: Target GC content
        include_cpg_islands: Whether to add CpG islands
        include_telomeres: Whether to add telomeric repeats
        include_centromeres: Whether to add centromeric regions
        
    Returns:
        SeqRecord object representing the chromosome
    """
    logging.info(f"Generating {name} ({length:,} bp)")
    
    # Generate base sequence
    sequence = generate_weighted_sequence(length, gc_content)
    
    # Add realistic genomic features
    if include_cpg_islands:
        sequence = add_cpg_islands(sequence, length)
    
    if include_telomeres:
        sequence = add_telomeric_repeats(sequence)
    
    if include_centromeres:
        sequence = add_centromeric_region(sequence, length)
    
    # Add repetitive regions
    sequence = add_repetitive_regions(sequence)
    
    # Add N regions
    sequence = add_n_regions(sequence)
    
    # Create chromosome number for header
    chr_num = name.replace('chr', '')
    if chr_num.isdigit():
        chr_description = f"Homo sapiens chromosome {chr_num}, Test Assembly"
    elif chr_num in ['X', 'Y', 'M']:
        chr_description = f"Homo sapiens chromosome {chr_num}, Test Assembly"
    else:
        chr_description = f"Homo sapiens chromosome {chr_num}, Test Assembly"
    
    # Create SeqRecord with hg38-like metadata
    seq_record = SeqRecord(
        Seq(sequence),
        id=name,
        description=f"{chr_description} (length={len(sequence)})"
    )
    
    return seq_record


def main():
    """Main function to generate synthetic reference genome."""
    parser = argparse.ArgumentParser(
        description="Generate a synthetic reference genome resembling hg38",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--num-chromosomes",
        type=int,
        default=3,
        help="Number of chromosomes to generate"
    )
    
    parser.add_argument(
        "--chr-length",
        type=int,
        default=100000,
        help="Base length for chromosome 1 in bp"
    )
    
    parser.add_argument(
        "--length-decay",
        type=float,
        default=0.9,
        help="Each chromosome is this fraction of the previous"
    )
    
    parser.add_argument(
        "--gc-content",
        type=float,
        default=0.41,
        help="Target GC content (hg38 is ~41%%)"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        default="data/test/reference.fa",
        help="Output FASTA file"
    )
    
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility"
    )
    
    parser.add_argument(
        "--add-cpg-islands",
        action="store_true",
        default=True,
        help="Add CpG islands with elevated GC content"
    )
    
    parser.add_argument(
        "--no-cpg-islands",
        action="store_false",
        dest="add_cpg_islands",
        help="Do not add CpG islands"
    )
    
    parser.add_argument(
        "--add-telomeres",
        action="store_true",
        default=True,
        help="Add telomeric repeats at chromosome ends"
    )
    
    parser.add_argument(
        "--no-telomeres",
        action="store_false",
        dest="add_telomeres",
        help="Do not add telomeric repeats"
    )
    
    parser.add_argument(
        "--add-centromeres",
        action="store_true",
        default=True,
        help="Add centromeric regions with tandem repeats"
    )
    
    parser.add_argument(
        "--no-centromeres",
        action="store_false",
        dest="add_centromeres",
        help="Do not add centromeric regions"
    )
    
    args = parser.parse_args()
    
    # Setup
    setup_logging()
    random.seed(args.seed)
    
    # Validate arguments
    if args.num_chromosomes < 1:
        logging.error("Number of chromosomes must be at least 1")
        sys.exit(1)
    
    if not 0.0 <= args.gc_content <= 1.0:
        logging.error("GC content must be between 0.0 and 1.0")
        sys.exit(1)
    
    if args.chr_length < 1000:
        logging.error("Chromosome length must be at least 1000 bp")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Generate chromosome specifications
    chr_names = generate_chromosome_names(args.num_chromosomes)
    chr_lengths = calculate_chromosome_lengths(args.num_chromosomes, args.chr_length, args.length_decay)
    
    logging.info(f"Generating {args.num_chromosomes} chromosomes")
    logging.info(f"Total genome size: {sum(chr_lengths):,} bp")
    logging.info(f"Target GC content: {args.gc_content:.1%}")
    logging.info(f"Features: CpG islands={args.add_cpg_islands}, Telomeres={args.add_telomeres}, Centromeres={args.add_centromeres}")
    logging.info(f"Output file: {args.output}")
    
    # Generate chromosomes
    chromosomes = []
    for name, length in zip(chr_names, chr_lengths):
        chromosome = generate_chromosome(
            name, length, args.gc_content,
            args.add_cpg_islands, args.add_telomeres, args.add_centromeres
        )
        chromosomes.append(chromosome)
    
    # Write to FASTA file
    logging.info(f"Writing FASTA file: {args.output}")
    with open(args.output, 'w') as output_file:
        SeqIO.write(chromosomes, output_file, "fasta")
    
    logging.info("Synthetic reference genome generation complete!")
    
    # Report statistics
    total_length = sum(len(chr.seq) for chr in chromosomes)
    total_gc = sum(chr.seq.count('G') + chr.seq.count('C') for chr in chromosomes)
    total_n = sum(chr.seq.count('N') for chr in chromosomes)
    actual_gc_content = total_gc / (total_length - total_n) if (total_length - total_n) > 0 else 0
    
    logging.info(f"Final statistics:")
    logging.info(f"  Total length: {total_length:,} bp")
    logging.info(f"  Actual GC content: {actual_gc_content:.1%}")
    logging.info(f"  N content: {total_n / total_length:.1%}")
    logging.info(f"  Chromosomes: {len(chromosomes)}")


if __name__ == "__main__":
    main()
