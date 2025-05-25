#!/usr/bin/env python3
"""
Generate random genetic variants from a reference genome.
"""

import argparse
import json
import logging
import os
import random
from pathlib import Path
from typing import Dict, List, Set, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from variant_models import SNV, Insertion, Deletion, Variant


def setup_logging() -> None:
    """Configure logging."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def load_reference_genome(fasta_path: str) -> Dict[str, SeqRecord]:
    """Load reference genome from FASTA file."""
    logging.info(f"Loading reference genome from {fasta_path}")
    
    genome = {}
    with open(fasta_path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            genome[record.id] = record
    
    logging.info(f"Loaded {len(genome)} chromosomes/contigs")
    return genome


def get_random_chromosome_position(genome: Dict[str, SeqRecord]) -> Tuple[str, int]:
    """Get a random chromosome and position from the genome."""
    # Weight by chromosome length
    chromosomes = list(genome.keys())
    weights = [len(genome[chrom].seq) for chrom in chromosomes]
    
    chromosome = random.choices(chromosomes, weights=weights)[0]
    position = random.randint(1, len(genome[chromosome].seq))
    
    return chromosome, position


def check_overlap(new_variant: Variant, existing_variants: List[Variant], buffer: int = 100) -> bool:
    """Check if a new variant overlaps with existing variants."""
    new_start, new_end = new_variant.get_affected_region()
    
    for existing in existing_variants:
        if existing.chromosome != new_variant.chromosome:
            continue
        
        existing_start, existing_end = existing.get_affected_region()
        
        # Check for overlap with buffer
        if (new_start - buffer <= existing_end + buffer and 
            new_end + buffer >= existing_start - buffer):
            return True
    
    return False


def generate_snv(genome: Dict[str, SeqRecord], existing_variants: List[Variant]) -> SNV:
    """Generate a random SNV."""
    max_attempts = 1000
    
    for _ in range(max_attempts):
        chromosome, position = get_random_chromosome_position(genome)
        
        # Get reference base
        ref_seq = str(genome[chromosome].seq)
        if position > len(ref_seq):
            continue
        
        ref_base = ref_seq[position - 1].upper()  # Convert to 0-based indexing
        
        # Skip if reference is N
        if ref_base == 'N':
            continue
        
        # Generate alternate base with transition/transversion ratio 2:1
        transitions = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
        transversions = {
            'A': ['C', 'T'], 'G': ['C', 'T'],
            'C': ['A', 'G'], 'T': ['A', 'G']
        }
        
        if random.random() < 2/3:  # Transition
            alt_base = transitions.get(ref_base, 'A')
        else:  # Transversion
            alt_base = random.choice(transversions.get(ref_base, ['A']))
        
        snv = SNV(
            chromosome=chromosome,
            position=position,
            ref_allele=ref_base,
            alt_allele=alt_base,
            variant_type="SNV"
        )
        
        if not check_overlap(snv, existing_variants):
            return snv
    
    raise RuntimeError("Could not generate non-overlapping SNV after maximum attempts")


def generate_insertion(genome: Dict[str, SeqRecord], existing_variants: List[Variant]) -> Insertion:
    """Generate a random insertion."""
    max_attempts = 1000
    
    for _ in range(max_attempts):
        chromosome, position = get_random_chromosome_position(genome)
        
        # Get reference base
        ref_seq = str(genome[chromosome].seq)
        if position > len(ref_seq):
            continue
        
        ref_base = ref_seq[position - 1].upper()
        
        # Skip if reference is N
        if ref_base == 'N':
            continue
        
        # Generate random insertion sequence (1-50 bp)
        insert_length = random.randint(1, 50)
        bases = ['A', 'T', 'G', 'C']
        insert_seq = ''.join(random.choices(bases, k=insert_length))
        
        insertion = Insertion(
            chromosome=chromosome,
            position=position,
            ref_allele=ref_base,
            alt_allele=ref_base + insert_seq,
            variant_type="Insertion",
            insert_sequence=insert_seq
        )
        
        if not check_overlap(insertion, existing_variants):
            return insertion
    
    raise RuntimeError("Could not generate non-overlapping insertion after maximum attempts")


def generate_deletion(genome: Dict[str, SeqRecord], existing_variants: List[Variant]) -> Deletion:
    """Generate a random deletion."""
    max_attempts = 1000
    
    for _ in range(max_attempts):
        chromosome, position = get_random_chromosome_position(genome)
        
        ref_seq = str(genome[chromosome].seq)
        
        # Generate deletion length (1-100 bp)
        max_del_length = min(100, len(ref_seq) - position + 1)
        if max_del_length < 1:
            continue
        
        del_length = random.randint(1, max_del_length)
        
        # Get reference sequence
        if position + del_length - 1 > len(ref_seq):
            continue
        
        ref_allele = ref_seq[position - 1:position + del_length - 1]
        
        # Skip if reference contains N
        if 'N' in ref_allele.upper():
            continue
        
        deletion = Deletion(
            chromosome=chromosome,
            position=position,
            ref_allele=ref_allele.upper(),
            alt_allele=ref_allele[0].upper(),  # Keep first base
            variant_type="Deletion",
            deleted_length=del_length
        )
        
        if not check_overlap(deletion, existing_variants):
            return deletion
    
    raise RuntimeError("Could not generate non-overlapping deletion after maximum attempts")


def apply_variant_to_sequence(sequence: str, variant: Variant) -> str:
    """Apply a variant to a sequence and return the modified sequence."""
    seq_list = list(sequence)
    pos = variant.position - 1  # Convert to 0-based
    
    if isinstance(variant, SNV):
        if pos < len(seq_list):
            seq_list[pos] = variant.alt_allele
    
    elif isinstance(variant, Insertion):
        # Insert after the reference position
        if pos < len(seq_list):
            seq_list[pos] = variant.ref_allele + variant.insert_sequence
    
    elif isinstance(variant, Deletion):
        # Replace deleted region with first base
        end_pos = pos + variant.deleted_length
        if pos < len(seq_list) and end_pos <= len(seq_list):
            seq_list[pos:end_pos] = [variant.alt_allele]
    
    return ''.join(seq_list)


def create_modified_genome(genome: Dict[str, SeqRecord], variant: Variant) -> Dict[str, SeqRecord]:
    """Create a complete modified genome with the variant applied."""
    modified_genome = {}
    
    for chrom_id, record in genome.items():
        if chrom_id == variant.chromosome:
            # Apply variant to this chromosome
            original_seq = str(record.seq)
            modified_seq = apply_variant_to_sequence(original_seq, variant)
            
            # Create modified record
            modified_record = SeqRecord(
                Seq(modified_seq),
                id=record.id,
                description=f"{record.description} [Modified with {variant.variant_type} at position {variant.position}]"
            )
            modified_genome[chrom_id] = modified_record
        else:
            # Keep original chromosome unchanged
            modified_genome[chrom_id] = record
    
    return modified_genome


def generate_variants(genome: Dict[str, SeqRecord], num_variants: int, 
                     variant_types: List[str]) -> List[Variant]:
    """Generate the specified number of variants."""
    variants = []
    
    # Calculate number of each variant type
    type_counts = {}
    base_count = num_variants // len(variant_types)
    remainder = num_variants % len(variant_types)
    
    for i, vtype in enumerate(variant_types):
        type_counts[vtype] = base_count + (1 if i < remainder else 0)
    
    logging.info(f"Generating variants: {type_counts}")
    
    # Generate variants
    for vtype, count in type_counts.items():
        for i in range(count):
            try:
                if vtype == 'snv':
                    variant = generate_snv(genome, variants)
                elif vtype == 'insertion':
                    variant = generate_insertion(genome, variants)
                elif vtype == 'deletion':
                    variant = generate_deletion(genome, variants)
                else:
                    logging.warning(f"Unknown variant type: {vtype}")
                    continue
                
                variants.append(variant)
                logging.info(f"Generated {vtype} {i+1}/{count}: {variant.to_english_description()}")
                
            except RuntimeError as e:
                logging.error(f"Failed to generate {vtype}: {e}")
    
    return variants


def write_variant_outputs(variant: Variant, genome: Dict[str, SeqRecord], output_dir: Path) -> None:
    """Write output files for a single variant."""
    variant_id = variant.variant_id[:8]  # Use first 8 chars of UUID
    
    # Write description
    desc_file = output_dir / f"{variant_id}_description.txt"
    with open(desc_file, 'w') as f:
        f.write(f"Variant ID: {variant.variant_id}\n")
        f.write(f"Type: {variant.variant_type}\n")
        f.write(f"Description: {variant.to_english_description()}\n")
        f.write(f"HGVS: {variant.to_hgvs_notation()}\n")
        f.write(f"Affected Region: {variant.get_affected_region()}\n")
    
    # Create complete modified genome
    modified_genome = create_modified_genome(genome, variant)
    
    # Write complete modified genome to FASTA
    fa_file = output_dir / f"{variant_id}_genome.fa"
    with open(fa_file, 'w') as f:
        # Write all chromosomes in the same order as original
        for chrom_id in genome.keys():
            SeqIO.write(modified_genome[chrom_id], f, 'fasta')
    
    logging.info(f"Written complete modified genome to {fa_file}")


def write_summary(variants: List[Variant], output_dir: Path) -> None:
    """Write summary JSON file."""
    summary_data = []
    
    for variant in variants:
        variant_data = {
            'variant_id': variant.variant_id,
            'variant_type': variant.variant_type,
            'chromosome': variant.chromosome,
            'position': variant.position,
            'ref_allele': variant.ref_allele,
            'alt_allele': variant.alt_allele,
            'description': variant.to_english_description(),
            'hgvs': variant.to_hgvs_notation(),
            'affected_region': variant.get_affected_region()
        }
        
        # Add type-specific fields
        if isinstance(variant, Insertion):
            variant_data['insert_sequence'] = variant.insert_sequence
        elif isinstance(variant, Deletion):
            variant_data['deleted_length'] = variant.deleted_length
        
        summary_data.append(variant_data)
    
    summary_file = output_dir / 'variants_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(summary_data, f, indent=2)
    
    logging.info(f"Summary written to {summary_file}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Generate random genetic variants')
    parser.add_argument('--reference', required=True, help='Input reference FASTA file')
    parser.add_argument('--num-variants', type=int, default=10, help='Total number of variants')
    parser.add_argument('--output-dir', default='data/test/variants/', help='Output directory')
    parser.add_argument('--variant-types', default='snv,insertion,deletion', 
                       help='Comma-separated list of variant types')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    
    args = parser.parse_args()
    
    setup_logging()
    
    # Set random seed
    random.seed(args.seed)
    
    # Parse variant types
    variant_types = [vtype.strip().lower() for vtype in args.variant_types.split(',')]
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load reference genome
    genome = load_reference_genome(args.reference)
    
    # Generate variants
    variants = generate_variants(genome, args.num_variants, variant_types)
    
    logging.info(f"Generated {len(variants)} variants")
    
    # Write outputs
    for variant in variants:
        write_variant_outputs(variant, genome, output_dir)
    
    write_summary(variants, output_dir)
    
    logging.info(f"All outputs written to {output_dir}")


if __name__ == '__main__':
    main()
