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
from dataclasses import asdict

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


def load_reference_genome(fasta_path: Path) -> Dict[str, SeqRecord]:
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
        # Note: This is a simplified application. For true VCF representation,
        # the insertion happens *between* pos-1 and pos.
        # Here, we are replacing the ref_allele at pos with ref_allele + insert_seq.
        # This might not be 100% accurate for all VCF interpretations but works for FASTA generation.
        if pos < len(seq_list):
            seq_list[pos] = variant.ref_allele + variant.insert_sequence
    
    elif isinstance(variant, Deletion):
        # Replace deleted region with first base
        end_pos = pos + variant.deleted_length
        if pos < len(seq_list) and end_pos <= len(seq_list):
            seq_list[pos:end_pos] = [variant.alt_allele]
    
    # Add other variant types if needed (Inversion, Duplication, ComplexVariant, Translocation)
    # For now, only SNV, Insertion, Deletion are handled by generate_* functions.
    
    return ''.join(seq_list)


def apply_multiple_variants_to_sequence(sequence: str, variants: List[Variant], chromosome: str) -> str:
    """Apply multiple variants to a sequence, handling position shifts from indels."""
    # Filter variants for this chromosome and sort by position (descending to avoid position shifts)
    chr_variants = [v for v in variants if v.chromosome == chromosome]
    chr_variants.sort(key=lambda v: v.position, reverse=True)
    
    modified_seq = sequence
    for variant in chr_variants:
        modified_seq = apply_variant_to_sequence(modified_seq, variant)
    
    return modified_seq


def create_modified_genome_with_variants(genome: Dict[str, SeqRecord], variants: List[Variant]) -> Dict[str, SeqRecord]:
    """Create a complete modified genome with multiple variants applied."""
    modified_genome = {}
    
    # Group variants by chromosome
    variants_by_chr = {}
    for variant in variants:
        if variant.chromosome not in variants_by_chr:
            variants_by_chr[variant.chromosome] = []
        variants_by_chr[variant.chromosome].append(variant)
    
    for chrom_id, record in genome.items():
        if chrom_id in variants_by_chr:
            # Apply all variants to this chromosome
            original_seq = str(record.seq)
            modified_seq = apply_multiple_variants_to_sequence(original_seq, variants_by_chr[chrom_id], chrom_id)
            
            # Create description with all variants
            variant_descriptions = []
            for variant in sorted(variants_by_chr[chrom_id], key=lambda v: v.position):
                variant_descriptions.append(f"{variant.variant_type} at {variant.position}")
            
            description = f"{record.description} [Modified with: {'; '.join(variant_descriptions)}]"
            
            # Create modified record
            modified_record = SeqRecord(
                Seq(modified_seq),
                id=record.id,
                description=description
            )
            modified_genome[chrom_id] = modified_record
        else:
            # Keep original chromosome unchanged
            modified_genome[chrom_id] = record
    
    return modified_genome


def generate_variant_set(genome: Dict[str, SeqRecord], variants_per_haplotype: int, 
                        variant_types: List[str]) -> List[Variant]:
    """Generate a set of variants for a single haplotype."""
    variants = []
    
    # Calculate number of each variant type
    type_counts = {}
    base_count = variants_per_haplotype // len(variant_types)
    remainder = variants_per_haplotype % len(variant_types)
    
    for i, vtype in enumerate(variant_types):
        type_counts[vtype] = base_count + (1 if i < remainder else 0)
    
    logging.info(f"Generating variant set with: {type_counts}")
    
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
                # Add other variant types here if their generate_* functions are implemented
                # elif vtype == 'inversion':
                #     variant = generate_inversion(genome, variants)
                # elif vtype == 'duplication':
                #     variant = generate_duplication(genome, variants)
                # elif vtype == 'complex':
                #     variant = generate_complex_variant(genome, variants)
                else:
                    logging.warning(f"Unknown or unimplemented variant type: {vtype}. Skipping.")
                    continue
                
                variants.append(variant)
                logging.info(f"Generated {vtype} {i+1}/{count}: {variant.to_english_description()}")
                
            except RuntimeError as e:
                logging.error(f"Failed to generate {vtype}: {e}")
    
    return variants


def generate_all_variant_sets_for_samples(genome: Dict[str, SeqRecord], num_samples: int, 
                                         variants_per_haplotype: int, variant_types: List[str],
                                         output_dir: Path) -> None:
    """
    Generates variant sets for multiple diploid samples and writes their FASTA files.
    Each sample will have two haplotypes, each with its own set of variants.
    """
    logging.info(f"Generating variants for {num_samples} diploid samples.")

    for i in range(num_samples):
        sample_name = f"sample_{i:03d}"
        sample_output_dir = output_dir / sample_name
        sample_output_dir.mkdir(parents=True, exist_ok=True)

        logging.info(f"Generating variants for {sample_name} (Haplotype 1)...")
        haplotype1_variants = generate_variant_set(genome, variants_per_haplotype, variant_types)
        
        logging.info(f"Generating variants for {sample_name} (Haplotype 2)...")
        haplotype2_variants = generate_variant_set(genome, variants_per_haplotype, variant_types)
        
        # Write Haplotype 1 FASTA
        modified_genome_hap1 = create_modified_genome_with_variants(genome, haplotype1_variants)
        hap1_fasta_path = sample_output_dir / f"{sample_name}_hap1.fasta"
        with open(hap1_fasta_path, "w") as f:
            # Write all chromosomes in the same order as original reference
            for chrom_id in genome.keys():
                SeqIO.write(modified_genome_hap1[chrom_id], f, 'fasta')
        logging.info(f"Modified genome for {sample_name}_hap1 written to {hap1_fasta_path}")

        # Write Haplotype 1 JSON (for debugging/record-keeping)
        hap1_json_path = sample_output_dir / f"{sample_name}_hap1_variants.json"
        with open(hap1_json_path, "w") as f:
            json.dump([asdict(v) for v in haplotype1_variants], f, indent=2)
        logging.info(f"Variant data for {sample_name}_hap1 written to {hap1_json_path}")

        # Write Haplotype 2 FASTA
        modified_genome_hap2 = create_modified_genome_with_variants(genome, haplotype2_variants)
        hap2_fasta_path = sample_output_dir / f"{sample_name}_hap2.fasta"
        with open(hap2_fasta_path, "w") as f:
            # Write all chromosomes in the same order as original reference
            for chrom_id in genome.keys():
                SeqIO.write(modified_genome_hap2[chrom_id], f, 'fasta')
        logging.info(f"Modified genome for {sample_name}_hap2 written to {hap2_fasta_path}")

        # Write Haplotype 2 JSON (for debugging/record-keeping)
        hap2_json_path = sample_output_dir / f"{sample_name}_hap2_variants.json"
        with open(hap2_json_path, "w") as f:
            json.dump([asdict(v) for v in haplotype2_variants], f, indent=2)
        logging.info(f"Variant data for {sample_name}_hap2 written to {hap2_json_path}")

    logging.info(f"All variant data for {num_samples} samples generated.")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Generate random genetic variants for diploid samples')
    parser.add_argument('--reference', type=Path, required=True, help='Input reference FASTA file')
    parser.add_argument('--num-samples', type=int, default=1, help='Number of diploid samples to generate')
    parser.add_argument('--variants-per-haplotype', type=int, default=10, help='Number of variants per haplotype')
    parser.add_argument('--output-dir', type=Path, default='data/test/variants/', help='Output directory')
    parser.add_argument('--variant-types', default='snv,insertion,deletion', 
                       help='Comma-separated list of variant types to generate (snv,insertion,deletion)')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    
    args = parser.parse_args()
    
    setup_logging()
    
    # Set random seed
    random.seed(args.seed)
    
    # Parse variant types
    variant_types = [vtype.strip().lower() for vtype in args.variant_types.split(',')]
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load reference genome
    genome = load_reference_genome(args.reference)
    
    # Generate all variant sets for samples
    generate_all_variant_sets_for_samples(
        genome, args.num_samples, args.variants_per_haplotype, variant_types, args.output_dir
    )
    
    logging.info(f"All outputs written to {args.output_dir}")


if __name__ == '__main__':
    main()
