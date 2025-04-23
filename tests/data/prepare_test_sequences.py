#!/usr/bin/env python3
"""
Script to prepare test sequences for minimap2 wrapper testing.

This script:
1. Downloads a small viral genome (PhiX174) from NCBI
2. Extracts small regions for testing
3. Creates artificial variants for testing different alignment scenarios
4. Saves these as separate FASTA files in the tests/data directory
"""

import os
import random
import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import requests
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Set your email for Entrez
Entrez.email = "your_email@example.com"  # Replace with your email

# Constants
PHIX174_ACCESSION = "NC_001422.1"  # PhiX174 genome
OUTPUT_DIR = Path(__file__).parent
SEED = 42  # For reproducibility


def download_genome(accession: str, output_file: str) -> str:
    """
    Download a genome from NCBI using Entrez.
    
    Args:
        accession: NCBI accession number
        output_file: Path to save the genome
        
    Returns:
        Path to the downloaded genome file
    """
    output_path = os.path.join(OUTPUT_DIR, output_file)
    
    # Check if file already exists
    if os.path.exists(output_path):
        print(f"Genome file already exists at {output_path}")
        return output_path
    
    print(f"Downloading genome {accession}...")
    
    try:
        # Fetch the record
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        
        # Save to file
        SeqIO.write(record, output_path, "fasta")
        print(f"Genome saved to {output_path}")
        
        return output_path
    except Exception as e:
        print(f"Error downloading genome: {e}")
        
        # Fallback to direct download if Entrez fails
        try:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession}&rettype=fasta&retmode=text"
            response = requests.get(url)
            response.raise_for_status()
            
            with open(output_path, "w") as f:
                f.write(response.text)
            
            print(f"Genome saved to {output_path} using direct download")
            return output_path
        except Exception as e2:
            print(f"Error with direct download: {e2}")
            raise


def extract_regions(genome_file: str, num_regions: int = 5, 
                   min_length: int = 100, max_length: int = 500) -> List[SeqRecord]:
    """
    Extract random regions from a genome.
    
    Args:
        genome_file: Path to the genome file
        num_regions: Number of regions to extract
        min_length: Minimum region length
        max_length: Maximum region length
        
    Returns:
        List of SeqRecord objects
    """
    # Set random seed for reproducibility
    random.seed(SEED)
    
    # Load the genome
    genome = SeqIO.read(genome_file, "fasta")
    genome_length = len(genome.seq)
    
    regions = []
    for i in range(num_regions):
        # Generate random region length
        region_length = random.randint(min_length, max_length)
        
        # Generate random start position
        max_start = genome_length - region_length
        start = random.randint(0, max_start)
        
        # Extract the region
        region_seq = genome.seq[start:start+region_length]
        
        # Create a SeqRecord
        region_record = SeqRecord(
            seq=region_seq,
            id=f"region_{i+1}",
            description=f"Region {i+1} from {genome.id} position {start}-{start+region_length}"
        )
        
        regions.append(region_record)
    
    return regions


def create_variants(regions: List[SeqRecord]) -> Dict[str, List[SeqRecord]]:
    """
    Create variants of each region for testing.
    
    Args:
        regions: List of original region SeqRecords
        
    Returns:
        Dictionary mapping variant type to list of variant SeqRecords
    """
    # Set random seed for reproducibility
    random.seed(SEED)
    
    variants = {
        "snp": [],
        "insertion": [],
        "deletion": [],
        "complex": []
    }
    
    for region in regions:
        # Create SNP variant (1% of bases changed)
        snp_seq = create_snp_variant(region.seq, rate=0.01)
        snp_record = SeqRecord(
            seq=snp_seq,
            id=f"{region.id}_snp",
            description=f"SNP variant of {region.id}"
        )
        variants["snp"].append(snp_record)
        
        # Create insertion variant (1-3 insertions of 1-10bp)
        ins_seq = create_insertion_variant(region.seq, num_insertions=random.randint(1, 3))
        ins_record = SeqRecord(
            seq=ins_seq,
            id=f"{region.id}_insertion",
            description=f"Insertion variant of {region.id}"
        )
        variants["insertion"].append(ins_record)
        
        # Create deletion variant (1-3 deletions of 1-10bp)
        del_seq = create_deletion_variant(region.seq, num_deletions=random.randint(1, 3))
        del_record = SeqRecord(
            seq=del_seq,
            id=f"{region.id}_deletion",
            description=f"Deletion variant of {region.id}"
        )
        variants["deletion"].append(del_record)
        
        # Create complex variant (combination of SNPs, insertions, and deletions)
        complex_seq = create_complex_variant(region.seq)
        complex_record = SeqRecord(
            seq=complex_seq,
            id=f"{region.id}_complex",
            description=f"Complex variant of {region.id}"
        )
        variants["complex"].append(complex_record)
    
    return variants


def create_snp_variant(seq: Seq, rate: float = 0.01) -> Seq:
    """
    Create a variant with single nucleotide polymorphisms.
    
    Args:
        seq: Original sequence
        rate: Rate of SNPs (0.01 = 1% of bases changed)
        
    Returns:
        Sequence with SNPs
    """
    bases = list(str(seq))
    num_snps = int(len(bases) * rate)
    
    # Ensure at least one SNP
    num_snps = max(1, num_snps)
    
    # Select positions for SNPs
    positions = random.sample(range(len(bases)), num_snps)
    
    # Apply SNPs
    for pos in positions:
        original_base = bases[pos]
        # Choose a different base
        new_base = random.choice([b for b in "ACGT" if b != original_base])
        bases[pos] = new_base
    
    return Seq("".join(bases))


def create_insertion_variant(seq: Seq, num_insertions: int = 2, 
                            min_length: int = 1, max_length: int = 10) -> Seq:
    """
    Create a variant with insertions.
    
    Args:
        seq: Original sequence
        num_insertions: Number of insertions to create
        min_length: Minimum insertion length
        max_length: Maximum insertion length
        
    Returns:
        Sequence with insertions
    """
    bases = list(str(seq))
    
    # Select positions for insertions
    positions = sorted(random.sample(range(len(bases)), num_insertions))
    
    # Apply insertions (starting from the end to avoid position shifts)
    for pos in reversed(positions):
        # Generate random insertion
        ins_length = random.randint(min_length, max_length)
        insertion = "".join(random.choices("ACGT", k=ins_length))
        
        # Insert at position
        bases.insert(pos, insertion)
    
    return Seq("".join(bases))


def create_deletion_variant(seq: Seq, num_deletions: int = 2,
                           min_length: int = 1, max_length: int = 10) -> Seq:
    """
    Create a variant with deletions.
    
    Args:
        seq: Original sequence
        num_deletions: Number of deletions to create
        min_length: Minimum deletion length
        max_length: Maximum deletion length
        
    Returns:
        Sequence with deletions
    """
    bases = list(str(seq))
    seq_length = len(bases)
    
    # Ensure we don't delete too much
    max_total_deletion = seq_length // 3
    max_length = min(max_length, max_total_deletion // num_deletions)
    
    # Create deletion regions (non-overlapping)
    regions = []
    for _ in range(num_deletions):
        while True:
            # Generate random deletion length
            del_length = random.randint(min_length, max_length)
            
            # Generate random start position
            max_start = seq_length - del_length
            if max_start <= 0:
                break
                
            start = random.randint(0, max_start)
            end = start + del_length
            
            # Check for overlap with existing regions
            if not any(r[0] < end and start < r[1] for r in regions):
                regions.append((start, end))
                break
    
    # Sort regions by start position (descending)
    regions.sort(reverse=True)
    
    # Apply deletions
    for start, end in regions:
        del bases[start:end]
    
    return Seq("".join(bases))


def create_complex_variant(seq: Seq) -> Seq:
    """
    Create a complex variant with SNPs, insertions, and deletions.
    
    Args:
        seq: Original sequence
        
    Returns:
        Sequence with complex changes
    """
    # Apply SNPs (0.5% of bases)
    variant_seq = create_snp_variant(seq, rate=0.005)
    
    # Apply insertions (1-2 insertions)
    variant_seq = create_insertion_variant(variant_seq, num_insertions=random.randint(1, 2))
    
    # Apply deletions (1-2 deletions)
    variant_seq = create_deletion_variant(variant_seq, num_deletions=random.randint(1, 2))
    
    return variant_seq


def save_sequences(sequences: List[SeqRecord], output_file: str):
    """
    Save sequences to a FASTA file.
    
    Args:
        sequences: List of SeqRecord objects
        output_file: Path to save the sequences
    """
    output_path = os.path.join(OUTPUT_DIR, output_file)
    SeqIO.write(sequences, output_path, "fasta")
    print(f"Saved {len(sequences)} sequences to {output_path}")


def save_variant_sequences(variants: Dict[str, List[SeqRecord]]):
    """
    Save variant sequences to separate FASTA files.
    
    Args:
        variants: Dictionary mapping variant type to list of variant SeqRecords
    """
    for variant_type, records in variants.items():
        output_file = f"{variant_type}_variants.fasta"
        save_sequences(records, output_file)


def main():
    """Main function to prepare test sequences."""
    parser = argparse.ArgumentParser(description="Prepare test sequences for minimap2 wrapper testing")
    parser.add_argument("--force", action="store_true", help="Force regeneration of test data")
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Check if test data already exists
    reference_file = os.path.join(OUTPUT_DIR, "reference.fasta")
    regions_file = os.path.join(OUTPUT_DIR, "regions.fasta")
    
    if os.path.exists(reference_file) and os.path.exists(regions_file) and not args.force:
        print("Test data already exists. Use --force to regenerate.")
        return
    
    # Download the PhiX174 genome
    genome_file = download_genome(PHIX174_ACCESSION, "reference.fasta")
    
    # Extract regions
    regions = extract_regions(genome_file)
    save_sequences(regions, "regions.fasta")
    
    # Create variants
    variants = create_variants(regions)
    save_variant_sequences(variants)
    
    print("Test data preparation complete.")


if __name__ == "__main__":
    main()
