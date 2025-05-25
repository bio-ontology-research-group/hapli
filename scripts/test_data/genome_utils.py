#!/usr/bin/env python3
"""
Utility functions for generating synthetic genomic data.

This module provides reusable functions for creating realistic synthetic
genomic sequences with various characteristics including:
- Weighted random DNA sequence generation
- Repetitive element insertion
- Pattern-based sequence creation
- Sequence statistics calculation
- Realistic chromosome length generation

These utilities can be used by other scripts in the test data generation
pipeline to create consistent and realistic synthetic genomic data.
"""

import random
import re
from typing import Dict, List, Optional, Tuple, Union
from collections import Counter


def get_random_sequence(length: int, gc_content: float = 0.5, seed: Optional[int] = None) -> str:
    """
    Generate a random DNA sequence with specified GC content.
    
    Args:
        length: Length of sequence to generate in base pairs
        gc_content: Target GC content as fraction (0.0 to 1.0)
        seed: Random seed for reproducibility (optional)
        
    Returns:
        Random DNA sequence string containing A, T, G, C
        
    Raises:
        ValueError: If length < 1 or gc_content not in [0.0, 1.0]
        
    Example:
        >>> seq = get_random_sequence(100, gc_content=0.6, seed=42)
        >>> len(seq)
        100
    """
    if length < 1:
        raise ValueError("Length must be at least 1")
    
    if not 0.0 <= gc_content <= 1.0:
        raise ValueError("GC content must be between 0.0 and 1.0")
    
    # Set seed if provided
    if seed is not None:
        random.seed(seed)
    
    # Calculate nucleotide probabilities
    gc_prob = gc_content / 2  # Split GC content between G and C
    at_prob = (1 - gc_content) / 2  # Split AT content between A and T
    
    nucleotides = ['A', 'T', 'G', 'C']
    weights = [at_prob, at_prob, gc_prob, gc_prob]
    
    sequence = random.choices(nucleotides, weights=weights, k=length)
    return ''.join(sequence)


def insert_repeats(sequence: str, repeat_unit: str, num_copies: int, 
                  positions: List[int]) -> str:
    """
    Insert repetitive elements at specified positions in a sequence.
    
    Args:
        sequence: Original DNA sequence
        repeat_unit: DNA sequence to repeat (e.g., "TTAGGG", "CA")
        num_copies: Number of copies of the repeat unit to insert
        positions: List of positions where to insert repeats
        
    Returns:
        Modified sequence with repeats inserted
        
    Raises:
        ValueError: If positions are out of bounds or repeat_unit is empty
        
    Example:
        >>> seq = "ATCGATCG"
        >>> result = insert_repeats(seq, "CA", 3, [2, 6])
        >>> "CACACA" in result
        True
    """
    if not repeat_unit:
        raise ValueError("Repeat unit cannot be empty")
    
    if not all(isinstance(pos, int) and 0 <= pos <= len(sequence) for pos in positions):
        raise ValueError("All positions must be valid indices within sequence")
    
    # Create the repeat sequence
    repeat_sequence = repeat_unit * num_copies
    
    # Sort positions in reverse order to maintain correct indices during insertion
    sorted_positions = sorted(positions, reverse=True)
    
    # Convert sequence to list for easier manipulation
    seq_list = list(sequence)
    
    # Insert repeats at each position
    for pos in sorted_positions:
        seq_list[pos:pos] = list(repeat_sequence)
    
    return ''.join(seq_list)


def create_patterned_sequence(length: int, patterns_dict: Dict[str, float]) -> str:
    """
    Create a DNA sequence with specific patterns at given frequencies.
    
    Args:
        length: Total length of sequence to generate
        patterns_dict: Dictionary mapping DNA patterns to their frequencies
                      Frequencies should sum to <= 1.0, remainder filled with random bases
        
    Returns:
        DNA sequence with specified patterns
        
    Raises:
        ValueError: If length < 1 or frequencies sum > 1.0 or patterns are invalid
        
    Example:
        >>> patterns = {"TTAGGG": 0.1, "CACA": 0.05, "NNNN": 0.02}
        >>> seq = create_patterned_sequence(1000, patterns)
        >>> len(seq)
        1000
    """
    if length < 1:
        raise ValueError("Length must be at least 1")
    
    # Validate patterns contain only valid DNA bases
    valid_bases = set('ATGCN')
    for pattern in patterns_dict.keys():
        if not pattern or not all(base.upper() in valid_bases for base in pattern):
            raise ValueError(f"Invalid pattern: {pattern}")
    
    # Check that frequencies sum to <= 1.0
    total_frequency = sum(patterns_dict.values())
    if total_frequency > 1.0:
        raise ValueError("Pattern frequencies cannot sum to more than 1.0")
    
    # Calculate how much sequence each pattern should occupy
    pattern_lengths = {}
    remaining_length = length
    
    for pattern, frequency in patterns_dict.items():
        pattern_length = int(length * frequency)
        pattern_lengths[pattern] = pattern_length
        remaining_length -= pattern_length
    
    # Build the sequence
    sequence_parts = []
    
    # Add patterned regions
    for pattern, pattern_length in pattern_lengths.items():
        if pattern_length > 0:
            # Calculate how many complete pattern units fit
            pattern_units = pattern_length // len(pattern)
            remainder = pattern_length % len(pattern)
            
            # Add complete pattern units
            sequence_parts.extend([pattern] * pattern_units)
            
            # Add partial pattern if needed
            if remainder > 0:
                sequence_parts.append(pattern[:remainder])
    
    # Fill remaining length with random sequence
    if remaining_length > 0:
        random_seq = get_random_sequence(remaining_length, gc_content=0.5)
        sequence_parts.append(random_seq)
    
    # Combine all parts and shuffle to distribute patterns
    full_sequence = ''.join(sequence_parts)
    
    # Convert to list, shuffle, and rejoin
    seq_list = list(full_sequence)
    random.shuffle(seq_list)
    
    return ''.join(seq_list)[:length]  # Ensure exact length


def get_sequence_stats(sequence: str) -> Dict[str, Union[float, int, Dict[str, int]]]:
    """
    Calculate comprehensive statistics for a DNA sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Dictionary containing:
        - length: Total sequence length
        - gc_content: GC content as percentage (0-100)
        - n_content: N content as percentage (0-100)
        - base_counts: Dictionary of individual base counts
        - repeat_content: Estimated repetitive content percentage
        - dinucleotide_freq: Frequency of all dinucleotides
        
    Example:
        >>> stats = get_sequence_stats("ATCGATCGNNNN")
        >>> stats['gc_content']
        50.0
        >>> stats['n_content']
        33.33
    """
    if not sequence:
        return {
            'length': 0,
            'gc_content': 0.0,
            'n_content': 0.0,
            'base_counts': {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0},
            'repeat_content': 0.0,
            'dinucleotide_freq': {}
        }
    
    sequence = sequence.upper()
    length = len(sequence)
    
    # Count individual bases
    base_counts = Counter(sequence)
    
    # Ensure all bases are represented
    for base in ['A', 'T', 'G', 'C', 'N']:
        if base not in base_counts:
            base_counts[base] = 0
    
    # Calculate GC content (excluding Ns)
    gc_count = base_counts['G'] + base_counts['C']
    non_n_count = length - base_counts['N']
    gc_content = (gc_count / non_n_count * 100) if non_n_count > 0 else 0.0
    
    # Calculate N content
    n_content = (base_counts['N'] / length * 100) if length > 0 else 0.0
    
    # Calculate dinucleotide frequencies
    dinucleotide_freq = {}
    if length > 1:
        dinucleotides = [sequence[i:i+2] for i in range(length-1)]
        dinuc_counts = Counter(dinucleotides)
        total_dinucs = len(dinucleotides)
        
        for dinuc, count in dinuc_counts.items():
            dinucleotide_freq[dinuc] = count / total_dinucs
    
    # Estimate repeat content by looking for simple tandem repeats
    repeat_content = _estimate_repeat_content(sequence)
    
    return {
        'length': length,
        'gc_content': round(gc_content, 2),
        'n_content': round(n_content, 2),
        'base_counts': dict(base_counts),
        'repeat_content': round(repeat_content, 2),
        'dinucleotide_freq': dinucleotide_freq
    }


def get_chromosome_lengths(num_chr: int, base_length: int, 
                          decay_factor: float = 0.9) -> List[Tuple[str, int]]:
    """
    Generate realistic chromosome lengths with special handling for sex chromosomes.
    
    Args:
        num_chr: Number of chromosomes to generate
        base_length: Length of chromosome 1 in base pairs
        decay_factor: Factor by which each chromosome is smaller than the previous
        
    Returns:
        List of tuples containing (chromosome_name, length)
        
    Raises:
        ValueError: If num_chr < 1 or base_length < 1000
        
    Example:
        >>> lengths = get_chromosome_lengths(5, 100000, 0.9)
        >>> lengths[0]
        ('chr1', 100000)
        >>> lengths[1][1] < lengths[0][1]  # chr2 shorter than chr1
        True
    """
    if num_chr < 1:
        raise ValueError("Number of chromosomes must be at least 1")
    
    if base_length < 1000:
        raise ValueError("Base length must be at least 1000 bp")
    
    chromosome_data = []
    current_length = base_length
    
    # Generate regular autosomes (chr1-chr22)
    regular_count = min(num_chr, 22)
    for i in range(1, regular_count + 1):
        chromosome_data.append((f"chr{i}", int(current_length)))
        current_length *= decay_factor
    
    remaining = num_chr - regular_count
    
    # Add special chromosomes if needed
    if remaining > 0:
        # chrX: typically ~85% the length of chr1
        chrx_length = int(base_length * 0.85)
        chromosome_data.append(("chrX", chrx_length))
        remaining -= 1
    
    if remaining > 0:
        # chrY: typically ~30% the length of chrX
        chry_length = int(base_length * 0.85 * 0.3)
        chromosome_data.append(("chrY", chry_length))
        remaining -= 1
    
    if remaining > 0:
        # chrM: mitochondrial chromosome, typically ~16kb
        chrm_length = 16569  # Actual human mitochondrial genome length
        chromosome_data.append(("chrM", chrm_length))
        remaining -= 1
    
    # Add any additional chromosomes with continued decay
    for i in range(remaining):
        chr_name = f"chr{23 + i}"
        chromosome_data.append((chr_name, int(current_length)))
        current_length *= decay_factor
    
    return chromosome_data


def _estimate_repeat_content(sequence: str) -> float:
    """
    Estimate the repetitive content of a sequence by detecting simple tandem repeats.
    
    This is a simplified estimation that looks for short tandem repeats (1-6 bp units).
    
    Args:
        sequence: DNA sequence to analyze
        
    Returns:
        Estimated percentage of sequence that is repetitive
    """
    if len(sequence) < 10:
        return 0.0
    
    repeat_bases = 0
    sequence = sequence.upper()
    
    # Look for tandem repeats of different unit sizes
    for unit_size in range(1, 7):  # 1-6 bp repeat units
        i = 0
        while i < len(sequence) - unit_size * 2:
            unit = sequence[i:i + unit_size]
            
            # Skip if unit contains N
            if 'N' in unit:
                i += 1
                continue
            
            # Count consecutive repeats of this unit
            repeat_count = 1
            j = i + unit_size
            
            while j + unit_size <= len(sequence) and sequence[j:j + unit_size] == unit:
                repeat_count += 1
                j += unit_size
            
            # If we found at least 3 consecutive repeats, count as repetitive
            if repeat_count >= 3:
                repeat_length = repeat_count * unit_size
                repeat_bases += repeat_length
                i = j  # Skip past this repeat region
            else:
                i += 1
    
    # Avoid double-counting overlapping repeats by capping at sequence length
    repeat_bases = min(repeat_bases, len(sequence))
    
    return (repeat_bases / len(sequence)) * 100 if len(sequence) > 0 else 0.0


# Convenience functions for common patterns
def create_telomeric_sequence(length: int, degenerate_rate: float = 0.05) -> str:
    """
    Create a telomeric sequence with TTAGGG repeats and some degenerate copies.
    
    Args:
        length: Length of telomeric sequence to generate
        degenerate_rate: Fraction of bases to mutate (0.0 to 1.0)
        
    Returns:
        Telomeric sequence string
    """
    telomere_motif = "TTAGGG"
    num_repeats = length // len(telomere_motif)
    base_sequence = telomere_motif * num_repeats
    
    # Add remainder
    remainder = length % len(telomere_motif)
    if remainder > 0:
        base_sequence += telomere_motif[:remainder]
    
    # Add degenerate mutations
    if degenerate_rate > 0:
        seq_list = list(base_sequence)
        for i in range(len(seq_list)):
            if random.random() < degenerate_rate:
                seq_list[i] = random.choice(['A', 'T', 'G', 'C'])
        base_sequence = ''.join(seq_list)
    
    return base_sequence


def create_cpg_island(length: int, gc_content: float = 0.75) -> str:
    """
    Create a CpG island sequence with elevated GC content.
    
    Args:
        length: Length of CpG island to generate
        gc_content: Target GC content (typically 0.7-0.8)
        
    Returns:
        CpG island sequence string
    """
    return get_random_sequence(length, gc_content=gc_content)


def create_centromeric_sequence(length: int) -> str:
    """
    Create a centromeric sequence with alpha satellite-like repeats and N regions.
    
    Args:
        length: Length of centromeric sequence to generate
        
    Returns:
        Centromeric sequence string
    """
    patterns = {
        "ATATACATAG": 0.15,
        "CATTCCATTC": 0.10,
        "AATCAACCC": 0.08,
        "TTCCATTCCATTC": 0.12,
        "NNNN": 0.05
    }
    
    return create_patterned_sequence(length, patterns)
