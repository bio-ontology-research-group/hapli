#!/usr/bin/env python3
"""
Utility functions for generating synthetic genomic data.
"""

import random
import logging
from typing import List, Tuple, Dict, Union, Optional
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Helper from generate_reference_genome.py

def generate_weighted_sequence(length: int, gc_content: float) -> str:
    """
    Generate a random DNA sequence with target GC content.
    """
    gc_prob = gc_content / 2
    at_prob = (1 - gc_content) / 2
    
    nucleotides = ['A', 'T', 'G', 'C']
    weights = [at_prob, at_prob, gc_prob, gc_prob]
    
    sequence = random.choices(nucleotides, weights=weights, k=length)
    return ''.join(sequence)

def add_cpg_islands(sequence: str, chr_length: int) -> str:
    seq_list = list(sequence)
    num_islands = random.randint(2, 5)
    search_region = int(chr_length * 0.1)
    
    for _ in range(num_islands):
        island_length = random.randint(500, 2000)
        if search_region < island_length: continue
            
        start_pos = random.randint(0, search_region - island_length)
        high_gc_content = random.uniform(0.70, 0.80)
        island_seq = generate_weighted_sequence(island_length, high_gc_content)
        
        for i, nucleotide in enumerate(island_seq):
            if start_pos + i < len(seq_list):
                seq_list[start_pos + i] = nucleotide
    return ''.join(seq_list)

def add_telomeric_repeats(sequence: str) -> str:
    seq_list = list(sequence)
    seq_length = len(sequence)
    num_repeats_start = random.randint(50, 200)
    num_repeats_end = random.randint(50, 200)
    telomere_motif = "TTAGGG"
    
    start_length = min(num_repeats_start * len(telomere_motif), seq_length // 4)
    start_seq = (telomere_motif * (start_length // len(telomere_motif) + 1))[:start_length]
    
    # Degenerate
    start_seq_list = list(start_seq)
    for i in range(len(start_seq_list)):
        if random.random() < 0.05:
            start_seq_list[i] = random.choice(['A', 'T', 'G', 'C'])
    start_seq = ''.join(start_seq_list)
    
    for i, nucleotide in enumerate(start_seq):
        if i < len(seq_list): seq_list[i] = nucleotide
    
    end_length = min(num_repeats_end * len(telomere_motif), seq_length // 4)
    end_seq = (telomere_motif * (end_length // len(telomere_motif) + 1))[:end_length]
    
    end_seq_list = list(end_seq)
    for i in range(len(end_seq_list)):
        if random.random() < 0.05:
            end_seq_list[i] = random.choice(['A', 'T', 'G', 'C'])
    end_seq = ''.join(end_seq_list)
    
    start_pos = seq_length - end_length
    for i, nucleotide in enumerate(end_seq):
        if start_pos + i < len(seq_list):
            seq_list[start_pos + i] = nucleotide
            
    return ''.join(seq_list)

def add_centromeric_region(sequence: str, chr_length: int) -> str:
    seq_list = list(sequence)
    centromere_fraction = random.uniform(0.01, 0.05)
    centromere_length = int(chr_length * centromere_fraction)
    center_position = random.uniform(0.40, 0.60)
    start_pos = int(chr_length * center_position - centromere_length // 2)
    start_pos = max(0, min(start_pos, chr_length - centromere_length))
    
    centromere_motifs = ["ATATACATAG", "CATTCCATTC", "AATCAACCC", "TTCCATTCCATTC"]
    centromere_seq = []
    pos = 0
    while pos < centromere_length:
        if random.random() < 0.3:
            n_length = min(random.randint(5, 50), centromere_length - pos)
            centromere_seq.extend(['N'] * n_length)
            pos += n_length
        else:
            motif = random.choice(centromere_motifs)
            repeat_length = min(random.randint(20, 200), centromere_length - pos)
            repeat_seq = (motif * (repeat_length // len(motif) + 1))[:repeat_length]
            centromere_seq.extend(list(repeat_seq))
            pos += repeat_length
            
    for i, nucleotide in enumerate(centromere_seq):
        if start_pos + i < len(seq_list):
            seq_list[start_pos + i] = nucleotide
            
    return ''.join(seq_list)

def add_repetitive_regions(sequence: str, repeat_fraction: float = 0.075) -> str:
    seq_list = list(sequence)
    seq_length = len(sequence)
    repeat_length = int(seq_length * repeat_fraction)
    repeat_motifs = ["CA", "GT", "ATCG", "TAGA", "CACA", "GTGT"]
    
    num_regions = random.randint(3, 8)
    for _ in range(num_regions):
        if repeat_length <= 0: break
        region_length = min(random.randint(50, 500), repeat_length)
        start_pos = random.randint(0, seq_length - region_length)
        motif = random.choice(repeat_motifs)
        repeat_seq = (motif * (region_length // len(motif) + 1))[:region_length]
        
        for i, nucleotide in enumerate(repeat_seq):
            if start_pos + i < seq_length: seq_list[start_pos + i] = nucleotide
        repeat_length -= region_length
    return ''.join(seq_list)

def add_n_regions(sequence: str, n_fraction: float = 0.015) -> str:
    seq_list = list(sequence)
    seq_length = len(sequence)
    n_length = int(seq_length * n_fraction)
    num_regions = random.randint(2, 6)
    for _ in range(num_regions):
        if n_length <= 0: break
        region_length = min(random.randint(10, 200), n_length)
        start_pos = random.randint(0, seq_length - region_length)
        for i in range(region_length):
            if start_pos + i < seq_length: seq_list[start_pos + i] = 'N'
        n_length -= region_length
    return ''.join(seq_list)

def generate_chromosome(name: str, length: int, gc_content: float, 
                       include_cpg_islands: bool = True, include_telomeres: bool = True, 
                       include_centromeres: bool = True) -> SeqRecord:
    logging.info(f"Generating {name} ({length:,} bp)")
    sequence = generate_weighted_sequence(length, gc_content)
    if include_cpg_islands: sequence = add_cpg_islands(sequence, length)
    if include_telomeres: sequence = add_telomeric_repeats(sequence)
    if include_centromeres: sequence = add_centromeric_region(sequence, length)
    sequence = add_repetitive_regions(sequence)
    sequence = add_n_regions(sequence)
    
    chr_num = name.replace('chr', '')
    desc = f"Homo sapiens chromosome {chr_num}, Test Assembly"
    
    return SeqRecord(Seq(sequence), id=name, description=f"{desc} (length={len(sequence)})")

def get_chromosome_lengths(num_chr: int, base_length: int, decay_factor: float = 0.9) -> List[int]:
    lengths = []
    current = base_length
    for i in range(num_chr):
        lengths.append(int(current))
        current *= decay_factor
    return lengths