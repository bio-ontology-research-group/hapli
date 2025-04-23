"""
Test data generator for alignment tests.

This module provides utilities to generate test data for alignment tests,
including mock minimap2 alignment objects and various alignment scenarios.
"""

import os
import random
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Any
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Mock class to simulate minimap2 alignment objects
@dataclass
class MockMinimapAlignment:
    """Mock class to simulate minimap2 alignment objects."""
    query_name: str
    ctg: str  # target/contig name
    q_st: int  # query start
    q_en: int  # query end
    strand: int  # strand: +1 for forward, -1 for reverse
    r_st: int  # target start
    r_en: int  # target end
    mlen: int  # number of matching bases
    blen: int  # alignment block length
    mapq: int  # mapping quality
    cigar_str: str  # CIGAR string
    
    @property
    def is_reverse(self) -> bool:
        """Return True if alignment is on reverse strand."""
        return self.strand == -1


class AlignmentTestData:
    """
    Generator for test alignment data.
    
    This class provides methods to create various test alignment scenarios
    for unit testing alignment result processing.
    """
    
    @staticmethod
    def generate_random_sequence(length: int) -> str:
        """
        Generate a random DNA sequence.
        
        Args:
            length: Length of the sequence to generate
            
        Returns:
            Random DNA sequence string
        """
        bases = ['A', 'C', 'G', 'T']
        return ''.join(random.choice(bases) for _ in range(length))
    
    @staticmethod
    def create_perfect_match_alignment(seq_length: int = 100) -> Tuple[MockMinimapAlignment, str, str]:
        """
        Create a perfect match alignment scenario.
        
        Args:
            seq_length: Length of the sequences
            
        Returns:
            Tuple of (alignment, query_sequence, target_sequence)
        """
        # Generate a random sequence
        sequence = AlignmentTestData.generate_random_sequence(seq_length)
        
        # Create alignment with perfect match
        alignment = MockMinimapAlignment(
            query_name="query1",
            ctg="target1",
            q_st=0,
            q_en=seq_length,
            strand=1,  # forward
            r_st=0,
            r_en=seq_length,
            mlen=seq_length,  # all bases match
            blen=seq_length,
            mapq=60,  # high mapping quality
            cigar_str=f"{seq_length}M"  # all matches
        )
        
        return alignment, sequence, sequence
    
    @staticmethod
    def create_snp_alignment(seq_length: int = 100, num_snps: int = 5) -> Tuple[MockMinimapAlignment, str, str]:
        """
        Create an alignment with SNPs.
        
        Args:
            seq_length: Length of the sequences
            num_snps: Number of SNPs to introduce
            
        Returns:
            Tuple of (alignment, query_sequence, target_sequence)
        """
        # Generate a random sequence
        target_seq = AlignmentTestData.generate_random_sequence(seq_length)
        
        # Create a copy with SNPs
        query_seq_list = list(target_seq)
        
        # Ensure num_snps doesn't exceed sequence length
        num_snps = min(num_snps, seq_length)
        
        # Select random positions for SNPs
        snp_positions = random.sample(range(seq_length), num_snps)
        
        # Introduce SNPs
        bases = {'A': ['C', 'G', 'T'], 'C': ['A', 'G', 'T'], 
                'G': ['A', 'C', 'T'], 'T': ['A', 'C', 'G']}
        
        for pos in snp_positions:
            original_base = query_seq_list[pos]
            query_seq_list[pos] = random.choice(bases[original_base])
        
        query_seq = ''.join(query_seq_list)
        
        # Create alignment with SNPs
        alignment = MockMinimapAlignment(
            query_name="query1",
            ctg="target1",
            q_st=0,
            q_en=seq_length,
            strand=1,  # forward
            r_st=0,
            r_en=seq_length,
            mlen=seq_length - num_snps,  # matches minus SNPs
            blen=seq_length,
            mapq=60 - num_snps,  # quality decreases with SNPs
            cigar_str=f"{seq_length}M"  # still all matches in CIGAR
        )
        
        return alignment, query_seq, target_seq
    
    @staticmethod
    def create_insertion_alignment(seq_length: int = 100, 
                                  ins_length: int = 10) -> Tuple[MockMinimapAlignment, str, str]:
        """
        Create an alignment with an insertion in the query.
        
        Args:
            seq_length: Length of the target sequence
            ins_length: Length of the insertion
            
        Returns:
            Tuple of (alignment, query_sequence, target_sequence)
        """
        # Generate a random sequence
        target_seq = AlignmentTestData.generate_random_sequence(seq_length)
        
        # Create a position for insertion (not at the very beginning or end)
        ins_pos = random.randint(10, seq_length - 10)
        
        # Generate insertion sequence
        insertion = AlignmentTestData.generate_random_sequence(ins_length)
        
        # Create query sequence with insertion
        query_seq = target_seq[:ins_pos] + insertion + target_seq[ins_pos:]
        
        # Create alignment with insertion
        alignment = MockMinimapAlignment(
            query_name="query1",
            ctg="target1",
            q_st=0,
            q_en=seq_length + ins_length,
            strand=1,  # forward
            r_st=0,
            r_en=seq_length,
            mlen=seq_length,  # matches without insertion
            blen=seq_length + ins_length,
            mapq=50,  # lower quality due to insertion
            cigar_str=f"{ins_pos}M{ins_length}I{seq_length - ins_pos}M"
        )
        
        return alignment, query_seq, target_seq
    
    @staticmethod
    def create_deletion_alignment(seq_length: int = 100, 
                                del_length: int = 10) -> Tuple[MockMinimapAlignment, str, str]:
        """
        Create an alignment with a deletion in the query.
        
        Args:
            seq_length: Length of the query sequence before deletion
            del_length: Length of the deletion
            
        Returns:
            Tuple of (alignment, query_sequence, target_sequence)
        """
        # Generate a random sequence for the full target
        target_seq = AlignmentTestData.generate_random_sequence(seq_length + del_length)
        
        # Create a position for deletion (not at the very beginning or end)
        del_pos = random.randint(10, seq_length - 10)
        
        # Create query sequence with deletion
        query_seq = target_seq[:del_pos] + target_seq[del_pos + del_length:]
        
        # Create alignment with deletion
        alignment = MockMinimapAlignment(
            query_name="query1",
            ctg="target1",
            q_st=0,
            q_en=seq_length,
            strand=1,  # forward
            r_st=0,
            r_en=seq_length + del_length,
            mlen=seq_length,  # matches without deletion
            blen=seq_length + del_length,
            mapq=50,  # lower quality due to deletion
            cigar_str=f"{del_pos}M{del_length}D{seq_length - del_pos}M"
        )
        
        return alignment, query_seq, target_seq
    
    @staticmethod
    def create_complex_alignment(seq_length: int = 100) -> Tuple[MockMinimapAlignment, str, str]:
        """
        Create a complex alignment with SNPs, insertions, and deletions.
        
        Args:
            seq_length: Base length of the sequences
            
        Returns:
            Tuple of (alignment, query_sequence, target_sequence)
        """
        # Generate a random sequence
        target_seq = AlignmentTestData.generate_random_sequence(seq_length)
        
        # Create a copy for modification
        query_seq_list = list(target_seq)
        
        # Introduce SNPs (5% of bases)
        num_snps = max(1, int(seq_length * 0.05))
        snp_positions = random.sample(range(seq_length), num_snps)
        
        bases = {'A': ['C', 'G', 'T'], 'C': ['A', 'G', 'T'], 
                'G': ['A', 'C', 'T'], 'T': ['A', 'C', 'G']}
        
        for pos in snp_positions:
            original_base = query_seq_list[pos]
            query_seq_list[pos] = random.choice(bases[original_base])
        
        # Create a small insertion
        ins_pos = random.randint(20, seq_length - 20)
        ins_length = random.randint(3, 8)
        insertion = AlignmentTestData.generate_random_sequence(ins_length)
        
        # Apply insertion
        query_seq_list = query_seq_list[:ins_pos] + list(insertion) + query_seq_list[ins_pos:]
        
        # Create a small deletion
        del_pos = random.randint(ins_pos + 10, min(len(query_seq_list) - 10, ins_pos + 30))
        del_length = random.randint(3, 8)
        
        # Apply deletion (if it doesn't go beyond the end)
        if del_pos + del_length < len(query_seq_list):
            query_seq_list = query_seq_list[:del_pos] + query_seq_list[del_pos + del_length:]
        
        query_seq = ''.join(query_seq_list)
        
        # Create a complex CIGAR string
        cigar = f"{ins_pos}M{ins_length}I"
        
        # Calculate the distance from insertion end to deletion start
        distance = del_pos - ins_pos - ins_length
        if distance > 0:
            cigar += f"{distance}M"
        
        if del_pos + del_length < len(query_seq_list) + del_length:
            cigar += f"{del_length}D"
            remaining = len(query_seq_list) - del_pos
            if remaining > 0:
                cigar += f"{remaining}M"
        else:
            # If deletion would go beyond end, just add remaining matches
            remaining = len(query_seq_list) - ins_pos - ins_length
            if remaining > 0:
                cigar += f"{remaining}M"
        
        # Create alignment with complex pattern
        alignment = MockMinimapAlignment(
            query_name="query1",
            ctg="target1",
            q_st=0,
            q_en=len(query_seq),
            strand=1,  # forward
            r_st=0,
            r_en=seq_length,
            mlen=seq_length - num_snps - del_length,  # approximate matches
            blen=max(len(query_seq), seq_length),
            mapq=40,  # lower quality due to complexity
            cigar_str=cigar
        )
        
        return alignment, query_seq, target_seq
    
    @staticmethod
    def create_reverse_complement_alignment(seq_length: int = 100) -> Tuple[MockMinimapAlignment, str, str]:
        """
        Create an alignment where the query is the reverse complement of the target.
        
        Args:
            seq_length: Length of the sequences
            
        Returns:
            Tuple of (alignment, query_sequence, target_sequence)
        """
        # Generate a random sequence
        target_seq = AlignmentTestData.generate_random_sequence(seq_length)
        
        # Create reverse complement
        query_seq = str(Seq(target_seq).reverse_complement())
        
        # Create alignment with reverse complement
        alignment = MockMinimapAlignment(
            query_name="query1",
            ctg="target1",
            q_st=0,
            q_en=seq_length,
            strand=-1,  # reverse
            r_st=0,
            r_en=seq_length,
            mlen=seq_length,  # all bases match (after reverse complement)
            blen=seq_length,
            mapq=60,  # high mapping quality
            cigar_str=f"{seq_length}M"  # all matches
        )
        
        return alignment, query_seq, target_seq
    
    @staticmethod
    def create_no_alignment() -> Tuple[None, str, str]:
        """
        Create a scenario with no alignment (completely different sequences).
        
        Returns:
            Tuple of (None, query_sequence, target_sequence)
        """
        # Generate two completely different random sequences
        query_seq = AlignmentTestData.generate_random_sequence(100)
        
        # Ensure target is very different (use different base distribution)
        bases = {'A': 0.7, 'C': 0.1, 'G': 0.1, 'T': 0.1}
        target_seq = ''.join(random.choices(
            population=list(bases.keys()),
            weights=list(bases.values()),
            k=100
        ))
        
        # No alignment object
        return None, query_seq, target_seq
    
    @staticmethod
    def get_expected_results() -> Dict[str, Dict[str, Any]]:
        """
        Get expected results for different alignment scenarios.
        
        Returns:
            Dictionary mapping scenario names to expected result values
        """
        return {
            "perfect_match": {
                "identity": 1.0,
                "coverage": 1.0,
                "has_snps": False,
                "has_indels": False
            },
            "snp": {
                "identity_range": (0.9, 0.99),  # 90-99% identity
                "coverage": 1.0,
                "has_snps": True,
                "has_indels": False
            },
            "insertion": {
                "identity_range": (0.8, 1.0),
                "coverage_range": (0.9, 1.1),
                "has_snps": False,
                "has_indels": True
            },
            "deletion": {
                "identity_range": (0.8, 1.0),
                "coverage_range": (0.9, 1.1),
                "has_snps": False,
                "has_indels": True
            },
            "complex": {
                "identity_range": (0.7, 0.95),
                "coverage_range": (0.8, 1.2),
                "has_snps": True,
                "has_indels": True
            },
            "reverse": {
                "identity": 1.0,
                "coverage": 1.0,
                "is_reverse": True
            },
            "no_alignment": {
                "identity": 0.0,
                "coverage": 0.0
            }
        }
    
    @classmethod
    def create_all_test_scenarios(cls) -> Dict[str, Tuple]:
        """
        Create all test alignment scenarios.
        
        Returns:
            Dictionary mapping scenario names to (alignment, query_seq, target_seq) tuples
        """
        return {
            "perfect_match": cls.create_perfect_match_alignment(),
            "snp": cls.create_snp_alignment(),
            "insertion": cls.create_insertion_alignment(),
            "deletion": cls.create_deletion_alignment(),
            "complex": cls.create_complex_alignment(),
            "reverse": cls.create_reverse_complement_alignment(),
            "no_alignment": cls.create_no_alignment()
        }
