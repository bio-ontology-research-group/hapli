"""
Alignment result processing module.

This module provides standardized data structures and utilities for handling
sequence alignment results, with a focus on minimap2 alignments.
"""

import json
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Tuple, Any, Union
import re
from enum import Enum
import math
from Bio.Seq import Seq  # For reverse complement operations

class AlignmentType(Enum):
    """Types of alignments based on their characteristics."""
    PERFECT = "perfect"
    HIGH_IDENTITY = "high_identity"
    MODERATE_IDENTITY = "moderate_identity"
    LOW_IDENTITY = "low_identity"
    NO_ALIGNMENT = "no_alignment"


@dataclass
class CigarOperation:
    """Represents a single CIGAR operation."""
    operation: str  # M, I, D, N, S, H, P, =, X
    length: int
    
    def __str__(self) -> str:
        return f"{self.length}{self.operation}"


@dataclass
class AlignmentStatistics:
    """Statistical metrics for an alignment."""
    identity: float = 0.0  # Percentage of matching bases
    coverage: float = 0.0  # Percentage of query covered by alignment
    matches: int = 0       # Number of matching positions
    mismatches: int = 0    # Number of mismatching positions
    insertions: int = 0    # Number of insertions
    deletions: int = 0     # Number of deletions
    gaps: int = 0          # Total number of gaps (insertions + deletions)
    alignment_length: int = 0  # Total length of alignment
    query_length: int = 0  # Length of query sequence
    target_length: int = 0 # Length of target sequence
    
    @property
    def alignment_type(self) -> AlignmentType:
        """Determine the type of alignment based on statistics."""
        if self.identity == 0 and self.coverage == 0:
            return AlignmentType.NO_ALIGNMENT
        elif self.identity >= 0.98:
            return AlignmentType.PERFECT
        elif self.identity >= 0.90:
            return AlignmentType.HIGH_IDENTITY
        elif self.identity >= 0.70:
            return AlignmentType.MODERATE_IDENTITY
        else:
            return AlignmentType.LOW_IDENTITY


@dataclass
class AlignmentResult:
    """
    Comprehensive representation of a sequence alignment result.
    
    This class stores all relevant information about an alignment between
    a query sequence and a target sequence, including positions, statistics,
    and the alignment details.
    """
    # Basic alignment information
    query_name: str
    target_name: str
    query_sequence: Optional[str] = None
    target_sequence: Optional[str] = None
    
    # Alignment coordinates
    query_start: int = 0
    query_end: int = 0
    target_start: int = 0
    target_end: int = 0
    
    # Alignment details
    score: float = 0.0
    cigar_string: str = ""
    cigar_operations: List[CigarOperation] = field(default_factory=list)
    is_reverse: bool = False
    
    # Computed statistics
    statistics: AlignmentStatistics = field(default_factory=AlignmentStatistics)
    
    # Alignment visualization
    aligned_query: Optional[str] = None
    aligned_target: Optional[str] = None
    alignment_indicator: Optional[str] = None
    
    @classmethod
    def from_minimap2(cls, alignment, query_seq: str, target_seq: str) -> 'AlignmentResult':
        """
        Create an AlignmentResult from a minimap2 alignment object.
        
        Args:
            alignment: A minimap2 alignment object
            query_seq: The query sequence string
            target_seq: The target sequence string
            
        Returns:
            An AlignmentResult object
        """
        # Extract basic information
        result = cls(
            query_name=alignment.query_name if hasattr(alignment, 'query_name') else "query",
            target_name=alignment.ctg if hasattr(alignment, 'ctg') else "target",
            query_sequence=query_seq,
            target_sequence=target_seq,
            query_start=alignment.q_st,
            query_end=alignment.q_en,
            target_start=alignment.r_st,
            target_end=alignment.r_en,
            score=alignment.mapq if hasattr(alignment, 'mapq') else 0,
            cigar_string=alignment.cigar_str if hasattr(alignment, 'cigar_str') else "",
            is_reverse=alignment.is_reverse if hasattr(alignment, 'is_reverse') else False
        )
        
        # Parse CIGAR string
        result.cigar_operations = cls._parse_cigar(result.cigar_string)
        
        # Calculate statistics
        result._calculate_statistics()
        
        # Generate alignment visualization
        result._generate_alignment_visualization()
        
        return result
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'AlignmentResult':
        """
        Create an AlignmentResult from a dictionary.
        
        Args:
            data: Dictionary containing alignment data
            
        Returns:
            An AlignmentResult object
        """
        # Handle nested objects
        if 'statistics' in data:
            data['statistics'] = AlignmentStatistics(**data['statistics'])
        
        if 'cigar_operations' in data:
            data['cigar_operations'] = [CigarOperation(**op) for op in data['cigar_operations']]
        
        return cls(**data)
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the AlignmentResult to a dictionary.
        
        Returns:
            Dictionary representation of the alignment result
        """
        result = asdict(self)
        return result
    
    def to_json(self, indent: int = 2) -> str:
        """
        Convert the AlignmentResult to a JSON string.
        
        Args:
            indent: Number of spaces for indentation
            
        Returns:
            JSON string representation of the alignment result
        """
        return json.dumps(self.to_dict(), indent=indent)
    
    @staticmethod
    def _parse_cigar(cigar_string: str) -> List[CigarOperation]:
        """
        Parse a CIGAR string into a list of CigarOperation objects.
        
        Args:
            cigar_string: CIGAR string (e.g., "10M2I5M")
            
        Returns:
            List of CigarOperation objects
        """
        operations = []
        if not cigar_string:
            return operations
        
        # Regular expression to match CIGAR operations
        pattern = r'(\d+)([MIDNSHP=X])'
        matches = re.findall(pattern, cigar_string)
        
        for length_str, op in matches:
            operations.append(CigarOperation(
                operation=op,
                length=int(length_str)
            ))
        
        return operations
    
    def _calculate_statistics(self) -> None:
        """Calculate alignment statistics from CIGAR operations."""
        stats = self.statistics
        
        # Set sequence lengths
        stats.query_length = len(self.query_sequence) if self.query_sequence else 0
        stats.target_length = len(self.target_sequence) if self.target_sequence else 0
        
        # Process CIGAR operations
        matches = 0
        mismatches = 0
        insertions = 0
        deletions = 0
        
        for op in self.cigar_operations:
            if op.operation == 'M':  # Match or mismatch
                # For 'M', we need to check the actual sequences to determine matches vs mismatches
                matches += op.length  # Simplified - actual matches would need sequence comparison
            elif op.operation == '=':  # Match
                matches += op.length
            elif op.operation == 'X':  # Mismatch
                mismatches += op.length
            elif op.operation == 'I':  # Insertion
                insertions += op.length
            elif op.operation == 'D':  # Deletion
                deletions += op.length
        
        # Update statistics
        stats.matches = matches
        stats.mismatches = mismatches
        stats.insertions = insertions
        stats.deletions = deletions
        stats.gaps = insertions + deletions
        
        # Calculate alignment length
        stats.alignment_length = matches + mismatches + deletions
        
        # Calculate identity and coverage
        if stats.alignment_length > 0:
            # Identity is the percentage of matching bases relative to alignment length
            # For SNP alignments, we need to be more accurate about matches vs mismatches
            if self.query_sequence and self.target_sequence and 'M' in self.cigar_string:
                # For 'M' operations, we need to count actual matches by comparing sequences
                actual_matches = 0
                actual_mismatches = 0
                
                q_pos = self.query_start
                t_pos = self.target_start
                
                for op in self.cigar_operations:
                    if op.operation == 'M':
                        # Compare the sequences base by base
                        for i in range(op.length):
                            if q_pos + i < len(self.query_sequence) and t_pos + i < len(self.target_sequence):
                                if self.query_sequence[q_pos + i] == self.target_sequence[t_pos + i]:
                                    actual_matches += 1
                                else:
                                    actual_mismatches += 1
                        q_pos += op.length
                        t_pos += op.length
                    elif op.operation == 'I':
                        q_pos += op.length
                    elif op.operation == 'D':
                        t_pos += op.length
                
                # Update matches and mismatches with actual counts
                stats.matches = actual_matches
                stats.mismatches = actual_mismatches
                
            # Calculate identity using updated match/mismatch counts
            total = stats.matches + stats.mismatches + stats.insertions + stats.deletions
            if total > 0:
                stats.identity = stats.matches / total
            else:
                stats.identity = 0.0
                
            # Special case for reverse complement alignments
            if self.is_reverse and 'M' in self.cigar_string and len(self.cigar_operations) == 1:
                # For pure reverse complement alignments, we should have high identity
                # Check if sequences are actually reverse complements
                if self.query_sequence and self.target_sequence:
                    from Bio.Seq import Seq
                    rev_comp = str(Seq(self.target_sequence).reverse_complement())
                    if self.query_sequence == rev_comp:
                        stats.identity = 1.0
        
        if stats.query_length > 0:
            stats.coverage = (self.query_end - self.query_start) / stats.query_length
    
    def _generate_alignment_visualization(self) -> None:
        """
        Generate a text-based visualization of the alignment.
        
        This creates three strings:
        - aligned_query: The query sequence with gaps
        - aligned_target: The target sequence with gaps
        - alignment_indicator: Shows matches (|), mismatches (.), and gaps (space)
        """
        if not self.query_sequence or not self.target_sequence:
            return
            
        # For SNP alignments with just an 'M' operation, we need to create a simple visualization
        if len(self.cigar_operations) == 1 and self.cigar_operations[0].operation == 'M':
            q_segment = self.query_sequence[self.query_start:self.query_end]
            t_segment = self.target_sequence[self.target_start:self.target_end]
            
            # Create match indicators
            indicator = []
            for q, t in zip(q_segment, t_segment):
                if q == t:
                    indicator.append('|')  # Match
                else:
                    indicator.append('.')  # Mismatch
            
            self.aligned_query = q_segment
            self.aligned_target = t_segment
            self.alignment_indicator = ''.join(indicator)
            return
            
        # For more complex alignments with multiple operations
        if not self.cigar_operations:
            # If no CIGAR operations, create a simple visualization
            length = min(len(self.query_sequence), len(self.target_sequence))
            self.aligned_query = self.query_sequence[:length]
            self.aligned_target = self.target_sequence[:length]
            
            # Create match indicators
            indicator = []
            for q, t in zip(self.aligned_query, self.aligned_target):
                if q == t:
                    indicator.append('|')  # Match
                else:
                    indicator.append('.')  # Mismatch
            
            self.alignment_indicator = ''.join(indicator)
            return
        
        query_aligned = []
        target_aligned = []
        indicator = []
        
        query_pos = self.query_start
        target_pos = self.target_start
        
        for op in self.cigar_operations:
            if op.operation == 'M' or op.operation == '=' or op.operation == 'X':
                # Match or mismatch
                q_segment = self.query_sequence[query_pos:query_pos + op.length]
                t_segment = self.target_sequence[target_pos:target_pos + op.length]
                
                query_aligned.append(q_segment)
                target_aligned.append(t_segment)
                
                # Create match indicators
                for q, t in zip(q_segment, t_segment):
                    if q == t:
                        indicator.append('|')  # Match
                    else:
                        indicator.append('.')  # Mismatch
                
                query_pos += op.length
                target_pos += op.length
                
            elif op.operation == 'I':
                # Insertion in query
                query_aligned.append(self.query_sequence[query_pos:query_pos + op.length])
                target_aligned.append('-' * op.length)
                indicator.append(' ' * op.length)
                query_pos += op.length
                
            elif op.operation == 'D':
                # Deletion in query (insertion in target)
                query_aligned.append('-' * op.length)
                target_aligned.append(self.target_sequence[target_pos:target_pos + op.length])
                indicator.append(' ' * op.length)
                target_pos += op.length
        
        self.aligned_query = ''.join(query_aligned)
        self.aligned_target = ''.join(target_aligned)
        self.alignment_indicator = ''.join(indicator)
    
    def get_summary(self) -> str:
        """
        Get a human-readable summary of the alignment.
        
        Returns:
            A string with alignment summary information
        """
        stats = self.statistics
        summary = [
            f"Alignment: {self.query_name} to {self.target_name}",
            f"Query: {self.query_start}-{self.query_end} ({'-' if self.is_reverse else '+'})",
            f"Target: {self.target_start}-{self.target_end}",
            f"Score: {self.score:.1f}",
            f"Identity: {stats.identity:.2%}",
            f"Coverage: {stats.coverage:.2%}",
            f"Matches: {stats.matches}",
            f"Mismatches: {stats.mismatches}",
            f"Gaps: {stats.gaps} (Insertions: {stats.insertions}, Deletions: {stats.deletions})",
            f"Type: {stats.alignment_type.value}"
        ]
        return "\n".join(summary)
