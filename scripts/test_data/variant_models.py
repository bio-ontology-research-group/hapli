"""
Data classes for representing different types of genetic variants.
"""

from dataclasses import dataclass
from typing import Optional, Tuple
import uuid


@dataclass
class Variant:
    """Base class for all genetic variants."""
    chromosome: str
    position: int  # 1-based position
    ref_allele: str
    alt_allele: str
    variant_type: str
    variant_id: str = ""
    
    def __post_init__(self):
        """Validate base variant fields."""
        if not self.variant_id:
            self.variant_id = str(uuid.uuid4())
        
        if self.position < 1:
            raise ValueError("Position must be 1-based (>= 1)")
        
        if not self.chromosome:
            raise ValueError("Chromosome cannot be empty")
        
        if not self.ref_allele:
            raise ValueError("Reference allele cannot be empty")
    
    def to_english_description(self) -> str:
        """Return a human-readable description of the variant."""
        return f"{self.variant_type} at {self.chromosome}:{self.position} ({self.ref_allele} -> {self.alt_allele})"
    
    def to_hgvs_notation(self) -> str:
        """Return HGVS notation for the variant (placeholder implementation)."""
        # Placeholder - would need proper HGVS formatting logic
        return f"{self.chromosome}:g.{self.position}{self.ref_allele}>{self.alt_allele}"
    
    def get_affected_region(self) -> Tuple[int, int]:
        """Return the genomic region affected by this variant as (start, end)."""
        return (self.position, self.position + len(self.ref_allele) - 1)


@dataclass
class SNV(Variant):
    """Single nucleotide variant."""
    
    def __post_init__(self):
        super().__post_init__()
        self.variant_type = "SNV"
        
        if len(self.ref_allele) != 1:
            raise ValueError("SNV reference allele must be single nucleotide")
        
        if len(self.alt_allele) != 1:
            raise ValueError("SNV alternate allele must be single nucleotide")
        
        if self.ref_allele == self.alt_allele:
            raise ValueError("Reference and alternate alleles cannot be the same")
    
    def to_english_description(self) -> str:
        return f"Single nucleotide variant at {self.chromosome}:{self.position} ({self.ref_allele} -> {self.alt_allele})"


@dataclass
class Insertion(Variant):
    """Insertion variant."""
    insert_sequence: str = ""
    
    def __post_init__(self):
        super().__post_init__()
        self.variant_type = "Insertion"
        
        if not self.insert_sequence:
            self.insert_sequence = self.alt_allele
        
        if not self.insert_sequence:
            raise ValueError("Insertion sequence cannot be empty")
    
    def to_english_description(self) -> str:
        return f"Insertion of {len(self.insert_sequence)} bp at {self.chromosome}:{self.position}"
    
    def get_affected_region(self) -> Tuple[int, int]:
        # Insertion affects the position between ref bases
        return (self.position, self.position)


@dataclass
class Deletion(Variant):
    """Deletion variant."""
    deleted_length: int = 0
    
    def __post_init__(self):
        super().__post_init__()
        self.variant_type = "Deletion"
        
        if self.deleted_length == 0:
            self.deleted_length = len(self.ref_allele)
        
        if self.deleted_length <= 0:
            raise ValueError("Deleted length must be positive")
    
    def to_english_description(self) -> str:
        return f"Deletion of {self.deleted_length} bp at {self.chromosome}:{self.position}"
    
    def get_affected_region(self) -> Tuple[int, int]:
        return (self.position, self.position + self.deleted_length - 1)


@dataclass
class Inversion(Variant):
    """Inversion variant."""
    end_position: int = 0
    
    def __post_init__(self):
        super().__post_init__()
        self.variant_type = "Inversion"
        
        if self.end_position == 0:
            raise ValueError("End position must be specified for inversion")
        
        if self.end_position <= self.position:
            raise ValueError("End position must be greater than start position")
    
    def to_english_description(self) -> str:
        length = self.end_position - self.position + 1
        return f"Inversion of {length} bp at {self.chromosome}:{self.position}-{self.end_position}"
    
    def get_affected_region(self) -> Tuple[int, int]:
        return (self.position, self.end_position)


@dataclass
class Duplication(Variant):
    """Duplication variant."""
    end_position: int = 0
    copy_number: int = 2
    
    def __post_init__(self):
        super().__post_init__()
        self.variant_type = "Duplication"
        
        if self.end_position == 0:
            raise ValueError("End position must be specified for duplication")
        
        if self.end_position <= self.position:
            raise ValueError("End position must be greater than start position")
        
        if self.copy_number < 2:
            raise ValueError("Copy number for duplication must be >= 2")
    
    def to_english_description(self) -> str:
        length = self.end_position - self.position + 1
        return f"Duplication of {length} bp at {self.chromosome}:{self.position}-{self.end_position} (copy number: {self.copy_number})"
    
    def get_affected_region(self) -> Tuple[int, int]:
        return (self.position, self.end_position)


@dataclass
class Translocation(Variant):
    """Translocation variant."""
    chr2: str = ""
    pos2: int = 0
    strand1: str = "+"
    strand2: str = "+"
    
    def __post_init__(self):
        super().__post_init__()
        self.variant_type = "Translocation"
        
        if not self.chr2:
            raise ValueError("Second chromosome must be specified for translocation")
        
        if self.pos2 < 1:
            raise ValueError("Second position must be 1-based (>= 1)")
        
        if self.strand1 not in ["+", "-"]:
            raise ValueError("Strand1 must be '+' or '-'")
        
        if self.strand2 not in ["+", "-"]:
            raise ValueError("Strand2 must be '+' or '-'")
    
    def to_english_description(self) -> str:
        return f"Translocation between {self.chromosome}:{self.position}({self.strand1}) and {self.chr2}:{self.pos2}({self.strand2})"
    
    def get_affected_region(self) -> Tuple[int, int]:
        # For translocation, return the breakpoint on the first chromosome
        return (self.position, self.position)
