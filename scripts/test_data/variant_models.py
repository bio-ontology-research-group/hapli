"""
Data classes for representing different types of genetic variants.
"""

from dataclasses import dataclass
from typing import Optional, Tuple, Dict
import uuid
import re


class HGVSFormatter:
    """
    Generates HGVS-compliant notation for test data variants.
    
    This is a simplified implementation suitable for test data generation.
    For production use, consider biocommons.hgvs with full sequence validation.
    Follows HGVS nomenclature v15.11 for genomic reference sequences.
    """
    
    def __init__(self, reference_fasta_path=None):
        """
        Initialize HGVSFormatter.
        
        Args:
            reference_fasta_path: Path to reference FASTA file (optional)
        """
        self.reference = None
        if reference_fasta_path:
            try:
                from pyfaidx import Fasta
                self.reference = Fasta(reference_fasta_path)
            except ImportError:
                print("Warning: pyfaidx not available. Reference sequence validation disabled.")
        
        # Standard chromosome to RefSeq accession mapping (GRCh38/hg38)
        self.chr_to_accession = {
            'chr1': 'NC_000001.11', 'chr2': 'NC_000002.12', 
            'chr3': 'NC_000003.12', 'chr4': 'NC_000004.12',
            'chr5': 'NC_000005.10', 'chr6': 'NC_000006.12',
            'chr7': 'NC_000007.14', 'chr8': 'NC_000008.11',
            'chr9': 'NC_000009.12', 'chr10': 'NC_000010.11',
            'chr11': 'NC_000011.10', 'chr12': 'NC_000012.12',
            'chr13': 'NC_000013.11', 'chr14': 'NC_000014.9',
            'chr15': 'NC_000015.10', 'chr16': 'NC_000016.10',
            'chr17': 'NC_000017.11', 'chr18': 'NC_000018.10',
            'chr19': 'NC_000019.10', 'chr20': 'NC_000020.11',
            'chr21': 'NC_000021.9', 'chr22': 'NC_000022.11',
            'chrX': 'NC_000023.11', 'chrY': 'NC_000024.10',
            'chrM': 'NC_012920.1',
            # Support both with and without 'chr' prefix
            '1': 'NC_000001.11', '2': 'NC_000002.12', 
            '3': 'NC_000003.12', '4': 'NC_000004.12',
            '5': 'NC_000005.10', '6': 'NC_000006.12',
            '7': 'NC_000007.14', '8': 'NC_000008.11',
            '9': 'NC_000009.12', '10': 'NC_000010.11',
            '11': 'NC_000011.10', '12': 'NC_000012.12',
            '13': 'NC_000013.11', '14': 'NC_000014.9',
            '15': 'NC_000015.10', '16': 'NC_000016.10',
            '17': 'NC_000017.11', '18': 'NC_000018.10',
            '19': 'NC_000019.10', '20': 'NC_000020.11',
            '21': 'NC_000021.9', '22': 'NC_000022.11',
            'X': 'NC_000023.11', 'Y': 'NC_000024.10',
            'M': 'NC_012920.1', 'MT': 'NC_012920.1'
        }
    
    def get_accession(self, chromosome: str) -> str:
        """Get RefSeq accession for chromosome."""
        # Normalize chromosome name
        chr_name = chromosome.lower()
        if chr_name.startswith('chr'):
            chr_key = chr_name
        else:
            chr_key = f'chr{chr_name}'
        
        return self.chr_to_accession.get(chr_key, chromosome)
    
    def format_genomic_substitution(self, chr: str, pos: int, ref: str, alt: str) -> str:
        """Format SNV/substitution: NC_000001.11:g.12345A>G"""
        acc = self.get_accession(chr)
        return f"{acc}:g.{pos}{ref}>{alt}"
    
    def format_genomic_deletion(self, chr: str, start: int, end: int, deleted_seq: Optional[str] = None) -> str:
        """Format deletion: NC_000001.11:g.12345_12350del or g.12345del"""
        acc = self.get_accession(chr)
        if start == end:
            return f"{acc}:g.{start}del"
        return f"{acc}:g.{start}_{end}del"
    
    def format_genomic_insertion(self, chr: str, pos_before: int, pos_after: int, inserted_seq: str) -> str:
        """Format insertion: NC_000001.11:g.12345_12346insATCG"""
        acc = self.get_accession(chr)
        return f"{acc}:g.{pos_before}_{pos_after}ins{inserted_seq}"
    
    def format_genomic_duplication(self, chr: str, start: int, end: int) -> str:
        """Format duplication: NC_000001.11:g.12345_12350dup"""
        acc = self.get_accession(chr)
        if start == end:
            return f"{acc}:g.{start}dup"
        return f"{acc}:g.{start}_{end}dup"
    
    def format_genomic_inversion(self, chr: str, start: int, end: int) -> str:
        """Format inversion: NC_000001.11:g.12345_12350inv"""
        acc = self.get_accession(chr)
        return f"{acc}:g.{start}_{end}inv"
    
    def format_genomic_delins(self, chr: str, start: int, end: int, inserted_seq: str) -> str:
        """Format deletion-insertion: NC_000001.11:g.12345_12350delinsATCG"""
        acc = self.get_accession(chr)
        if start == end:
            return f"{acc}:g.{start}delins{inserted_seq}"
        return f"{acc}:g.{start}_{end}delins{inserted_seq}"
    
    def format_genomic_copy_number(self, chr: str, start: int, end: int, copy_number: int) -> str:
        """Format copy number variant: NC_000001.11:g.12345_12350[5]"""
        acc = self.get_accession(chr)
        return f"{acc}:g.{start}_{end}[{copy_number}]"
    
    def validate_hgvs_syntax(self, hgvs_string: str) -> bool:
        """Basic syntax validation using regex"""
        # Pattern for genomic HGVS notation
        patterns = [
            r'^NC_\d+\.\d+:g\.\d+[ATCGN]+>[ATCGN]+$',  # substitution
            r'^NC_\d+\.\d+:g\.\d+(_\d+)?del$',  # deletion
            r'^NC_\d+\.\d+:g\.\d+_\d+ins[ATCGN]+$',  # insertion
            r'^NC_\d+\.\d+:g\.\d+(_\d+)?dup$',  # duplication
            r'^NC_\d+\.\d+:g\.\d+_\d+inv$',  # inversion
            r'^NC_\d+\.\d+:g\.\d+(_\d+)?delins[ATCGN]+$',  # deletion-insertion
            r'^NC_\d+\.\d+:g\.\d+_\d+\[\d+\]$',  # copy number
        ]
        
        return any(re.match(pattern, hgvs_string) for pattern in patterns)
    
    def get_reference_sequence(self, chr: str, start: int, end: int) -> Optional[str]:
        """Get reference sequence for validation (if reference provided)"""
        if self.reference and chr in self.reference:
            # Convert to 0-based for pyfaidx
            return str(self.reference[chr][start-1:end])
        return None


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
    
    def to_hgvs_notation(self, formatter: Optional[HGVSFormatter] = None) -> str:
        """Return HGVS notation for the variant."""
        if not formatter:
            formatter = HGVSFormatter()
        return formatter.format_genomic_substitution(
            self.chromosome, self.position, self.ref_allele, self.alt_allele
        )
    
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
    
    def to_hgvs_notation(self, formatter: Optional[HGVSFormatter] = None) -> str:
        """Return HGVS notation for SNV."""
        if not formatter:
            formatter = HGVSFormatter()
        return formatter.format_genomic_substitution(
            self.chromosome, self.position, self.ref_allele, self.alt_allele
        )


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
    
    def to_hgvs_notation(self, formatter: Optional[HGVSFormatter] = None) -> str:
        """Return HGVS notation for insertion."""
        if not formatter:
            formatter = HGVSFormatter()
        # For insertions, HGVS uses the positions flanking the insertion
        return formatter.format_genomic_insertion(
            self.chromosome, self.position, self.position + 1, self.insert_sequence
        )
    
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
    
    def to_hgvs_notation(self, formatter: Optional[HGVSFormatter] = None) -> str:
        """Return HGVS notation for deletion."""
        if not formatter:
            formatter = HGVSFormatter()
        end_position = self.position + self.deleted_length - 1
        return formatter.format_genomic_deletion(
            self.chromosome, self.position, end_position
        )
    
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
    
    def to_hgvs_notation(self, formatter: Optional[HGVSFormatter] = None) -> str:
        """Return HGVS notation for inversion."""
        if not formatter:
            formatter = HGVSFormatter()
        return formatter.format_genomic_inversion(
            self.chromosome, self.position, self.end_position
        )
    
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
    
    def to_hgvs_notation(self, formatter: Optional[HGVSFormatter] = None) -> str:
        """Return HGVS notation for duplication."""
        if not formatter:
            formatter = HGVSFormatter()
        
        # For simple duplications (copy number = 2), use standard dup notation
        if self.copy_number == 2:
            return formatter.format_genomic_duplication(
                self.chromosome, self.position, self.end_position
            )
        else:
            # For higher copy numbers, use copy number notation
            return formatter.format_genomic_copy_number(
                self.chromosome, self.position, self.end_position, self.copy_number
            )
    
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
    
    def to_hgvs_notation(self, formatter: Optional[HGVSFormatter] = None) -> str:
        """Return HGVS notation for translocation."""
        if not formatter:
            formatter = HGVSFormatter()
        
        # Translocations are complex and don't have a standard simple HGVS format
        # Return descriptive text with note about complexity
        acc1 = formatter.get_accession(self.chromosome)
        acc2 = formatter.get_accession(self.chr2)
        return f"Complex translocation: {acc1}:g.{self.position} and {acc2}:g.{self.pos2} (see structural variant description)"
    
    def get_affected_region(self) -> Tuple[int, int]:
        # For translocation, return the breakpoint on the first chromosome
        return (self.position, self.position)


@dataclass
class ComplexVariant(Variant):
    """Complex variant with deletion and insertion."""
    deleted_length: int = 0
    insert_sequence: str = ""
    
    def __post_init__(self):
        super().__post_init__()
        self.variant_type = "Complex"
        
        if self.deleted_length == 0:
            self.deleted_length = len(self.ref_allele)
        
        if not self.insert_sequence:
            self.insert_sequence = self.alt_allele
        
        if self.deleted_length <= 0:
            raise ValueError("Deleted length must be positive")
        
        if not self.insert_sequence:
            raise ValueError("Insert sequence cannot be empty")
    
    def to_english_description(self) -> str:
        return f"Complex variant: deletion of {self.deleted_length} bp and insertion of {len(self.insert_sequence)} bp at {self.chromosome}:{self.position}"
    
    def to_hgvs_notation(self, formatter: Optional[HGVSFormatter] = None) -> str:
        """Return HGVS notation for complex variant (deletion-insertion)."""
        if not formatter:
            formatter = HGVSFormatter()
        
        end_position = self.position + self.deleted_length - 1
        return formatter.format_genomic_delins(
            self.chromosome, self.position, end_position, self.insert_sequence
        )
    
    def get_affected_region(self) -> Tuple[int, int]:
        return (self.position, self.position + self.deleted_length - 1)
