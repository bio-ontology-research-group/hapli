"""
Parser for FASTA sequence files.
Uses Biopython's SeqIO module to handle FASTA format parsing.
"""
import os
import logging
from typing import Dict, List, Optional, Union, Iterator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class FastaParser:
    """
    Parser for FASTA sequence files.
    Provides methods to access sequences by ID and retrieve sequence information.
    """
    
    def __init__(self):
        """Initialize the FASTA parser."""
        self.logger = logging.getLogger(__name__)
        self.sequences = {}  # Dictionary of SeqRecord objects indexed by ID
        
    def parse(self, filepath: str) -> Dict[str, SeqRecord]:
        """
        Parse a FASTA file and return a dictionary of sequences.
        
        Args:
            filepath: Path to the FASTA file
            
        Returns:
            Dictionary of SeqRecord objects indexed by sequence ID
            
        Raises:
            FileNotFoundError: If the file doesn't exist
            ValueError: If the file is empty or malformed
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"FASTA file not found: {filepath}")
            
        self.logger.info(f"Parsing FASTA file: {filepath}")
        
        # Check if file is empty
        if os.path.getsize(filepath) == 0:
            raise ValueError("Empty FASTA file")
            
        try:
            # Use a context manager to properly close the file handle
            with open(filepath, 'r') as fasta_file:
                # Check for valid FASTA format (should start with '>')
                first_line = fasta_file.readline().strip()
                if not first_line.startswith('>'):
                    raise ValueError("File does not appear to be in FASTA format")
                
                # Reset file pointer to beginning
                fasta_file.seek(0)
                
                # Parse the file
                records = SeqIO.parse(fasta_file, "fasta")
                self.sequences = SeqIO.to_dict(records)
            
            if not self.sequences:
                raise ValueError("No sequences found in FASTA file")
                
            self.logger.info(f"Successfully parsed FASTA with {len(self.sequences)} sequences")
            return self.sequences
            
        except ValueError as e:
            self.logger.error(f"Failed to parse FASTA file: {e}")
            self.sequences = {}  # Reset sequences to empty dict
            raise
        except Exception as e:
            self.logger.error(f"Failed to parse FASTA file: {e}")
            self.sequences = {}  # Reset sequences to empty dict
            raise
    
    def get_sequence(self, seq_id: str) -> Optional[SeqRecord]:
        """
        Get a sequence by its ID.
        
        Args:
            seq_id: ID of the sequence to retrieve
            
        Returns:
            SeqRecord object for the sequence, or None if not found
        """
        return self.sequences.get(seq_id)
    
    def get_sequence_string(self, seq_id: str) -> Optional[str]:
        """
        Get a sequence as a string by its ID.
        
        Args:
            seq_id: ID of the sequence to retrieve
            
        Returns:
            String representation of the sequence, or None if not found
        """
        seq_record = self.get_sequence(seq_id)
        if seq_record:
            return str(seq_record.seq)
        return None
    
    def get_all_sequences(self) -> Dict[str, SeqRecord]:
        """
        Get all sequences.
        
        Returns:
            Dictionary of all sequences indexed by ID
        """
        return self.sequences
    
    def get_sequence_length(self, seq_id: str) -> Optional[int]:
        """
        Get the length of a sequence by its ID.
        
        Args:
            seq_id: ID of the sequence
            
        Returns:
            Length of the sequence, or None if the sequence is not found
        """
        seq_record = self.get_sequence(seq_id)
        if seq_record:
            return len(seq_record.seq)
        return None
