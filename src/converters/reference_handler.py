# src/converters/reference_handler.py
import logging
from typing import Optional
import pyfaidx
import os

logger = logging.getLogger(__name__)

class ReferenceHandlerError(Exception):
    """Custom exception for ReferenceHandler errors."""
    pass

class ReferenceHandler:
    """
    Handles access to reference genome sequences using pyfaidx.

    Provides methods to efficiently fetch subsequences from a FASTA file.
    """

    def __init__(self, fasta_filepath: str):
        """
        Initializes the ReferenceHandler.

        Args:
            fasta_filepath: Path to the reference FASTA file (can be .fa or .fa.gz).

        Raises:
            ReferenceHandlerError: If the FASTA file cannot be opened or indexed.
            FileNotFoundError: If the FASTA file does not exist.
        """
        self.fasta_filepath = fasta_filepath
        self._fasta = None
        if not os.path.exists(fasta_filepath):
             logger.error(f"Reference FASTA file not found: {fasta_filepath}")
             raise FileNotFoundError(f"Reference FASTA file not found: {fasta_filepath}")
        try:
            # pyfaidx automatically handles .fai index creation/checking
            # It might need write permissions in the directory for index creation
            logger.info(f"Attempting to open and index reference FASTA: {fasta_filepath}")
            self._fasta = pyfaidx.Fasta(fasta_filepath, sequence_always_upper=True)
            logger.info(f"Successfully opened and indexed reference FASTA: {fasta_filepath}")
            logger.info(f"Reference sequences available: {list(self._fasta.keys())}")
        except pyfaidx.FastaIndexingError as e:
            logger.error(f"Error indexing FASTA file {fasta_filepath}: {e}. Check file integrity and permissions.")
            raise ReferenceHandlerError(f"Failed to index FASTA file: {e}") from e
        except Exception as e:
            logger.error(f"An unexpected error occurred opening FASTA {fasta_filepath}: {e}")
            raise ReferenceHandlerError(f"Failed to open FASTA file: {e}") from e

    def get_sequence(self, chrom: str, start: int, end: int) -> Optional[str]:
        """
        Retrieves a subsequence from the reference genome.

        Args:
            chrom: The chromosome or contig name.
            start: The 0-based start coordinate (inclusive).
            end: The 0-based end coordinate (exclusive).

        Returns:
            The subsequence as an uppercase string, or None if the chromosome
            does not exist. Returns an empty string if coordinates are invalid
            or outside chromosome bounds.

        Raises:
            ReferenceHandlerError: If there's an issue accessing the sequence data.
        """
        if self._fasta is None:
            raise ReferenceHandlerError("FASTA reference is not loaded.")
        if start < 0 or end < start:
            logger.warning(f"Invalid coordinates requested: {chrom}:{start}-{end}")
            return "" # Return empty string for invalid coordinates
        if chrom not in self._fasta:
            logger.warning(f"Chromosome '{chrom}' not found in reference FASTA.")
            return None

        try:
            # pyfaidx uses 0-based coordinates, matching Python slicing
            # Ensure coordinates are within bounds
            chrom_len = len(self._fasta[chrom])
            if start >= chrom_len:
                logger.warning(f"Start coordinate {start} is beyond chromosome {chrom} length {chrom_len}.")
                return "" # Return empty string if start is out of bounds
            # Adjust end if it exceeds chromosome length
            safe_end = min(end, chrom_len)
            if start >= safe_end: # Handle cases where start >= end after adjustment
                 logger.warning(f"Start coordinate {start} >= end coordinate {safe_end} for {chrom}. Returning empty string.")
                 return ""

            # Fetch sequence using pyfaidx FastaRecord object slicing
            sequence = self._fasta[chrom][start:safe_end].seq
            # logger.debug(f"Fetched sequence for {chrom}:{start}-{end} (len {len(sequence)})")
            return sequence
        except KeyError:
            # This should be caught by the 'chrom not in self._fasta' check, but added for safety
            logger.warning(f"Chromosome '{chrom}' not found during sequence retrieval.")
            return None
        except Exception as e:
            logger.error(f"Error retrieving sequence for {chrom}:{start}-{end}: {e}")
            raise ReferenceHandlerError(f"Failed to retrieve sequence: {e}") from e

    def get_chrom_length(self, chrom: str) -> Optional[int]:
        """
        Gets the length of a specific chromosome or contig.

        Args:
            chrom: The chromosome or contig name.

        Returns:
            The length of the chromosome, or None if not found.
        """
        if self._fasta is None:
            raise ReferenceHandlerError("FASTA reference is not loaded.")
        if chrom not in self._fasta:
            logger.warning(f"Chromosome '{chrom}' not found when requesting length.")
            return None
        try:
            return len(self._fasta[chrom])
        except Exception as e:
            logger.error(f"Error getting length for chromosome {chrom}: {e}")
            raise ReferenceHandlerError(f"Failed to get chromosome length: {e}") from e

    def close(self):
        """Closes the FASTA file handle."""
        if self._fasta:
            try:
                self._fasta.close()
                logger.info(f"Closed reference FASTA file: {self.fasta_filepath}")
            except Exception as e:
                logger.error(f"Error closing FASTA file {self.fasta_filepath}: {e}")
        self._fasta = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def list_chromosomes(self) -> list[str]:
        """Returns a list of chromosome/contig names in the reference."""
        if self._fasta is None:
            raise ReferenceHandlerError("FASTA reference is not loaded.")
        return list(self._fasta.keys())
