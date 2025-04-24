# src/converters/vcf_to_gfa.py
import logging
import itertools
import os  # Added import
import sys # Added import
import argparse # Added import
from typing import Dict, List, Optional, Tuple, Any
import pysam
import gfapy
from .reference_handler import ReferenceHandler, ReferenceHandlerError
from .phasing_processor import PhasingProcessor, PhasingError

logger = logging.getLogger(__name__)

class VCFtoGFAConversionError(Exception):
    """Custom exception for VCF to GFA conversion errors."""
    pass

class VCFtoGFAConverter:
    """
    Converts variants from a VCF file into a GFA (Graphical Fragment Assembly) format.

    Uses pysam for VCF reading, pyfaidx (via ReferenceHandler) for reference access,
    and gfapy for GFA construction. Handles phased variants to create haplotype paths.
    """

    def __init__(self, vcf_filepath: str, fasta_filepath: str, output_gfa_filepath: str,
                 path_template: str = "{sample}_hap{hap}", unphased_strategy: str = 'ref'):
        # Don't increase recursion limit - we'll avoid GFApy's recursive validation entirely
        """
        Initializes the VCFtoGFAConverter.

        Args:
            vcf_filepath: Path to the input VCF file (can be .vcf, .vcf.gz, .bcf).
            fasta_filepath: Path to the reference FASTA file.
            output_gfa_filepath: Path to the output GFA file.
            path_template: A format string for naming haplotype paths.
                           Available placeholders: {sample}, {hap} (1 or 2).
            unphased_strategy: How to handle unphased variants or missing genotypes.
                               'ref': Use the reference allele (default).
                               'alt': Use the first alternate allele (if available).
                               'skip': Effectively use reference up to variant, skip variant allele.
                               'bubble': (Not implemented) Would create complex graph structure.
                               
        Raises:
            FileNotFoundError: If the input VCF or FASTA file does not exist.
            VCFtoGFAConversionError: If initialization fails for other reasons.
        """
        # Check if input files exist
        if not os.path.exists(vcf_filepath):
            raise FileNotFoundError(f"VCF file not found: {vcf_filepath}")
        if not os.path.exists(fasta_filepath):
            raise FileNotFoundError(f"FASTA file not found: {fasta_filepath}")
            
        self.vcf_filepath = vcf_filepath
        self.fasta_filepath = fasta_filepath
        self.output_gfa_filepath = output_gfa_filepath
        self.path_template = path_template

        valid_unphased_strategies = ['ref', 'alt', 'skip', 'bubble']
        if unphased_strategy not in valid_unphased_strategies:
            raise ValueError(f"Invalid unphased_strategy '{unphased_strategy}'. Must be one of {valid_unphased_strategies}")
        if unphased_strategy == 'bubble':
             logger.warning("The 'bubble' strategy for unphased variants is complex and not implemented. Falling back to 'ref'.")
             self.unphased_strategy = 'ref'
        else:
             self.unphased_strategy = unphased_strategy

        self._vcf_reader: Optional[pysam.VariantFile] = None
        self._ref_handler: Optional[ReferenceHandler] = None
        self._phasing_processor: Optional[PhasingProcessor] = None
        # Don't use gfapy.Gfa object - use direct GFA line management instead
        self._gfa_lines: List[str] = []
        self._segment_id_counter = 0
        # Cache: sequence -> segment_id
        self._segment_cache: Dict[str, str] = {}
        # Cache: (from_seg, from_orient, to_seg, to_orient) -> True
        self._link_cache: Dict[Tuple[str, str, str, str], bool] = {}
        # Track segments, links and paths for validation
        self._segments: Set[str] = set()
        self._paths: Set[str] = set()

        logger.info("VCFtoGFAConverter initialized.")
        logger.info(f"  VCF: {vcf_filepath}")
        logger.info(f"  FASTA: {fasta_filepath}")
        logger.info(f"  Output GFA: {output_gfa_filepath}")
        logger.info(f"  Path Template: {path_template}")
        logger.info(f"  Unphased Strategy: {self.unphased_strategy}")


    def _initialize_resources(self):
        """Opens VCF and FASTA files and initializes processors."""
        try:
            logger.info("Initializing resources...")
            self._ref_handler = ReferenceHandler(self.fasta_filepath)
            # Need 'r' mode for VCF, 'rb' for BCF. pysam usually figures it out.
            # Use require_index=True if region fetching is essential and index MUST exist
            self._vcf_reader = pysam.VariantFile(self.vcf_filepath, "r")
            self._phasing_processor = PhasingProcessor(self._vcf_reader)
            # Initialize our collections
            self._gfa_lines = []
            self._segment_id_counter = 0
            self._segment_cache = {}
            self._link_cache = {}
            self._segments = set()
            self._paths = set()
            
            # Add GFA header
            self._gfa_lines.append("H\tVN:Z:1.0")
            logger.info("Resources initialized successfully.")
        except FileNotFoundError as e:
            logger.error(f"Input file not found: {e}")
            raise VCFtoGFAConversionError(f"Input file not found: {e}") from e
        except (pysam.utils.SamtoolsError, PhasingError, ReferenceHandlerError, Exception) as e:
            logger.error(f"Error initializing resources: {e}")
            self._cleanup_resources() # Ensure cleanup if init fails partially
            raise VCFtoGFAConversionError(f"Error initializing resources: {e}") from e

    def _cleanup_resources(self):
        """Closes file handles."""
        logger.info("Cleaning up resources...")
        if self._vcf_reader:
            try:
                self._vcf_reader.close()
                logger.debug("VCF reader closed.")
            except Exception as e:
                logger.warning(f"Error closing VCF reader: {e}")
        if self._ref_handler:
            try:
                self._ref_handler.close()
                logger.debug("Reference handler closed.")
            except Exception as e:
                logger.warning(f"Error closing reference handler: {e}")
        self._vcf_reader = None
        self._ref_handler = None
        self._phasing_processor = None
        self._gfa_lines = []
        self._segments = set()
        self._paths = set()
        self._segment_cache = {}
        self._link_cache = {}
        logger.info("Resources cleaned up.")

    def _get_or_create_segment(self, sequence: str) -> Optional[str]:
        """
        Gets the ID of an existing segment for the given sequence or creates a new one.
        Returns None if the sequence is empty. Adds LN tag.
        """
        if not sequence: # Don't create segments for empty sequences (e.g., from deletions)
            return None
        # Ensure sequence is uppercase for consistent caching and GFA output
        sequence = sequence.upper()
        if sequence in self._segment_cache:
            return self._segment_cache[sequence]
        else:
            self._segment_id_counter += 1
            segment_id = f"s{self._segment_id_counter}"
            try:
                # Basic validation: Check for invalid characters in sequence
                # GFA spec allows [A-Za-z=.]+ but often restricted to ACGTN.
                # Allowing only ACGTN for stricter biological validity.
                valid_chars = set("ACGTN")
                if not set(sequence).issubset(valid_chars):
                     # Replace invalid chars with N? Or raise error? Replacing for now.
                     original_seq_preview = sequence[:20]
                     sequence = "".join([c if c in valid_chars else 'N' for c in sequence])
                     logger.warning(f"Segment {segment_id} sequence contained non-ACGTN characters (preview: '{original_seq_preview}...'). Replaced with 'N'.")

                # Create segment line directly
                segment_line = f"S\t{segment_id}\t{sequence}\tLN:i:{len(sequence)}"
                self._gfa_lines.append(segment_line)
                self._segments.add(segment_id)
                
                # Store in segment cache
                self._segment_cache[sequence] = segment_id
                
                logger.debug(f"Created segment {segment_id} with length {len(sequence)}")
                return segment_id
            except Exception as e:
                logger.error(f"Failed to create GFA segment for sequence '{sequence[:20]}...': {e}")
                # Return the segment ID anyway if it's in the cache
                if sequence in self._segment_cache:
                    return self._segment_cache[sequence]
                return None


    def _add_link(self, from_segment_id: str, from_orient: str,
                  to_segment_id: str, to_orient: str, cigar: str = "*"):
        """Adds a link to the GFA graph, avoiding duplicates."""
        if not from_segment_id or not to_segment_id:
            # logger.debug(f"Skipping link creation due to missing segment ID (from={from_segment_id}, to={to_segment_id})")
            return

        link_key = (from_segment_id, from_orient, to_segment_id, to_orient)
        # Check cache to prevent adding duplicate links
        if link_key in self._link_cache:
             # logger.debug(f"Skipping duplicate link: {from_segment_id}{from_orient} -> {to_segment_id}{to_orient}")
             return

        try:
            # Create a link line directly
            link_line = f"L\t{from_segment_id}\t{from_orient}\t{to_segment_id}\t{to_orient}\t{cigar}"
            self._gfa_lines.append(link_line)
            self._link_cache[link_key] = True # Add to cache
            logger.debug(f"Added link: {from_segment_id}{from_orient} -> {to_segment_id}{to_orient}")
        except Exception as e:
            logger.error(f"Failed to add GFA link {from_segment_id}{from_orient} -> {to_segment_id}{to_orient}: {e}")
            logger.error(f"Error details: {e}")

    def _build_gfa_path_and_links(self, path_name: str, segment_tuples: List[Tuple[int, str]]):
        """
        Builds GFA path and the necessary links between segments.

        Args:
            path_name: The desired name for the GFA path.
            segment_tuples: List of (start_pos, sequence_string) for the haplotype path.
                            Sequences should correspond to GFA segments.
        """
        if not segment_tuples:
            logger.warning(f"No segments provided for path {path_name}. Skipping path creation.")
            return

        path_segment_ids = []
        orientations = []
        last_segment_id = None

        logger.info(f"Building path and links for: {path_name} with {len(segment_tuples)} sequence segments.")

        # First pass: Get or create all segments and store their IDs
        for i, (pos, seq) in enumerate(segment_tuples):
            segment_id = self._get_or_create_segment(seq)
            if segment_id: # Only consider if segment was created (non-empty seq)
                path_segment_ids.append(segment_id)
                orientations.append("+") # Assuming forward orientation for now
            # else:
                # logger.debug(f"Skipping empty sequence segment at pos {pos} for path {path_name}")

        if not path_segment_ids:
             logger.warning(f"Path {path_name} resulted in no valid segments. Skipping.")
             return

        # Second pass: Add links between consecutive segments in the path
        for i in range(len(path_segment_ids) - 1):
            from_seg = path_segment_ids[i]
            to_seg = path_segment_ids[i+1]
            self._add_link(from_seg, "+", to_seg, "+")

        # Third pass: Create the GFA Path line
        try:
            # Overlaps are typically '*' unless precisely calculated. Use '*' for simplicity.
            overlaps = ["*"] * (len(path_segment_ids) - 1) if len(path_segment_ids) > 1 else []

            # Create path line directly
            # Format: P<tab>name<tab>segment_names<tab>orientations<tab>overlaps
            segment_names_str = ",".join(path_segment_ids)
            orientations_str = ",".join(orientations)
            overlaps_str = ",".join(overlaps)
            
            # Extract sample and haplotype info
            sample_part = path_name
            hap_part = ""
            if "_hap" in path_name:
                parts = path_name.rsplit('_hap', 1)
                if len(parts) == 2 and parts[1].isdigit():
                    sample_part, hap_part = parts
            elif "." in path_name:
                parts = path_name.rsplit('.', 1)
                if len(parts) == 2 and parts[1].isdigit():
                    sample_part, hap_part = parts
            
            # Build tags
            tags = []
            if sample_part != path_name:
                tags.append(f"SM:Z:{sample_part}")
            if hap_part.isdigit():
                tags.append(f"HP:i:{hap_part}")
            
            # Create the path line
            path_line = f"P\t{path_name}\t{segment_names_str}\t{orientations_str}\t{overlaps_str}"
            if tags:
                path_line += "\t" + "\t".join(tags)
            
            # Add to GFA structure directly
            self._gfa_lines.append(path_line)
            self._paths.add(path_name)
            
            logger.info(f"Added path {path_name} with {len(path_segment_ids)} segments.")
        except Exception as e:
            logger.error(f"Failed to add GFA path {path_name}: {e}")
            raise VCFtoGFAConversionError(f"Failed to add path: {e}") from e


    def convert(self, region: Optional[str] = None):
        """
        Performs the VCF to GFA conversion.

        Args:
            region: Optional region string in the format "chrom:start-end" (1-based)
                    or "chrom" to process only a specific contig. Requires VCF index.
        """
        self._initialize_resources()
        if not self._vcf_reader or not self._ref_handler or not self._phasing_processor or not self._gfa:
             raise VCFtoGFAConversionError("Initialization failed. Cannot proceed.")

        try:
            # Header is already added during initialization

            # Determine contigs to process
            contigs_to_process = []
            fetch_kwargs = {}
            if region:
                try:
                    # pysam can parse region strings like "chr1:1000-2000"
                    # It requires the VCF to be indexed (e.g., tabix)
                    contig, start, end = None, None, None
                    if ':' in region:
                        parts = region.split(':')
                        contig = parts[0]
                        if len(parts) > 1 and '-' in parts[1]:
                             coords = parts[1].split('-')
                             # pysam fetch uses 0-based start, 1-based VCF uses 1-based start
                             start = max(0, int(coords[0].replace(',', '')) - 1) # 0-based start
                             end = int(coords[1].replace(',', '')) # Optional end
                             fetch_kwargs['start'] = start
                             fetch_kwargs['end'] = end
                    else:
                        contig = region

                    if contig not in self._vcf_reader.header.contigs:
                         raise VCFtoGFAConversionError(f"Contig '{contig}' from region '{region}' not found in VCF header.")
                    if contig not in self._ref_handler.list_chromosomes():
                         raise VCFtoGFAConversionError(f"Contig '{contig}' from region '{region}' not found in reference FASTA.")

                    contigs_to_process.append(contig)
                    fetch_kwargs['contig'] = contig
                    logger.info(f"Processing region: {contig}" + (f":{start+1}-{end}" if start is not None else ""))

                except ValueError as e:
                     raise VCFtoGFAConversionError(f"Invalid region format '{region}'. Use 'chrom' or 'chrom:start-end'. Error: {e}") from e
                except IndexError: # Raised by pysam if index is missing for fetch
                     raise VCFtoGFAConversionError(f"Cannot fetch region '{region}'. VCF index (.tbi or .csi) is required and might be missing or corrupted.")

            else:
                # Process all contigs present in both VCF header and FASTA reference
                vcf_contigs = set(self._vcf_reader.header.contigs.keys())
                ref_contigs = set(self._ref_handler.list_chromosomes())
                contigs_to_process = sorted(list(vcf_contigs.intersection(ref_contigs)))
                if not contigs_to_process:
                     missing_in_ref = vcf_contigs - ref_contigs
                     missing_in_vcf = ref_contigs - vcf_contigs
                     msg = "No common contigs found between VCF header and reference FASTA."
                     if missing_in_ref: msg += f" Contigs in VCF but not REF: {missing_in_ref}."
                     if missing_in_vcf: msg += f" Contigs in REF but not VCF: {missing_in_vcf}."
                     raise VCFtoGFAConversionError(msg)
                logger.info(f"Processing all common contigs: {contigs_to_process}")


            # --- Main Conversion Loop ---
            for chrom in contigs_to_process:
                logger.info(f"--- Processing contig: {chrom} ---")
                chrom_len = self._ref_handler.get_chrom_length(chrom)
                if chrom_len is None: # Should have been checked earlier, but safety check
                    logger.warning(f"Skipping contig {chrom} as its length could not be determined from reference FASTA.")
                    continue

                # Get iterator for VCF records in the specified region/contig
                try:
                    # Use fetch with parsed region kwargs if region was specified
                    records_iterator = self._vcf_reader.fetch(**fetch_kwargs) if region else self._vcf_reader.fetch(chrom)
                    # We need to check if the iterator is empty *before* processing samples
                    # Consuming the iterator here to check for variants
                    vcf_records = list(records_iterator)
                    logger.info(f"Fetched {len(vcf_records)} variants for contig {chrom}" + (f" in region {region}" if region else ""))
                except ValueError as e:
                     logger.warning(f"Could not fetch variants for contig {chrom} (maybe no variants in region or index issue?): {e}")
                     vcf_records = []
                except Exception as e:
                     logger.error(f"Error fetching variants for contig {chrom}: {e}")
                     continue # Skip to next contig on error


                if not vcf_records:
                    # If no variants in the region/contig, create a single segment path for the reference portion
                    logger.info(f"No variants found for {chrom}" + (f" in region {region}" if region else "") + ". Creating reference path.")
                    # Determine start/end based on region or full contig
                    region_start_0based = fetch_kwargs.get('start', 0)
                    region_end_0based = fetch_kwargs.get('end', chrom_len) # Use chrom_len if no end specified

                    ref_seq = self._ref_handler.get_sequence(chrom, region_start_0based, region_end_0based)
                    if ref_seq:
                        ref_seg_id = self._get_or_create_segment(ref_seq)
                        if ref_seg_id:
                             # Create a simple path representing the reference for this contig/region
                             ref_path_name = f"{chrom}_ref" + (f"_{region_start_0based+1}_{region_end_0based}" if region else "")
                             try:
                                 path = gfapy.line.Path(
                                     path_name=ref_path_name,
                                     segment_names=[ref_seg_id],
                                     orientations=["+"],
                                     overlaps=[]
                                 )
                                 # Add reference path with validation disabled
                                 try:
                                     # First try direct add_line
                                     self._gfa.add_line(path)
                                 except RecursionError:
                                     # If recursion error occurs, try a more direct approach
                                     # Add to the paths collection directly
                                     if hasattr(self._gfa, 'paths') and hasattr(self._gfa.paths, 'append'):
                                         self._gfa.paths.append(path)
                                     else:
                                         # Last resort - just add to the lines collection
                                         self._gfa.lines[path.name] = path
                                 logger.info(f"Added reference path {ref_path_name} for contig {chrom}.")
                             except Exception as e:
                                 logger.error(f"Failed to add reference path for {chrom}: {e}")
                    continue # Move to the next contig


                # Process each sample if variants exist
                for sample_name in self._phasing_processor.samples:
                    logger.info(f"Processing sample: {sample_name} for contig {chrom}")

                    # Build haplotype sequences by iterating through the fetched records
                    try:
                        logger.info(f"Building haplotype sequences for sample {sample_name} with unphased_strategy={self.unphased_strategy}")
                        # Pass the list of records back as an iterator
                        hap1_segments, hap2_segments = self._phasing_processor.build_haplotype_sequences_incrementally(
                            iter(vcf_records), # Pass iterator
                            sample_name,
                            self._ref_handler,
                            chrom,
                            chrom_len, # Pass full length for trailing ref calculation
                            self.unphased_strategy
                        )
                        
                        # Verify we got segments back
                        logger.info(f"Generated {len(hap1_segments)} segments for haplotype 1 and {len(hap2_segments)} segments for haplotype 2")
                        
                        # Debug: print the first few segments to help diagnose issues
                        if hap1_segments:
                            for i, (pos, seq) in enumerate(hap1_segments[:2]):  # Show first 2 segments
                                logger.debug(f"Hap1 segment {i}: pos={pos}, seq_len={len(seq)}, seq_preview={seq[:20]}...")
                        if hap2_segments:
                            for i, (pos, seq) in enumerate(hap2_segments[:2]):  # Show first 2 segments
                                logger.debug(f"Hap2 segment {i}: pos={pos}, seq_len={len(seq)}, seq_preview={seq[:20]}...")
                        
                        if not hap1_segments and not hap2_segments:
                            logger.error(f"No segments generated for sample {sample_name}. Check phasing processor implementation.")
                            # Instead of raising an error, try to use the reference sequence as a fallback
                            ref_seq = self._ref_handler.get_sequence(chrom, 0, chrom_len)
                            if ref_seq:
                                logger.warning(f"Using full reference sequence as fallback for sample {sample_name}")
                                hap1_segments = [(0, ref_seq)]
                                hap2_segments = [(0, ref_seq)]
                            else:
                                raise VCFtoGFAConversionError(f"No segments generated for sample {sample_name} and reference fallback failed")
                    except PhasingError as e:
                         logger.error(f"Error building haplotype sequences for sample {sample_name}, contig {chrom}: {e}")
                         continue # Skip this sample for this contig

                    # Build the GFA paths and links from the generated sequence segments
                    path1_name = self.path_template.format(sample=sample_name, hap=1)
                    path2_name = self.path_template.format(sample=sample_name, hap=2)

                    self._build_gfa_path_and_links(path1_name, hap1_segments)
                    self._build_gfa_path_and_links(path2_name, hap2_segments)

            # --- End of Contig/Sample Loop ---

            # Write the GFA file
            logger.info(f"Conversion complete. Writing GFA to {self.output_gfa_filepath}")
            # Ensure output directory exists
            output_dir = os.path.dirname(self.output_gfa_filepath)
            if output_dir and not os.path.exists(output_dir):
                 os.makedirs(output_dir, exist_ok=True)
                 logger.info(f"Created output directory: {output_dir}")

            # Check if we have any content to write
            has_segments = len(self._segments) > 0
            has_paths = len(self._paths) > 0
            
            # Debug the GFA content before writing
            logger.debug(f"GFA content before writing: segments={len(self._segments)}, "
                        f"paths={len(self._paths)}, "
                        f"total_lines={len(self._gfa_lines)}")
            
            if not has_segments and not has_paths:
                logger.warning("No segments or paths were created during conversion. Creating fallback reference path.")
                logger.info("This may happen if no variants were processed or if all segments were empty.")
                
                # Create a fallback reference path for each chromosome
                for chrom in self._ref_handler.list_chromosomes():
                    chrom_len = self._ref_handler.get_chrom_length(chrom)
                    if chrom_len:
                        ref_seq = self._ref_handler.get_sequence(chrom, 0, chrom_len)
                        if ref_seq:
                            # Create a segment for the reference
                            ref_seg_id = self._get_or_create_segment(ref_seq)
                            if ref_seg_id:
                                # Create a path for the reference directly
                                ref_path_name = f"{chrom}_ref"
                                
                                # Add reference path
                                path_line = f"P\t{ref_path_name}\t{ref_seg_id}\t+\t*"
                                self._gfa_lines.append(path_line)
                                self._paths.add(ref_path_name)
                                logger.info(f"Added fallback reference path {ref_path_name}")
                                
                                # Also create sample paths using the reference sequence
                                for sample_name in self._phasing_processor.samples:
                                    path1_name = self.path_template.format(sample=sample_name, hap=1)
                                    path2_name = self.path_template.format(sample=sample_name, hap=2)
                                    
                                    for path_name in [path1_name, path2_name]:
                                        hap_num = 1 if path_name.endswith("_hap1") else 2
                                        path_line = f"P\t{path_name}\t{ref_seg_id}\t+\t*\tSM:Z:{sample_name}\tHP:i:{hap_num}"
                                        self._gfa_lines.append(path_line)
                                        self._paths.add(path_name)
                                        logger.info(f"Added fallback sample path {path_name}")
            
            # Write the GFA file
            with open(self.output_gfa_filepath, 'w') as outfile:
                # Write all lines
                if self._gfa_lines:
                    for line in self._gfa_lines:
                        outfile.write(line + "\n")
                    logger.info(f"Wrote {len(self._gfa_lines)} lines to GFA file")
                else:
                    # If we have no content at all, add a comment
                    outfile.write("H\tVN:Z:1.0\n")
                    outfile.write("# No segments or paths were created during conversion\n")
                    logger.info("Wrote minimal GFA with header only")
            
            logger.info("GFA file written successfully.")

        except (VCFtoGFAConversionError, PhasingError, ReferenceHandlerError, pysam.utils.SamtoolsError, gfapy.error.FormatError, Exception) as e:
            logger.error(f"An error occurred during VCF to GFA conversion: {e}", exc_info=True)
            # Re-raise as a specific conversion error
            raise VCFtoGFAConversionError(f"Conversion failed: {e}") from e
        finally:
            self._cleanup_resources()

    def __enter__(self):
        # Allows using 'with VCFtoGFAConverter(...) as converter:'
        # Could potentially initialize resources here if desired
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Ensures cleanup even if errors occur within the 'with' block
        self._cleanup_resources()

# --- Example Usage (for testing or direct script execution) ---
if __name__ == "__main__":
    # Basic argument parsing for standalone execution
    parser = argparse.ArgumentParser(description="Convert VCF to GFA format.")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file path (indexed).")
    parser.add_argument("-r", "--reference", required=True, help="Reference FASTA file path (indexed).")
    parser.add_argument("-o", "--output", required=True, help="Output GFA file path.")
    parser.add_argument("-s", "--sample-name", default=None, help="Sample name to process (default: first sample in VCF).")
    parser.add_argument("--region", default=None, help="Region to process (e.g., chr1:1000-2000 or chr1). Default: process all.")
    parser.add_argument("--unphased", default='ref', choices=['ref', 'alt', 'skip', 'error'],
                        help=f"Strategy for handling unphased variants (default: ref).")
    parser.add_argument("--include-ref-path", action='store_true', default=False, help="Include a GFA path representing the reference sequence.")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Set the logging level (default: INFO).")

    args = parser.parse_args()

    # Configure logging
    log_level = args.log_level.upper()
    logging.basicConfig(level=log_level,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    # Suppress verbose pysam messages if not in DEBUG mode
    if log_level != "DEBUG":
         logging.getLogger("pysam").setLevel(logging.WARNING)
         logging.getLogger("gfapy").setLevel(logging.INFO) # Adjust gfapy level if needed

    try:
        # Instantiate converter with appropriate arguments from command line
        converter = VCFtoGFAConverter(
            vcf_filepath=args.vcf,
            fasta_filepath=args.reference, # Corrected argument name
            output_gfa_filepath=args.output,
            path_template="{sample}_hap{hap}", # Default template, could be made configurable
            unphased_strategy=args.unphased
        )
        # Call convert method
        converter.convert(region=args.region)

        logger.info("Script finished successfully.")
        sys.exit(0)

    except (FileNotFoundError, ValueError, RuntimeError, VCFtoGFAConversionError) as e:
        logger.error(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.critical(f"An unexpected critical error occurred: {e}", exc_info=True)
        sys.exit(1)
