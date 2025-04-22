# src/converters/vcf_to_gfa.py
import logging
import itertools
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
        """
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
        self._gfa: Optional[gfapy.Gfa] = None
        self._segment_id_counter = 0
        # Cache: sequence -> segment_id
        self._segment_cache: Dict[str, str] = {}
        # Cache: (from_seg, from_orient, to_seg, to_orient) -> True
        self._link_cache: Dict[Tuple[str, str, str, str], bool] = {}

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
            self._gfa = gfapy.Gfa()
            self._segment_id_counter = 0
            self._segment_cache = {}
            self._link_cache = {}
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
        self._gfa = None
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

                segment = gfapy.line.Segment(
                    name=segment_id,
                    sequence=sequence,
                    LN=len(sequence) # Add LN tag explicitly
                )
                self._gfa.add_line(segment)
                self._segment_cache[sequence] = segment_id
                # logger.debug(f"Created segment {segment_id} (len {len(sequence)})")
                return segment_id
            except Exception as e:
                logger.error(f"Failed to create GFA segment for sequence '{sequence[:20]}...': {e}")
                # Decide how to handle: re-raise, return None? Re-raising for clarity.
                raise VCFtoGFAConversionError(f"Failed to add segment: {e}") from e


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
            link = gfapy.line.Link(
                from_segment=from_segment_id,
                from_orient=from_orient,
                to_segment=to_segment_id,
                to_orient=to_orient,
                overlap=cigar # Use '*' as default overlap CIGAR
            )
            self._gfa.add_line(link)
            self._link_cache[link_key] = True # Add to cache
            # logger.debug(f"Added link: {from_segment_id}{from_orient} -> {to_segment_id}{to_orient}")
        except Exception as e:
            logger.error(f"Failed to add GFA link {from_segment_id}{from_orient} -> {to_segment_id}{to_orient}: {e}")
            raise VCFtoGFAConversionError(f"Failed to add link: {e}") from e

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

            path = gfapy.line.Path(
                path_name=path_name,
                segment_names=path_segment_ids,
                orientations=orientations,
                overlaps=overlaps
            )
            # Add sample and haplotype info as tags (optional but good practice)
            try:
                 # Basic parsing assuming template like "{sample}_hap{hap}" or "{sample}.{hap}" etc.
                 # This is a simple heuristic and might fail for complex templates.
                 hap_part = ""
                 sample_part = path_name
                 if "_hap" in path_name:
                     parts = path_name.rsplit('_hap', 1)
                     if len(parts) == 2 and parts[1].isdigit():
                         sample_part, hap_part = parts
                 elif "." in path_name:
                      parts = path_name.rsplit('.', 1)
                      if len(parts) == 2 and parts[1].isdigit():
                          sample_part, hap_part = parts

                 if sample_part != path_name: # Check if parsing likely succeeded
                     path.set_datatype("SM", "Z") # Sample Name
                     path.set_tag("SM", sample_part)
                 if hap_part.isdigit():
                     path.set_datatype("HP", "i") # Haplotype ID
                     path.set_tag("HP", int(hap_part))
                 else:
                      # If hap couldn't be parsed, maybe set default?
                      # path.set_datatype("HP", "i")
                      # path.set_tag("HP", 1) # Or skip tag
                      pass

            except Exception as tag_ex:
                 logger.warning(f"Could not parse sample/haplotype from path name '{path_name}' for tagging: {tag_ex}")

            self._gfa.add_line(path)
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
            # Add GFA header (optional but recommended)
            # Using GFA1 for broader compatibility unless GFA2 features are needed
            self._gfa.add_line(gfapy.line.Header(VN="1.0"))

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
                                 self._gfa.add_line(path)
                                 logger.info(f"Added reference path {ref_path_name} for contig {chrom}.")
                             except Exception as e:
                                 logger.error(f"Failed to add reference path for {chrom}: {e}")
                    continue # Move to the next contig


                # Process each sample if variants exist
                for sample_name in self._phasing_processor.samples:
                    logger.info(f"Processing sample: {sample_name} for contig {chrom}")

                    # Build haplotype sequences by iterating through the fetched records
                    try:
                        # Pass the list of records back as an iterator
                        hap1_segments, hap2_segments = self._phasing_processor.build_haplotype_sequences_incrementally(
                            iter(vcf_records), # Pass iterator
                            sample_name,
                            self._ref_handler,
                            chrom,
                            chrom_len, # Pass full length for trailing ref calculation
                            self.unphased_strategy
                        )
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

            with open(self.output_gfa_filepath, 'w') as outfile:
                 # Use gfapy's __str__ method for output
                 outfile.write(str(self._gfa))
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
