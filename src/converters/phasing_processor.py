# src/converters/phasing_processor.py
import logging
from typing import Dict, List, Optional, Tuple, Set
import pysam
from .reference_handler import ReferenceHandler # Use forward reference for type hint

logger = logging.getLogger(__name__)

class PhasingError(Exception):
    """Custom exception for phasing processing errors."""
    pass

class PhasingProcessor:
    """
    Processes phasing information from VCF records for samples.

    Extracts haplotype alleles based on GT and PS fields using pysam.
    Groups variants by phase set.
    """

    def __init__(self, vcf_reader: pysam.VariantFile):
        """
        Initializes the PhasingProcessor.

        Args:
            vcf_reader: An opened pysam.VariantFile object.
        """
        if not isinstance(vcf_reader, pysam.VariantFile):
            raise TypeError("vcf_reader must be a pysam.VariantFile instance.")
        self.vcf_reader = vcf_reader
        self.samples = list(vcf_reader.header.samples)
        logger.info(f"PhasingProcessor initialized for samples: {self.samples}")

    def get_sample_haplotype_alleles(self, record: pysam.VariantRecord, sample_name: str) -> Optional[Tuple[Optional[int], Optional[int]]]:
        """
        Gets the allele indices for the two haplotypes of a sample at a given record.

        Handles phased (e.g., 0|1) and unphased (e.g., 0/1) genotypes.
        Returns None for missing genotypes (./.) or if the sample is not found.

        Args:
            record: A pysam.VariantRecord object.
            sample_name: The name of the sample.

        Returns:
            A tuple containing two elements (haplotype1_allele_index, haplotype2_allele_index).
            Allele indices are 0 for REF, 1 for ALT1, 2 for ALT2, etc.
            Returns (None, None) if genotype is missing or invalid (e.g., non-diploid).
            Returns None if the sample is not in the VCF record's samples.
        """
        if sample_name not in record.samples:
            logger.warning(f"Sample '{sample_name}' not found in record {record.id or record.pos}.")
            return None

        sample_data = record.samples[sample_name]
        genotype = sample_data.get('GT')

        # Check for missing genotype (tuple like (None,) or (None, None))
        if genotype is None or any(a is None for a in genotype):
            # logger.debug(f"Missing genotype for sample '{sample_name}' at {record.chrom}:{record.pos}")
            return (None, None)

        # Check for ploidy (currently assuming diploid)
        if len(genotype) != 2:
            logger.warning(f"Non-diploid genotype {genotype} found for sample '{sample_name}' at {record.chrom}:{record.pos}. Treating as missing.")
            return (None, None)

        hap1_allele_idx, hap2_allele_idx = genotype
        # phased = sample_data.phased # We don't strictly need 'phased' here, just the GT

        # logger.debug(f"Sample {sample_name} at {record.chrom}:{record.pos}: GT={genotype}, Phased={sample_data.phased}")
        return (hap1_allele_idx, hap2_allele_idx)


    def get_phase_set(self, record: pysam.VariantRecord, sample_name: str) -> Optional[str]:
        """
        Gets the phase set identifier (PS tag) for a sample at a given record.

        Args:
            record: A pysam.VariantRecord object.
            sample_name: The name of the sample.

        Returns:
            The phase set identifier (usually an integer as a string) or None if
            the PS tag is not present, the sample is not found, or the genotype is unphased.
        """
        if sample_name not in record.samples:
            return None

        sample_data = record.samples[sample_name]
        # Only consider PS tag if the genotype is actually phased
        if not sample_data.phased:
            return None

        try:
            # Use sample.get() which handles missing tags gracefully
            ps_tag = sample_data.get('PS')
            if ps_tag is not None:
                # PS is typically an integer, convert to string for consistent key usage
                return str(ps_tag)
            else:
                # Phased genotype but no PS tag. This is unusual but possible.
                # Treat as unphased for block grouping purposes.
                # logger.debug(f"Phased genotype for {sample_name} at {record.chrom}:{record.pos} lacks PS tag.")
                return None
        except KeyError:
             # Should not happen with sample.get(), but as safety
             logger.warning(f"Could not retrieve PS tag for sample {sample_name} at {record.chrom}:{record.pos}, though sample exists.")
             return None


    def get_allele_sequence(self, record: pysam.VariantRecord, allele_index: Optional[int]) -> Optional[str]:
        """
        Gets the actual sequence for a given allele index (0=REF, 1=ALT1, etc.).

        Args:
            record: A pysam.VariantRecord object.
            allele_index: The index of the allele (0, 1, 2, ...).

        Returns:
            The allele sequence string (uppercase), or None if the index is invalid or None.
        """
        if allele_index is None:
            return None

        # record.alleles provides a tuple of (REF, ALT1, ALT2, ...)
        try:
            if 0 <= allele_index < len(record.alleles):
                allele_seq = record.alleles[allele_index]
                # Ensure sequence is returned, not None if allele is missing/invalid in VCF record itself
                return str(allele_seq).upper() if allele_seq is not None else None
            else:
                logger.warning(f"Invalid allele index {allele_index} requested for record at {record.chrom}:{record.pos} with {len(record.alleles)} alleles: {record.alleles}")
                return None
        except IndexError:
             logger.warning(f"IndexError accessing allele index {allele_index} for record at {record.chrom}:{record.pos} with alleles {record.alleles}")
             return None


    def group_variants_by_phase_set(self, records: List[pysam.VariantRecord], sample_name: str) -> Dict[Optional[str], List[pysam.VariantRecord]]:
        """
        Groups a list of VCF records by phase set (PS tag) for a specific sample.

        Variants without a PS tag or that are unphased for the sample are grouped
        under the key None.

        Args:
            records: A list of pysam.VariantRecord objects (presumably from the same contig).
            sample_name: The name of the sample to process phasing for.

        Returns:
            A dictionary where keys are phase set identifiers (strings, or None for unphased/no PS)
            and values are lists of records belonging to that phase set for the sample.
            Records within each list are sorted by position.
        """
        phased_groups: Dict[Optional[str], List[pysam.VariantRecord]] = {}
        for record in records:
            if sample_name not in record.samples:
                continue # Skip records where the sample has no data

            ps_id = self.get_phase_set(record, sample_name)
            # If ps_id is None, it means unphased or no PS tag for this sample

            if ps_id not in phased_groups:
                phased_groups[ps_id] = []
            phased_groups[ps_id].append(record)

        # Sort records within each group by position
        for ps_id in phased_groups:
            phased_groups[ps_id].sort(key=lambda r: r.pos)

        logger.info(f"Grouped {len(records)} records for sample '{sample_name}' into {len(phased_groups)} phase sets (including None).")
        # logger.debug(f"Phase set keys for sample '{sample_name}': {list(phased_groups.keys())}")

        return phased_groups

    def build_haplotype_sequences_incrementally(
        self,
        records_iterator: iter, # Iterator over pysam.VariantRecord
        sample_name: str,
        reference_handler: ReferenceHandler,
        chrom: str,
        chrom_len: int,
        unphased_strategy: str = 'ref'
    ) -> Tuple[List[Tuple[int, str]], List[Tuple[int, str]]]:
        """
        Constructs the two haplotype sequences by iterating through VCF records.

        This method processes records one by one, adding intervening reference
        segments and variant alleles based on the sample's genotype.

        Args:
            records_iterator: An iterator yielding pysam.VariantRecord objects for the chromosome, sorted by position.
            sample_name: The sample name.
            reference_handler: An initialized ReferenceHandler instance.
            chrom: The chromosome/contig name.
            chrom_len: The length of the chromosome.
            unphased_strategy: How to handle unphased/missing genotypes ('ref', 'alt', 'skip').

        Returns:
            A tuple containing two lists: (hap1_segments, hap2_segments).
            Each list contains tuples of (start_pos, sequence_string) for that haplotype.
            Positions are 0-based relative to the chromosome.
            Returns empty lists if errors occur.

        Raises:
            PhasingError: If reference sequence fetching fails.
        """
        hap1_segments: List[Tuple[int, str]] = []
        hap2_segments: List[Tuple[int, str]] = []
        last_pos_processed = 0 # 0-based coordinate tracking the end of the last processed segment/variant

        logger.info(f"Building haplotype sequences incrementally for sample '{sample_name}', chrom '{chrom}', unphased_strategy='{unphased_strategy}'")

        for record in records_iterator:
            # VCF position is 1-based, convert to 0-based start
            record_start_0based = record.start # pysam record.start is 0-based
            ref_allele = record.ref.upper() if record.ref else ""
            ref_len = len(ref_allele)
            # The reference region affected by the variant ends at start + ref_len
            record_ref_end_0based = record_start_0based + ref_len

            # --- Sanity checks ---
            if record_start_0based < last_pos_processed:
                 logger.warning(f"VCF record at {chrom}:{record.pos} starts before the end of the previous record ({last_pos_processed}). This might indicate overlapping variants or unsorted VCF. Skipping record.")
                 continue # Skip overlapping/unsorted variants

            # 1. Add intervening reference segment
            if record_start_0based > last_pos_processed:
                intervening_ref = reference_handler.get_sequence(chrom, last_pos_processed, record_start_0based)
                if intervening_ref is None:
                    msg = f"Failed to get intervening reference sequence for {chrom}:{last_pos_processed}-{record_start_0based}"
                    logger.error(msg)
                    raise PhasingError(msg)
                if intervening_ref: # Only add if not empty
                    hap1_segments.append((last_pos_processed, intervening_ref))
                    hap2_segments.append((last_pos_processed, intervening_ref))
                # Update last_pos_processed even if intervening_ref is empty
                last_pos_processed = record_start_0based


            # 2. Get haplotype alleles for this variant
            hap_alleles_idx = self.get_sample_haplotype_alleles(record, sample_name)
            hap1_seq, hap2_seq = None, None
            sample_data = record.samples[sample_name]
            is_phased = sample_data.phased
            
            # Log variant details for debugging
            logger.debug(f"Processing variant at {chrom}:{record.pos} for {sample_name}: "
                         f"REF={record.ref}, ALT={record.alts}, "
                         f"GT={sample_data.get('GT')}, Phased={is_phased}")

            if hap_alleles_idx is None or hap_alleles_idx == (None, None):
                # Handle missing/invalid genotype based on strategy
                logger.debug(f"Missing/invalid GT for {sample_name} at {chrom}:{record.pos}. Strategy: {unphased_strategy}")
                if unphased_strategy == 'ref':
                    hap1_seq = ref_allele
                    hap2_seq = ref_allele
                elif unphased_strategy == 'alt':
                    # Use first ALT if available, otherwise REF
                    alt1_seq = self.get_allele_sequence(record, 1)
                    hap1_seq = alt1_seq if alt1_seq is not None else ref_allele
                    hap2_seq = alt1_seq if alt1_seq is not None else ref_allele
                # else 'skip': sequences remain None, effectively using reference up to this point

            else:
                # Handle phased or unphased genotype
                hap1_idx, hap2_idx = hap_alleles_idx

                if not is_phased:
                    logger.debug(f"Unphased GT {hap1_idx}/{hap2_idx} for {sample_name} at {chrom}:{record.pos}. Strategy: {unphased_strategy}")
                    if unphased_strategy == 'ref':
                        hap1_seq = self.get_allele_sequence(record, 0) # REF
                        hap2_seq = self.get_allele_sequence(record, 0) # REF
                    elif unphased_strategy == 'alt':
                        # Use the actual alleles from the unphased genotype
                        # This is the key fix - we should use the actual GT values, not just ALT1
                        hap1_seq = self.get_allele_sequence(record, hap1_idx)
                        hap2_seq = self.get_allele_sequence(record, hap2_idx)
                        
                        # If either allele is missing, use the first ALT as fallback
                        if hap1_seq is None or hap2_seq is None:
                            alt1_seq = self.get_allele_sequence(record, 1)
                            hap1_seq = alt1_seq if hap1_seq is None else hap1_seq
                            hap2_seq = alt1_seq if hap2_seq is None else hap2_seq
                            
                            # If ALT is also None, fall back to REF
                            hap1_seq = ref_allele if hap1_seq is None else hap1_seq
                            hap2_seq = ref_allele if hap2_seq is None else hap2_seq
                    elif unphased_strategy == 'skip':
                         hap1_seq = None
                         hap2_seq = None
                    else: # Default to REF for safety
                         hap1_seq = ref_allele
                         hap2_seq = ref_allele

                    # Ensure sequences are not None if strategy wasn't skip
                    if unphased_strategy != 'skip':
                         hap1_seq = ref_allele if hap1_seq is None else hap1_seq
                         hap2_seq = ref_allele if hap2_seq is None else hap2_seq

                else: # Phased genotype
                    hap1_seq = self.get_allele_sequence(record, hap1_idx)
                    hap2_seq = self.get_allele_sequence(record, hap2_idx)

                    if hap1_seq is None or hap2_seq is None:
                        logger.warning(f"Could not get allele sequence for phased GT {hap1_idx}|{hap2_idx} for {sample_name} at {chrom}:{record.pos}. Using REF as fallback.")
                        hap1_seq = ref_allele if hap1_seq is None else hap1_seq
                        hap2_seq = ref_allele if hap2_seq is None else hap2_seq


            # 3. Add allele sequences if not skipped
            # Ensure sequences are strings (might be None if skipped or error)
            hap1_seq_str = hap1_seq if isinstance(hap1_seq, str) else ""
            hap2_seq_str = hap2_seq if isinstance(hap2_seq, str) else ""

            # Only add non-empty sequences. Deletions result in empty strings here.
            if hap1_seq_str:
                hap1_segments.append((record_start_0based, hap1_seq_str))
            if hap2_seq_str:
                hap2_segments.append((record_start_0based, hap2_seq_str))

            # 4. Update last processed position to the end of the reference allele span
            last_pos_processed = record_ref_end_0based


        # Add final reference segment from last variant to end of chromosome
        if last_pos_processed < chrom_len:
            trailing_ref = reference_handler.get_sequence(chrom, last_pos_processed, chrom_len)
            if trailing_ref is None:
                 msg = f"Failed to get trailing reference sequence for {chrom}:{last_pos_processed}-{chrom_len}"
                 logger.error(msg)
                 raise PhasingError(msg)
            if trailing_ref:
                hap1_segments.append((last_pos_processed, trailing_ref))
                hap2_segments.append((last_pos_processed, trailing_ref))

        logger.info(f"Finished building sequences for sample '{sample_name}', chrom '{chrom}'. Hap1 segments: {len(hap1_segments)}, Hap2 segments: {len(hap2_segments)}")
        return hap1_segments, hap2_segments
