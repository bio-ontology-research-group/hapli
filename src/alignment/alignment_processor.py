"""
Processor for the two-phase feature alignment strategy.
"""
import logging
from typing import Dict, List, Optional, Set, Tuple, Union

import networkx as nx
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq # Import Seq for reverse_complement

# Assuming these parsers and graph classes are correctly imported if needed
# from src.parsers.feature_graph import FeatureGraph
# from src.parsers.gfa_parser import GFAParser
# from src.parsers.gff_parser import GFF3Parser
# from src.parsers.fasta_parser import FastaParser
from src.alignment.minimap_wrapper import MinimapAligner
from src.parsers.feature_graph import FeatureGraph # Needed for type hint

logger = logging.getLogger(__name__)

class AlignmentProcessor:
    """
    Processor for aligning genomic features to paths in a variation graph.

    This class implements a two-phase alignment strategy:
    1. First phase: Aligns parent features (e.g., genes) to target paths
    2. Second phase: Aligns child features (e.g., exons) within parent boundaries

    It handles feature duplications and multiple alignments.
    """

    def __init__(self, gfa_parser, gff_parser, fasta_parser, feature_graph: FeatureGraph, minimap_preset: str = "splice"):
        """
        Initialize the alignment processor with pre-loaded parsers and graph.

        Args:
            gfa_parser: An initialized and loaded GFAParser instance.
            gff_parser: An initialized and loaded GFF3Parser instance.
            fasta_parser: An initialized and loaded FastaParser instance.
            feature_graph: A FeatureGraph instance built from GFF features.
            minimap_preset: The preset to use for minimap2 alignment (e.g., 'splice', 'map-ont').
        """
        self.gfa_parser = gfa_parser
        self.gff_parser = gff_parser
        self.fasta_parser = fasta_parser
        self.feature_graph = feature_graph
        self.aligner = MinimapAligner(preset=minimap_preset)
        logger.info(f"AlignmentProcessor initialized with minimap2 preset: {minimap_preset}")

        # Results storage
        self.path_sequences: Dict[str, str] = {}  # Path ID -> sequence
        self.aligned_features: Dict[str, Dict[str, List[SeqFeature]]] = {}  # Path ID -> {feature_id: [aligned_features]}

    def extract_path_sequences(self, path_ids: List[str]):
        """
        Extract sequences for specified paths using the loaded GFA parser.

        Args:
            path_ids: List of path IDs to extract sequences for
        """
        if not self.gfa_parser:
            logger.error("GFA parser not initialized. Cannot extract path sequences.")
            return

        logger.info(f"Extracting sequences for {len(path_ids)} paths...")
        extracted_count = 0
        for path_id in path_ids:
            # Get the path from the GFA
            paths = self.gfa_parser.get_paths()
            if path_id not in paths:
                logger.warning(f"Path {path_id} not found in GFA")
                continue

            # Extract segment IDs from the path
            path = paths[path_id]
            segments = [] # List of segment IDs (str)

            # Handle different GFA versions/representations from GFApy
            try:
                # GFApy representation (often list of segment names/objects)
                if hasattr(path, 'segment_names'): # GFA1-like attribute
                    for seg_ref in path.segment_names:
                        if isinstance(seg_ref, str):
                            segments.append(seg_ref.strip('+').strip('-'))
                        elif hasattr(seg_ref, 'name'): # GFApy segment object
                            segments.append(str(seg_ref.name))
                        else: # Fallback: try converting to string
                            seg_str = str(seg_ref).strip('+').strip('-')
                            if seg_str: segments.append(seg_str)
                elif hasattr(path, 'items'): # GFA2-like attribute
                     for item in path.items:
                         if hasattr(item, 'name'): # GFApy segment object
                             segments.append(str(item.name).strip('+').strip('-'))
                         else: # Fallback: try converting to string
                             item_str = str(item).strip('+').strip('-')
                             if item_str: segments.append(item_str)
                elif isinstance(path, str) and ',' in path: # Simple comma-separated string path
                     for seg in path.split(','):
                         segments.append(seg.strip('+').strip('-'))
                else:
                     # Attempt to get segments via gfapy's path object methods if available
                     if hasattr(path, 'segments'):
                         segs_list = path.segments # This might be the list we need
                         if isinstance(segs_list, list):
                              for seg_obj in segs_list:
                                   if hasattr(seg_obj, 'name'):
                                       segments.append(str(seg_obj.name))
                                   else:
                                       seg_str = str(seg_obj).strip('+').strip('-')
                                       if seg_str: segments.append(seg_str)

            except Exception as e:
                logger.warning(f"Error extracting segments from path {path_id} (type: {type(path)}): {e}")
                # Continue trying to get sequence below, maybe segments were partially extracted

            if not segments:
                logger.warning(f"Could not determine segments for path {path_id}. Path representation: {path}")
                continue

            # Concatenate segment sequences
            path_seq_list = []
            valid_sequence = True
            for seg_id in segments:
                seg_seq = None
                try:
                    # Use the provided GFA parser instance
                    seg_seq = self.gfa_parser.get_segment_sequence(seg_id)
                except Exception as e:
                    logger.warning(f"Error getting sequence for segment {seg_id} using get_segment_sequence: {e}")

                if seg_seq:
                    path_seq_list.append(seg_seq)
                else:
                    logger.warning(f"No sequence found for segment {seg_id} in path {path_id}")
                    valid_sequence = False
                    break # Stop processing this path if a segment is missing

            if valid_sequence:
                full_path_seq = "".join(path_seq_list)
                if full_path_seq:
                    self.path_sequences[path_id] = full_path_seq
                    logger.debug(f"Extracted sequence of length {len(full_path_seq)} for path {path_id}")
                    extracted_count += 1
                else:
                    logger.warning(f"Failed to construct sequence for path {path_id} (empty result)")
            else:
                 logger.warning(f"Failed to extract full sequence for path {path_id} due to missing segments.")

        logger.info(f"Successfully extracted sequences for {extracted_count}/{len(path_ids)} requested paths.")


    def align_features_to_paths(self,
                              path_ids: List[str],
                              feature_types: Optional[List[str]] = None,
                              min_identity: float = 0.8,
                              min_coverage: float = 0.8):
        """
        Align features to specified paths using the two-phase strategy.

        Args:
            path_ids: List of path IDs to align features to. Sequences must be pre-extracted.
            feature_types: Optional list of top-level feature types to align (default: genes).
            min_identity: Minimum sequence identity for an alignment to be considered valid.
            min_coverage: Minimum coverage of the feature sequence for an alignment to be considered valid.

        Returns:
            Dictionary mapping path IDs to dictionaries of feature IDs to lists of aligned features.
            The aligned features (SeqFeature objects) contain alignment details in their qualifiers.
        """
        if not feature_types:
            feature_types = ["gene"] # Default to aligning genes as top-level features

        # Ensure path sequences are available for all requested paths
        missing_seq_paths = [pid for pid in path_ids if pid not in self.path_sequences]
        if missing_seq_paths:
            logger.warning(f"Missing sequences for paths: {missing_seq_paths}. Extract sequences first or these paths will be skipped.")
            # Filter path_ids to only include those with sequences
            path_ids = [pid for pid in path_ids if pid in self.path_sequences]
            if not path_ids:
                 logger.error("No paths with sequences available for alignment.")
                 return {}

        # Reset results for the current run
        self.aligned_features = {path_id: {} for path_id in path_ids}

        # --- Phase 1: Align Parent Features ---
        logger.info("--- Starting Alignment Phase 1: Parent Features ---")
        self._align_parent_features(path_ids, feature_types, min_identity, min_coverage)
        logger.info("--- Finished Alignment Phase 1 ---")

        # --- Phase 2: Align Child Features ---
        logger.info("--- Starting Alignment Phase 2: Child Features ---")
        self._align_child_features(path_ids, min_identity, min_coverage)
        logger.info("--- Finished Alignment Phase 2 ---")

        return self.aligned_features

    def _get_feature_sequence(self, feature: SeqFeature) -> Optional[Tuple[str, str]]:
        """Helper to extract reference sequence ID and sequence for a feature."""
        ref_id = None
        # Try standard GFF attributes first (often stored in qualifiers by parsers)
        if 'seqid' in feature.qualifiers:
            ref_id = feature.qualifiers['seqid'][0]
        elif hasattr(feature, 'ref'): # Biopython often uses 'ref'
            ref_id = feature.ref
        elif hasattr(feature, 'seq_id'): # Or 'seq_id'
            ref_id = feature.seq_id

        # Fallback: check common qualifier keys
        if not ref_id and hasattr(feature, 'qualifiers'):
            for key in ['chromosome', 'contig', 'sequence_id', 'ref_seq']:
                if key in feature.qualifiers:
                    ref_id = feature.qualifiers[key][0]
                    break

        if not ref_id:
            logger.debug(f"Could not determine reference sequence ID for feature {getattr(feature, 'id', 'UnknownID')}")
            return None

        # Get reference sequence using the provided FastaParser
        ref_seq_record = self.fasta_parser.get_sequence(ref_id)
        if not ref_seq_record:
            logger.warning(f"Reference sequence '{ref_id}' not found in FASTA for feature {getattr(feature, 'id', 'UnknownID')}")
            return None

        # Extract the feature sequence
        try:
            start = int(feature.location.start)
            end = int(feature.location.end)
            ref_len = len(ref_seq_record.seq)

            # Basic bounds check
            if start < 0 or end > ref_len or start >= end:
                 logger.warning(f"Invalid coordinates for feature {getattr(feature, 'id', 'UnknownID')} on ref '{ref_id}': start={start}, end={end}, ref_len={ref_len}. Skipping.")
                 return None

            feature_seq = str(ref_seq_record.seq[start:end])

            # Handle strand if necessary (e.g., for reverse complement)
            if feature.location.strand == -1:
                 feature_seq = str(Seq(feature_seq).reverse_complement()) # Requires from Bio.Seq import Seq

            if not feature_seq:
                 logger.warning(f"Extracted empty sequence for feature {getattr(feature, 'id', 'UnknownID')} at {ref_id}:{start}-{end}")
                 return None

            return ref_id, feature_seq

        except Exception as e:
            logger.error(f"Error extracting sequence for feature {getattr(feature, 'id', 'UnknownID')} from ref '{ref_id}': {e}", exc_info=True)
            return None


    def _align_parent_features(self, path_ids: List[str], feature_types: List[str], min_identity: float, min_coverage: float):
        """
        Phase 1: Align top-level features (e.g., genes) to paths.
        """
        parent_features = []
        for feature_type in feature_types:
            # Use the provided GFF parser instance
            features_of_type = self.gff_parser.get_features_by_type(feature_type)
            # Filter for features that are likely parents (no parent in the graph or specified types)
            for feature in features_of_type:
                 feature_id = getattr(feature, 'id', None)
                 if not feature_id: continue # Skip features without ID

                 # Check if it's a root in the feature graph or explicitly listed type
                 is_root = self.feature_graph.graph.in_degree(feature_id) == 0 if feature_id in self.feature_graph.graph else True
                 if is_root:
                      parent_features.append(feature)

        logger.info(f"Identified {len(parent_features)} potential parent features of types {feature_types} for Phase 1 alignment.")

        # Align each parent feature to each path
        for path_id in path_ids:
            path_seq = self.path_sequences.get(path_id)
            if not path_seq:
                logger.warning(f"No sequence for path {path_id}, skipping parent alignment.")
                continue

            logger.debug(f"Aligning parent features to path {path_id} (length {len(path_seq)})")
            # Load path sequence into aligner
            self.aligner.load_reference(path_seq)

            for feature in parent_features:
                feature_id = getattr(feature, 'id', None)
                if not feature_id:
                    # Already warned during feature collection, but double-check
                    continue

                # Get feature sequence
                seq_info = self._get_feature_sequence(feature)
                if not seq_info:
                    continue # Error logged in helper
                ref_id, feature_seq = seq_info
                original_location_str = f"{ref_id}:{int(feature.location.start)}-{int(feature.location.end)}"

                # Align the feature sequence to the path sequence
                try:
                    # Use aligner's method - adjust min_score/min_len if needed
                    # mappy returns list of Hit objects
                    alignments = self.aligner.align_sequence(feature_seq, min_score=40, min_len=50)
                except Exception as e:
                    logger.error(f"Alignment failed for feature {feature_id} on path {path_id}: {e}", exc_info=True)
                    alignments = []

                if not alignments:
                    logger.debug(f"No significant alignments found for parent feature {feature_id} on path {path_id}")
                    continue

                # Process valid alignments
                aligned_features_list = []
                for aln in alignments:
                    # Basic filtering based on identity and coverage (mappy doesn't directly provide these)
                    # We can estimate identity from NM (edit distance) and query length
                    # Coverage needs alignment length on query vs query length
                    # For now, store raw alignment info. Downstream analysis needs to calculate/filter.
                    q_len = len(feature_seq)
                    # Identity approximation: (aln.mlen - aln.NM) / aln.mlen if aln.mlen > 0 else 0
                    # Coverage approximation: aln.q_en - aln.q_st / q_len if q_len > 0 else 0
                    # Let's store the CIGAR and score, downstream can refine

                    # Create a new SeqFeature for the aligned location on the path
                    # Use alignment coordinates on the reference (path sequence)
                    aligned_feature = SeqFeature(
                        location=FeatureLocation(aln.r_st, aln.r_en, strand=feature.location.strand), # Use original strand
                        type=feature.type,
                        id=feature_id,
                        qualifiers=feature.qualifiers.copy() # Copy original qualifiers
                    )

                    # Add alignment-specific details to qualifiers
                    aligned_feature.qualifiers['alignment_score'] = [str(aln.mapq)] # Mapping quality
                    aligned_feature.qualifiers['alignment_cigar'] = [aln.cigar_str]
                    aligned_feature.qualifiers['alignment_target_span'] = [f"{aln.r_st}-{aln.r_en}"]
                    aligned_feature.qualifiers['alignment_query_span'] = [f"{aln.q_st}-{aln.q_en}"]
                    aligned_feature.qualifiers['alignment_strand'] = [aln.strand] # Strand relative to target
                    aligned_feature.qualifiers['original_location'] = [original_location_str]
                    aligned_feature.qualifiers['is_primary'] = [aln.is_primary]
                    # Add NM (edit distance) and blen (alignment length) if needed for downstream calcs
                    aligned_feature.qualifiers['edit_distance_NM'] = [aln.NM]
                    aligned_feature.qualifiers['alignment_length_blen'] = [aln.blen]


                    aligned_features_list.append(aligned_feature)

                if aligned_features_list:
                    if feature_id not in self.aligned_features[path_id]:
                         self.aligned_features[path_id][feature_id] = []
                    self.aligned_features[path_id][feature_id].extend(aligned_features_list)
                    logger.debug(f"Stored {len(aligned_features_list)} alignments for parent feature {feature_id} on path {path_id}")


    def _align_child_features(self, path_ids: List[str], min_identity: float, min_coverage: float):
        """
        Phase 2: Align child features within the aligned boundaries of their parents.
        """
        # Iterate through paths and the parent features already aligned to them
        for path_id in path_ids:
            if path_id not in self.aligned_features or not self.aligned_features[path_id]:
                logger.debug(f"No parent features aligned to path {path_id}, skipping child alignment.")
                continue

            path_seq = self.path_sequences.get(path_id)
            if not path_seq:
                logger.warning(f"Sequence missing for path {path_id} during child alignment phase.")
                continue

            logger.debug(f"Aligning child features within parent boundaries on path {path_id}")

            # Get IDs of parents that were successfully aligned to this path
            aligned_parent_ids = list(self.aligned_features[path_id].keys())

            for parent_id in aligned_parent_ids:
                # Get the children of this parent from the feature graph
                try:
                    # Use the provided FeatureGraph instance
                    children_ids = list(self.feature_graph.get_children(parent_id))
                except Exception as e:
                    logger.warning(f"Could not get children for parent {parent_id}: {e}")
                    children_ids = []

                if not children_ids:
                    # This is normal for features like genes that might not have annotated children (e.g., non-coding)
                    # logger.debug(f"Parent feature {parent_id} has no children in the graph.")
                    continue

                logger.debug(f"Processing {len(children_ids)} children for parent {parent_id}")

                # Iterate through each *instance* of the aligned parent on this path
                # A parent might align multiple times
                parent_alignment_instances = self.aligned_features[path_id][parent_id]
                for parent_instance_idx, aligned_parent_feature in enumerate(parent_alignment_instances):

                    parent_path_start = int(aligned_parent_feature.location.start)
                    parent_path_end = int(aligned_parent_feature.location.end)
                    parent_instance_label = f"{parent_id}_instance_{parent_instance_idx+1}" # Unique label

                    # Extract the sequence of the aligned parent region from the path sequence
                    if parent_path_start < 0 or parent_path_end > len(path_seq) or parent_path_start >= parent_path_end:
                         logger.warning(f"Invalid parent instance coordinates for {parent_instance_label} on path {path_id}: {parent_path_start}-{parent_path_end}. Skipping child alignment for this instance.")
                         continue

                    parent_region_seq = path_seq[parent_path_start:parent_path_end]
                    if not parent_region_seq:
                         logger.warning(f"Extracted empty parent region sequence for {parent_instance_label} ({parent_path_start}-{parent_path_end}) on path {path_id}. Skipping.")
                         continue

                    # Load this specific parent region sequence as the reference for child alignment
                    try:
                        self.aligner.load_reference(parent_region_seq)
                    except Exception as e:
                         logger.error(f"Failed to load parent region sequence ({len(parent_region_seq)}bp) into aligner for {parent_instance_label}: {e}. Skipping children for this instance.")
                         continue

                    # Now, align each child feature to this parent region sequence
                    for child_id in children_ids:
                        # Get the original child feature definition from the GFF parser
                        child_feature = self.gff_parser.get_feature_by_id(child_id)
                        if not child_feature:
                            logger.warning(f"Child feature definition {child_id} (parent {parent_id}) not found in GFF parser.")
                            continue

                        # Get the child's original sequence from the reference genome
                        seq_info = self._get_feature_sequence(child_feature)
                        if not seq_info:
                            continue # Error logged in helper
                        ref_id, child_seq = seq_info
                        original_location_str = f"{ref_id}:{int(child_feature.location.start)}-{int(child_feature.location.end)}"

                        # Align child sequence to the extracted parent region sequence
                        try:
                            # Use stricter parameters for children? Maybe lower min_score?
                            child_alignments = self.aligner.align_sequence(child_seq, min_score=20, min_len=20)
                        except Exception as e:
                            logger.error(f"Alignment failed for child {child_id} within {parent_instance_label} on path {path_id}: {e}", exc_info=True)
                            child_alignments = []

                        if not child_alignments:
                            logger.debug(f"No significant alignments found for child {child_id} within {parent_instance_label} on path {path_id}")
                            continue

                        # Process valid child alignments
                        for aln in child_alignments:
                            # Coordinates are relative to the parent_region_seq
                            # Map them back to the full path coordinates
                            child_path_start = parent_path_start + aln.r_st
                            child_path_end = parent_path_start + aln.r_en

                            # Create a new SeqFeature for the aligned child on the path
                            aligned_child_feature = SeqFeature(
                                location=FeatureLocation(child_path_start, child_path_end, strand=child_feature.location.strand), # Use original strand
                                type=child_feature.type,
                                id=child_id,
                                qualifiers=child_feature.qualifiers.copy() # Copy original qualifiers
                            )

                            # Add alignment-specific details
                            aligned_child_feature.qualifiers['alignment_score'] = [str(aln.mapq)]
                            aligned_child_feature.qualifiers['alignment_cigar'] = [aln.cigar_str]
                            aligned_child_feature.qualifiers['alignment_target_span'] = [f"{aln.r_st}-{aln.r_en}"] # Relative to parent region
                            aligned_child_feature.qualifiers['alignment_query_span'] = [f"{aln.q_st}-{aln.q_en}"]
                            aligned_child_feature.qualifiers['alignment_strand'] = [aln.strand]
                            aligned_child_feature.qualifiers['original_location'] = [original_location_str]
                            aligned_child_feature.qualifiers['parent_feature'] = [parent_id] # Link to parent ID
                            aligned_child_feature.qualifiers['parent_instance'] = [parent_instance_label] # Link to specific parent alignment instance
                            aligned_child_feature.qualifiers['parent_path_location'] = [f"{parent_path_start}-{parent_path_end}"] # Parent location on path
                            aligned_child_feature.qualifiers['is_primary'] = [aln.is_primary]
                            aligned_child_feature.qualifiers['edit_distance_NM'] = [aln.NM]
                            aligned_child_feature.qualifiers['alignment_length_blen'] = [aln.blen]


                            # Store the aligned child feature, associated with the path
                            if child_id not in self.aligned_features[path_id]:
                                self.aligned_features[path_id][child_id] = []
                            self.aligned_features[path_id][child_id].append(aligned_child_feature)

                        logger.debug(f"Stored {len(child_alignments)} alignments for child {child_id} within {parent_instance_label} on path {path_id}")


    def get_alignment_results(self) -> Dict[str, Dict[str, List[SeqFeature]]]:
        """
        Get all alignment results.

        Returns:
            Dictionary mapping path IDs to dictionaries of feature IDs to lists of aligned SeqFeature objects.
            Each SeqFeature represents one alignment of a feature on that path.
        """
        return self.aligned_features

    def get_path_alignments(self, path_id: str) -> Dict[str, List[SeqFeature]]:
        """
        Get alignment results for a specific path.

        Args:
            path_id: Path ID to get alignments for.

        Returns:
            Dictionary mapping feature IDs to lists of aligned SeqFeature objects for the specified path.
            Returns an empty dict if the path was not processed or had no alignments.
        """
        return self.aligned_features.get(path_id, {})

    def get_feature_alignments(self, path_id: str, feature_id: str) -> List[SeqFeature]:
        """
        Get alignment results for a specific feature on a specific path.

        Args:
            path_id: Path ID to get alignments for.
            feature_id: Feature ID to get alignments for.

        Returns:
            List of aligned SeqFeature objects for the specified feature on the specified path.
            Returns an empty list if the feature was not aligned to the path.
        """
        return self.aligned_features.get(path_id, {}).get(feature_id, [])

    def export_alignments_tsv(self, output_file: str):
        """
        Export alignment results to a TSV file.

        Args:
            output_file: Path to output file.
        """
        logger.info(f"Exporting alignment results to {output_file}...")
        count = 0
        try:
            # Ensure directory exists
            output_dir = os.path.dirname(output_file)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)

            with open(output_file, 'w', newline='') as f: # Use newline='' for csv writer
                # Define header - include new fields
                header = [
                    "path_id", "feature_id", "type", "start", "end", "strand",
                    "parent_feature", "parent_instance", "parent_path_location",
                    "score", "cigar", "target_span", "query_span", "aln_strand",
                    "edit_distance_NM", "aln_length_blen", "is_primary",
                    "original_location"
                ]
                # Use csv.writer for proper TSV handling
                import csv
                writer = csv.writer(f, delimiter='\t', lineterminator='\n')
                writer.writerow(header)

                # Write alignments
                for path_id, features in self.aligned_features.items():
                    for feature_id, alignments in features.items():
                        for aln_feature in alignments:
                            # Get feature info from the aligned SeqFeature
                            feature_type = aln_feature.type
                            start = int(aln_feature.location.start)
                            end = int(aln_feature.location.end)
                            strand_val = aln_feature.location.strand
                            strand = '+' if strand_val == 1 else '-' if strand_val == -1 else '.'

                            # Get alignment info from qualifiers (use .get with default)
                            qual = aln_feature.qualifiers
                            parent = qual.get('parent_feature', ['.'])[0]
                            parent_instance = qual.get('parent_instance', ['N/A'])[0]
                            parent_loc = qual.get('parent_path_location', ['N/A'])[0]
                            score = qual.get('alignment_score', [''])[0]
                            cigar = qual.get('alignment_cigar', [''])[0]
                            target_span = qual.get('alignment_target_span', [''])[0]
                            query_span = qual.get('alignment_query_span', [''])[0]
                            aln_strand = qual.get('alignment_strand', [''])[0]
                            nm = qual.get('edit_distance_NM', [''])[0]
                            blen = qual.get('alignment_length_blen', [''])[0]
                            primary = qual.get('is_primary', [''])[0]
                            orig_loc = qual.get('original_location', [''])[0]

                            # Write row using csv writer
                            row = [
                                path_id, feature_id, feature_type, start, end, strand,
                                parent, parent_instance, parent_loc,
                                score, cigar, target_span, query_span, aln_strand,
                                nm, blen, primary,
                                orig_loc
                            ]
                            writer.writerow(row)
                            count += 1

            logger.info(f"Successfully exported {count} alignment records to {output_file}")
        except Exception as e:
            logger.error(f"Failed to export alignments to TSV: {e}", exc_info=True)

```
src/config.py
