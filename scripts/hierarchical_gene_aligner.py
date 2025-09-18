#!/usr/bin/env python3

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

import gffutils
import pysam
from Bio import pairwise2

# --- Setup ---

def setup_logging(verbose: bool = False):
    """Set up logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

# --- Data Models ---

class AlignmentResult:
    """Stores the result of a single alignment."""
    def __init__(self, feature_id: str, feature_type: str, score: int, target_start: int, target_end: int, identity: float):
        self.feature_id = feature_id
        self.feature_type = feature_type
        self.score = score
        self.target_start = target_start
        self.target_end = target_end
        self.identity = identity
        self.children: List[AlignmentResult] = []

    def __repr__(self) -> str:
        return f"AlignmentResult(id={self.feature_id}, type={self.feature_type}, score={self.score})"

# --- Core Classes ---

class GFFProcessor:
    """Handles GFF parsing and feature hierarchy."""
    def __init__(self, gff_file: Path):
        self.db_path = gff_file.with_suffix(".db")
        if not self.db_path.exists():
            logging.info(f"Creating GFF database at {self.db_path}...")
            self.db = gffutils.create_db(str(gff_file), dbfn=str(self.db_path), force=True, keep_order=True, merge_strategy='merge', id_spec=['ID', 'Name'])
        else:
            logging.info(f"Loading existing GFF database from {self.db_path}...")
            self.db = gffutils.FeatureDB(str(self.db_path))

    def find_gene(self, gene_identifier: str) -> Optional[gffutils.Feature]:
        """Finds a gene by its ID or Name attribute."""
        try:
            # Try fetching by ID first
            gene = self.db[gene_identifier]
            if gene.featuretype == 'gene':
                return gene
        except gffutils.exceptions.FeatureNotFoundError:
            pass

        # If not found by ID, search by Name
        for gene in self.db.features_of_type('gene'):
            if gene.attributes.get('Name', [None])[0] == gene_identifier:
                return gene
        
        logging.error(f"Gene '{gene_identifier}' not found in GFF file.")
        return None

    def get_children(self, feature: gffutils.Feature) -> Iterator[gffutils.Feature]:
        """Yields all direct children of a feature, sorted by position."""
        return self.db.children(feature, order_by='start')


class SequenceExtractor:
    """Extracts feature sequences from a reference FASTA."""
    def __init__(self, reference_fasta: Path):
        self.fasta = pysam.FastaFile(str(reference_fasta))

    def get_sequence(self, feature: gffutils.Feature) -> Optional[str]:
        """Extracts a feature's sequence, handling strand."""
        try:
            seq = self.fasta.fetch(feature.chrom, feature.start - 1, feature.end)
        except (KeyError, ValueError) as e:
            logging.warning(f"Could not fetch sequence for {feature.id} ({feature.chrom}:{feature.start}-{feature.end}): {e}")
            return None

        if feature.strand == '-':
            return self._reverse_complement(seq)
        return seq

    @staticmethod
    def _reverse_complement(dna: str) -> str:
        """Returns the reverse complement of a DNA sequence."""
        complement = str.maketrans('ATCGNRY', 'TAGCNYR')
        return dna.upper().translate(complement)[::-1]


class HierarchicalAligner:
    """Performs hierarchical alignment of features against target sequences."""
    def __init__(self, gff_proc: GFFProcessor, seq_ext: SequenceExtractor):
        self.gff = gff_proc
        self.seq = seq_ext
        # Alignment parameters for pairwise2.align.localms:
        # (match, mismatch, open_gap, extend_gap)
        self.align_params = (5, -4, -7, -2)

    def align_gene_to_haplotypes(self, gene: gffutils.Feature, haplotypes_fasta: Path) -> Dict[str, AlignmentResult]:
        """Aligns a gene and its sub-features to all sequences in a multi-sample FASTA."""
        gene_seq = self.seq.get_sequence(gene)
        if not gene_seq:
            logging.error(f"Could not extract sequence for gene '{gene.id}'. Aborting.")
            return {}

        results = {}
        with pysam.FastaFile(str(haplotypes_fasta)) as multi_fasta:
            for hap_name in multi_fasta.references:
                logging.info(f"--- Aligning to haplotype: {hap_name} ---")
                hap_seq = multi_fasta.fetch(hap_name)
                
                # 1. Align the parent gene
                # Use local alignment to find the best match for the gene sequence in the haplotype
                alignments = pairwise2.align.localms(hap_seq, gene_seq, *self.align_params)
                
                if not alignments:
                    logging.warning(f"Gene '{gene.id}' did not align to haplotype '{hap_name}'.")
                    continue
                
                # Take the best alignment
                best_aln = alignments[0]
                aligned_hap_seq, aligned_gene_seq, score, target_start, target_end = best_aln

                # Calculate identity
                matches = sum(1 for a, b in zip(aligned_hap_seq, aligned_gene_seq) if a == b)
                align_len = len(aligned_hap_seq)
                identity = matches / align_len if align_len > 0 else 0

                gene_result = AlignmentResult(
                    feature_id=gene.id,
                    feature_type=gene.featuretype,
                    score=score,
                    target_start=target_start,
                    target_end=target_end,
                    identity=identity
                )
                
                # 2. Recursively align children within the parent's aligned region
                hap_region_seq = hap_seq[target_start:target_end]
                
                for child_feature in self.gff.get_children(gene):
                    child_result = self._align_feature_recursively(child_feature, hap_region_seq, offset=target_start)
                    if child_result:
                        gene_result.children.append(child_result)
                
                results[hap_name] = gene_result
        return results

    def _align_feature_recursively(self, feature: gffutils.Feature, target_region_seq: str, offset: int) -> Optional[AlignmentResult]:
        """Aligns a feature and its children within a constrained target region."""
        feature_seq = self.seq.get_sequence(feature)
        if not feature_seq:
            return None

        # Use local alignment to find the feature within the parent's region
        alignments = pairwise2.align.localms(target_region_seq, feature_seq, *self.align_params)

        if not alignments:
            logging.debug(f"Feature '{feature.id}' did not align within parent region.")
            return None

        best_aln = alignments[0]
        aligned_target_seq, aligned_feature_seq, score, start, end = best_aln

        identity = sum(1 for a, b in zip(aligned_target_seq, aligned_feature_seq) if a == b) / len(aligned_target_seq) if len(aligned_target_seq) > 0 else 0
        
        target_start = start
        target_end = end

        result = AlignmentResult(
            feature_id=feature.id,
            feature_type=feature.featuretype,
            score=score,
            target_start=target_start + offset,
            target_end=target_end + offset,
            identity=identity
        )
        
        # Recursively align children
        child_target_region = target_region_seq[target_start:target_end]
        for child_feature in self.gff.get_children(feature):
            child_result = self._align_feature_recursively(child_feature, child_target_region, offset=offset + target_start)
            if child_result:
                result.children.append(child_result)
        
        return result

# --- Output ---

def print_summary(results: Dict[str, AlignmentResult]):
    """Prints a formatted summary of the alignment results."""
    print("\n" + "="*80)
    print("HIERARCHICAL ALIGNMENT SUMMARY")
    print("="*80)

    for hap_name, gene_result in results.items():
        print(f"\n--- Haplotype: {hap_name} ---")
        _print_result_recursively(gene_result, indent=0)

def _print_result_recursively(result: AlignmentResult, indent: int):
    """Helper to print results hierarchically."""
    prefix = "  " * indent
    print(f"{prefix}- {result.feature_type} '{result.feature_id}':")
    print(f"{prefix}  - Aligned to region: {result.target_start}-{result.target_end}")
    print(f"{prefix}  - Score: {result.score}, Identity: {result.identity:.2%}")

    # Sort children by their start position for logical output
    sorted_children = sorted(result.children, key=lambda r: r.target_start)
    for child in sorted_children:
        _print_result_recursively(child, indent + 1)

# --- Main Execution ---

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Hierarchically align a gene's features against multiple sequences.")
    parser.add_argument("--multi-fasta", required=True, type=Path, help="Multi-sample FASTA file of haplotypes.")
    parser.add_argument("--gff", required=True, type=Path, help="GFF3 annotation file.")
    parser.add_argument("--reference", required=True, type=Path, help="Reference genome FASTA file.")
    parser.add_argument("--gene", required=True, type=str, help="Name or ID of the gene to align.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")
    args = parser.parse_args()

    setup_logging(args.verbose)

    try:
        # 1. Initialize tools
        gff_processor = GFFProcessor(args.gff)
        seq_extractor = SequenceExtractor(args.reference)
        aligner = HierarchicalAligner(gff_processor, seq_extractor)

        # 2. Find the gene
        gene_feature = gff_processor.find_gene(args.gene)
        if not gene_feature:
            sys.exit(1)
        
        logging.info(f"Found gene '{gene_feature.id}' ({gene_feature.attributes.get('Name', ['-'])[0]}) at {gene_feature.chrom}:{gene_feature.start}-{gene_feature.end}")

        # 3. Run alignment
        alignment_results = aligner.align_gene_to_haplotypes(gene_feature, args.multi_fasta)

        # 4. Print summary
        print_summary(alignment_results)

    except FileNotFoundError as e:
        logging.error(f"Error: Input file not found - {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}", exc_info=args.verbose)
        sys.exit(1)

if __name__ == "__main__":
    main()
