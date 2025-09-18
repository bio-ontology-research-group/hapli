#!/usr/bin/env python3

import argparse
import logging
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

import gffutils
import pysam
from tqdm import tqdm

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
    def __init__(self, feature_id: str, feature_type: str, mapq: int, target_start: int, target_end: int, identity: float):
        self.feature_id = feature_id
        self.feature_type = feature_type
        self.mapq = mapq
        self.target_start = target_start
        self.target_end = target_end
        self.identity = identity
        self.children: List[AlignmentResult] = []

    def __repr__(self) -> str:
        return f"AlignmentResult(id={self.feature_id}, type={self.feature_type}, mapq={self.mapq})"

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
            # Try fetching by ID first. This is the most reliable way.
            feature = self.db[gene_identifier]
            if feature.featuretype == 'gene':
                return feature
            # If an ID matches but it's not a gene, we'll proceed to search by Name.
        except gffutils.exceptions.FeatureNotFoundError:
            pass  # Not found by ID, will search by Name below.

        # If not found by ID, or if the ID matched a non-gene feature, search by Name.
        for gene in self.db.features_of_type('gene'):
            # The 'Name' attribute can be a list of names.
            if gene_identifier in gene.attributes.get('Name', []):
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
    """Performs hierarchical alignment of features against target sequences using minimap2."""
    def __init__(self, gff_proc: GFFProcessor, seq_ext: SequenceExtractor, threads: int = 1):
        self.gff = gff_proc
        self.seq = seq_ext
        self.threads = threads

    def _get_all_features_recursive(self, feature: gffutils.Feature, feature_list: List[gffutils.Feature]):
        """Recursively gather a feature and all its descendants."""
        feature_list.append(feature)
        for child in self.gff.get_children(feature):
            self._get_all_features_recursive(child, feature_list)

    def _run_minimap2(self, target_fasta: Path, query_fasta: Path) -> Dict[str, List[Dict[str, Any]]]:
        """Runs minimap2 and parses PAF output."""
        logging.info(f"Running minimap2 with {self.threads} threads...")
        cmd = [
            "minimap2",
            "-x", "asm20",  # Preset for accurate assembly-to-assembly alignment
            "-c",           # Output CIGAR string in PAF
            "--paf-no-hit", # Suppress output for sequences with no hits
            "-t", str(self.threads),
            str(target_fasta),
            str(query_fasta)
        ]
            
        try:
            process = subprocess.run(cmd, capture_output=True, text=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logging.error(f"minimap2 execution failed: {e}")
            if isinstance(e, subprocess.CalledProcessError):
                logging.error(f"minimap2 stderr:\n{e.stderr}")
            sys.exit(1)

        alignments = defaultdict(list)
        for line in process.stdout.strip().split('\n'):
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 12:
                continue
                
            paf_record = {
                "query_name": parts[0],
                "target_name": parts[5],
                "target_start": int(parts[7]),
                "target_end": int(parts[8]),
                "matches": int(parts[9]),
                "align_len": int(parts[10]),
                "mapq": int(parts[11]),
            }
            alignments[paf_record["target_name"]].append(paf_record)
        return alignments

    def align_gene_to_haplotypes(self, gene: gffutils.Feature, haplotypes_fasta: Path) -> Dict[str, AlignmentResult]:
        """Aligns a gene and its sub-features to all sequences in a multi-sample FASTA."""
        # 1. Gather all features and their sequences
        all_features = []
        self._get_all_features_recursive(gene, all_features)
            
        logging.info(f"Extracted {len(all_features)} features (gene and descendants) to align.")

        with tempfile.NamedTemporaryFile(mode='w+', delete=True, suffix=".fa") as query_fasta:
            for feature in tqdm(all_features, desc="Extracting feature sequences"):
                seq = self.seq.get_sequence(feature)
                if seq:
                    query_fasta.write(f">{feature.id}\n{seq}\n")
            query_fasta.flush()

            # 2. Run minimap2 to get all alignments at once
            all_alignments = self._run_minimap2(haplotypes_fasta, Path(query_fasta.name))

        # 3. Build hierarchical results from flat alignment list
        results = {}
        hap_names = list(all_alignments.keys())
        for hap_name in tqdm(hap_names, desc="Building hierarchical results"):
            # Group alignments by query name and select the best one (highest MAPQ)
            # This handles cases where a feature aligns to multiple places on a haplotype.
            grouped_by_query = defaultdict(list)
            for aln in all_alignments[hap_name]:
                grouped_by_query[aln['query_name']].append(aln)
            
            alignments_for_hap = {}
            for query_name, aln_list in grouped_by_query.items():
                if aln_list:
                    alignments_for_hap[query_name] = max(aln_list, key=lambda x: x['mapq'])

            # Find the top-level gene alignment
            gene_aln_data = alignments_for_hap.get(gene.id)
            if not gene_aln_data:
                continue

            identity = gene_aln_data['matches'] / gene_aln_data['align_len'] if gene_aln_data['align_len'] > 0 else 0
            gene_result = AlignmentResult(
                feature_id=gene.id,
                feature_type=gene.featuretype,
                mapq=gene_aln_data['mapq'],
                target_start=gene_aln_data['target_start'],
                target_end=gene_aln_data['target_end'],
                identity=identity
            )

            # Recursively build the tree for children
            self._build_result_tree(gene, gene_result, alignments_for_hap)
            results[hap_name] = gene_result
                
        return results

    def _build_result_tree(self, parent_feature: gffutils.Feature, parent_result: AlignmentResult, alignments: Dict[str, Dict]):
        """Recursively constructs the alignment result tree, applying hierarchical constraints."""
        for child_feature in self.gff.get_children(parent_feature):
            child_aln_data = alignments.get(child_feature.id)
            if not child_aln_data:
                continue

            # Hierarchical constraint: child must be within parent's aligned region
            if not (child_aln_data['target_start'] >= parent_result.target_start and child_aln_data['target_end'] <= parent_result.target_end):
                logging.debug(f"Skipping child {child_feature.id} as its alignment is outside parent {parent_feature.id}'s region.")
                continue

            identity = child_aln_data['matches'] / child_aln_data['align_len'] if child_aln_data['align_len'] > 0 else 0
            child_result = AlignmentResult(
                feature_id=child_feature.id,
                feature_type=child_feature.featuretype,
                mapq=child_aln_data['mapq'],
                target_start=child_aln_data['target_start'],
                target_end=child_aln_data['target_end'],
                identity=identity
            )
            parent_result.children.append(child_result)
                
            # Recurse
            self._build_result_tree(child_feature, child_result, alignments)

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
    print(f"{prefix}  - MAPQ: {result.mapq}, Identity: {result.identity:.2%}")

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
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for minimap2.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")
    args = parser.parse_args()

    setup_logging(args.verbose)

    try:
        # 1. Initialize tools
        gff_processor = GFFProcessor(args.gff)
        seq_extractor = SequenceExtractor(args.reference)
        aligner = HierarchicalAligner(gff_processor, seq_extractor, threads=args.threads)

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
        logging.error(f"An unexpected error occurred: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
