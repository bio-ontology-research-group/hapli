#!/usr/bin/env python3

import argparse
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pysam

# --- Impact Classification ---
LOSS_OF_FUNCTION = "Loss of Function"
CHANGE_OF_FUNCTION = "Change in Function"
NO_CHANGE = "No Change"
UNKNOWN = "Unknown"

# Severity order for impact aggregation
IMPACT_SEVERITY = {
    LOSS_OF_FUNCTION: 3,
    CHANGE_OF_FUNCTION: 2,
    NO_CHANGE: 1,
    UNKNOWN: 0,
}


def setup_logging(verbose: bool = False):
    """Set up logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="[%(asctime)s] %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S")


class GFFParser:
    """Parses a GFF file to build a feature hierarchy."""

    def __init__(self, gff_file: Path):
        self.gff_file = gff_file
        self.features: Dict[str, Dict[str, Any]] = {}
        self.children: Dict[str, List[str]] = defaultdict(list)
        self.roots: List[str] = []
        self._parse()

    def _parse(self):
        """Reads the GFF file and builds the feature tree."""
        logging.info(f"Parsing GFF hierarchy from {self.gff_file}...")
        # First pass: collect all features
        with open(self.gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                attributes = {k: v for k, v in (re.split(r'=', attr, 1) for attr in parts[8].split(';'))}
                feature_id = attributes.get("ID")
                if feature_id:
                    self.features[feature_id] = {
                        "id": feature_id,
                        "type": parts[2],
                        "chrom": parts[0],
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[6],
                        "parent": attributes.get("Parent"),
                        "attributes": attributes,
                    }

        # Second pass: build hierarchy
        for feature_id, feature_data in self.features.items():
            parent_id = feature_data.get("parent")
            if parent_id:
                self.children[parent_id].append(feature_id)
            else:
                self.roots.append(feature_id)
        logging.info(f"Found {len(self.roots)} root features and {len(self.features)} total features.")

    def get_feature(self, feature_id: str) -> Optional[Dict[str, Any]]:
        return self.features.get(feature_id)

    def get_children(self, feature_id: str) -> List[str]:
        return self.children.get(feature_id, [])


class PAFParser:
    """Parses a PAF file to store alignment information."""

    def __init__(self, paf_file: Path):
        self.paf_file = paf_file
        # {query_name: {target_name: [alignment_dict]}}
        self.alignments: Dict[str, Dict[str, List[Dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
        self._parse()

    def _parse(self):
        """Reads the PAF file and stores alignments."""
        logging.info(f"Parsing PAF alignments from {self.paf_file}...")
        with open(self.paf_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                
                alignment = {
                    "query_name": parts[0],
                    "query_len": int(parts[1]),
                    "query_start": int(parts[2]),
                    "query_end": int(parts[3]),
                    "strand": parts[4],
                    "target_name": parts[5],
                    "target_len": int(parts[6]),
                    "target_start": int(parts[7]),
                    "target_end": int(parts[8]),
                    "matches": int(parts[9]),
                    "align_len": int(parts[10]),
                    "mapq": int(parts[11]),
                    "cigar": None,
                }

                # Find and store CIGAR string
                for tag in parts[12:]:
                    if tag.startswith("cg:Z:"):
                        alignment["cigar"] = tag[5:]
                        break
                
                self.alignments[alignment["query_name"]][alignment["target_name"]].append(alignment)
        logging.info(f"Loaded alignments for {len(self.alignments)} unique features.")

    def get_alignment(self, query_name: str, target_name: str) -> Optional[Dict[str, Any]]:
        """Returns the best alignment for a given query and target."""
        alns = self.alignments.get(query_name, {}).get(target_name)
        if not alns:
            return None
        # Return the one with the highest mapping quality
        return max(alns, key=lambda x: x['mapq'])


class ImpactAnalyzer:
    """Analyzes feature alignments to determine functional impact."""

    def __init__(self, gff_parser: GFFParser, paf_parser: PAFParser, reference_fasta: Path, path_fasta: Path, target_genes: List[str], use_color: bool, show_alignments: bool):
        self.gff = gff_parser
        self.paf = paf_parser
        self.target_genes = set(target_genes)
        self.use_color = use_color
        self.show_alignments = show_alignments
        self.results: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(lambda: defaultdict(dict))

        logging.info("Loading reference genome sequence...")
        self.reference_genome = pysam.FastaFile(str(reference_fasta))
        logging.info("Loading path sequences...")
        self.path_seqs = pysam.FastaFile(str(path_fasta))

        # ANSI colors
        self.colors = {
            "RED": '\033[91m' if use_color else '',
            "GREEN": '\033[92m' if use_color else '',
            "YELLOW": '\033[93m' if use_color else '',
            "BLUE": '\033[94m' if use_color else '',
            "BOLD": '\033[1m' if use_color else '',
            "ENDC": '\033[0m' if use_color else ''
        }

    def analyze(self):
        """Performs impact analysis for all features on all paths."""
        all_paths = self._get_all_paths()
        logging.info(f"Analyzing impacts across {len(all_paths)} paths...")

        root_features = self.gff.roots
        if self.target_genes:
            logging.info(f"Filtering for genes: {', '.join(self.target_genes)}")
            root_features = [
                r for r in self.gff.roots
                if self.gff.get_feature(r)['attributes'].get('gene_name') in self.target_genes
            ]
            logging.info(f"Found {len(root_features)} matching root features to analyze.")

        for path_name in sorted(list(all_paths)):
            logging.debug(f"Processing path: {path_name}")
            for root_id in root_features:
                # Only analyze features that are actually in the PAF file
                if root_id in self.paf.alignments:
                    self._analyze_feature_recursively(root_id, path_name)
        
        return self.results

    def _get_all_paths(self) -> Set[str]:
        """Collect all unique target path names from the PAF file."""
        paths = set()
        for query_alns in self.paf.alignments.values():
            paths.update(query_alns.keys())
        return paths

    def _analyze_feature_recursively(self, feature_id: str, path_name: str) -> Dict[str, Any]:
        """Analyzes a feature and its children, returning the aggregated impact."""
        if feature_id in self.results[path_name]:
            return self.results[path_name][feature_id]

        child_ids = self.gff.get_children(feature_id)
        
        if not child_ids:
            # This is a leaf node (e.g., exon), perform direct analysis
            impact = self._get_direct_impact(feature_id, path_name)
        else:
            # This is a parent node (e.g., gene), aggregate from children
            child_impacts = [self._analyze_feature_recursively(child_id, path_name) for child_id in child_ids]
            impact = self._aggregate_impacts(child_impacts, feature_id, path_name)

        self.results[path_name][feature_id] = impact
        return impact

    def _get_direct_impact(self, feature_id: str, path_name: str) -> Dict[str, Any]:
        """Determines the impact on a single leaf feature based on its alignment."""
        feature = self.gff.get_feature(feature_id)
        aln = self.paf.get_alignment(feature_id, path_name)

        if not aln or not aln.get("cigar"):
            return {"impact": LOSS_OF_FUNCTION, "details": ["Feature not aligned or alignment is partial"], "identity": 0.0, "mismatches": 0, "insertions": 0, "deletions": 0, "alignment": None}

        identity = (aln['matches'] / aln['align_len'] * 100) if aln['align_len'] > 0 else 0

        cigar_tuples = re.findall(r'(\d+)([MIDNSHPX=])', aln["cigar"])
        insertions = sum(int(length) for length, op in cigar_tuples if op == 'I')
        deletions = sum(int(length) for length, op in cigar_tuples if op == 'D')
        
        # Approximate mismatches (substitutions) from PAF info. For exact counts, `cs` tag is needed.
        mismatches = 0
        if aln['align_len'] > 0:
            mismatches = aln['align_len'] - aln['matches'] - insertions - deletions
            mismatches = max(0, mismatches)

        details = []
        impact = NO_CHANGE

        if feature["type"] == "CDS":
            total_indels = insertions + deletions
            if total_indels > 0 and total_indels % 3 != 0:
                impact = LOSS_OF_FUNCTION
                details.append(f"Frameshift indel (length {total_indels})")
            elif total_indels > 0:
                impact = CHANGE_OF_FUNCTION
                details.append(f"In-frame indel (length {total_indels})")
            
            if mismatches > 0:
                if impact == NO_CHANGE:
                    impact = CHANGE_OF_FUNCTION
                details.append(f"{mismatches} mismatch(es) found (potential non-synonymous change)")

        else:  # Non-coding features (UTR, intron, etc.)
            if deletions > 20 or insertions > 20:
                impact = CHANGE_OF_FUNCTION
                details.append(f"Large indel found (I:{insertions}, D:{deletions})")
            elif mismatches > 5:
                impact = CHANGE_OF_FUNCTION
                details.append(f"{mismatches} mismatch(es) found")

        if not details:
            details.append("No significant changes detected")

        return {
            "impact": impact,
            "details": details,
            "identity": identity,
            "mismatches": mismatches,
            "insertions": insertions,
            "deletions": deletions,
            "alignment": aln,
        }

    def _aggregate_impacts(self, child_impacts: List[Dict[str, Any]], parent_id: str, path_name: str) -> Dict[str, Any]:
        """Aggregates impacts from children to determine parent's impact."""
        if not child_impacts:
            # This is a parent with no children in the GFF (e.g. a gene with no transcripts). Treat as a leaf.
            return self._get_direct_impact(parent_id, path_name)

        most_severe_impact = UNKNOWN
        all_details = []
        best_child_result = None

        for impact_result in child_impacts:
            if not impact_result:
                continue
            impact = impact_result["impact"]
            if IMPACT_SEVERITY[impact] > IMPACT_SEVERITY[most_severe_impact]:
                most_severe_impact = impact
                best_child_result = impact_result
            all_details.extend(impact_result["details"])
        
        if best_child_result:
            unique_details = sorted(list(set(all_details)))
            return {
                "impact": most_severe_impact,
                "details": unique_details,
                "identity": best_child_result.get('identity', 0),
                "mismatches": best_child_result.get('mismatches', 0),
                "insertions": best_child_result.get('insertions', 0),
                "deletions": best_child_result.get('deletions', 0),
                "alignment": None  # No single alignment for parent
            }
        else:
            return {"impact": UNKNOWN, "details": ["No child impacts found"], "identity": 0.0, "mismatches": 0, "insertions": 0, "deletions": 0, "alignment": None}

    @staticmethod
    def _reverse_complement(dna: str) -> str:
        """Returns the reverse complement of a DNA sequence."""
        complement = str.maketrans('ATCGNRY', 'TAGCNYR')
        return dna.upper().translate(complement)[::-1]

    def get_feature_sequence(self, feature_id: str) -> Optional[str]:
        """Extracts a feature's sequence from the reference genome."""
        feature = self.gff.get_feature(feature_id)
        if not feature:
            return None
        
        try:
            # GFF is 1-based, pysam is 0-based
            seq = self.reference_genome.fetch(feature['chrom'], feature['start'] - 1, feature['end'])
        except (KeyError, ValueError) as e:
            logging.warning(f"Could not fetch sequence for {feature_id} ({feature['chrom']}:{feature['start']}-{feature['end']}): {e}")
            return None

        if feature['strand'] == '-':
            return self._reverse_complement(seq)
        return seq

    def _get_alignment_visualization(self, feature_id: str, aln: Dict[str, Any]) -> str:
        """Generates a colorized string showing the base-level alignment."""
        query_full_seq = self.get_feature_sequence(aln['query_name'])
        if not query_full_seq:
            return f"  (Sequence for feature '{aln['query_name']}' could not be extracted, cannot display alignment)"

        try:
            target_full_seq = self.path_seqs.fetch(aln['target_name'])
        except (KeyError, ValueError) as e:
            return f"  (Sequence not found for path '{aln['target_name']}', cannot display alignment)"

        cigar_tuples = re.findall(r'(\d+)([MIDNSHPX=])', aln["cigar"])
        
        target_pos = aln['target_start']
        query_pos = aln['query_start']
        
        if cigar_tuples and cigar_tuples[0][1] == 'S':
            query_pos += int(cigar_tuples[0][0])
            cigar_tuples.pop(0)

        viz = {'target': '', 'mid': '', 'query': ''}
        
        for length_str, op in cigar_tuples:
            length = int(length_str)
            if op == 'M':
                for _ in range(length):
                    t = target_full_seq[target_pos]
                    q = query_full_seq[query_pos]
                    viz['target'] += t
                    viz['query'] += q
                    if t.upper() == q.upper():
                        viz['mid'] += '|'
                    else:
                        viz['mid'] += ' '
                        viz['target'] = viz['target'][:-1] + self.colors['RED'] + t + self.colors['ENDC']
                        viz['query'] = viz['query'][:-1] + self.colors['RED'] + q + self.colors['ENDC']
                    target_pos += 1
                    query_pos += 1
            elif op == 'I':
                for _ in range(length):
                    q = query_full_seq[query_pos]
                    viz['target'] += self.colors['YELLOW'] + '-' + self.colors['ENDC']
                    viz['mid'] += ' '
                    viz['query'] += self.colors['YELLOW'] + q + self.colors['ENDC']
                    query_pos += 1
            elif op == 'D':
                for _ in range(length):
                    t = target_full_seq[target_pos]
                    viz['target'] += self.colors['YELLOW'] + t + self.colors['ENDC']
                    viz['mid'] += ' '
                    viz['query'] += self.colors['YELLOW'] + '-' + self.colors['ENDC']
                    target_pos += 1
            elif op in ['S', 'H']:
                pass

        lines = []
        chunk_size = 80
        for i in range(0, len(viz['target']), chunk_size):
            lines.append(f"  Target: {viz['target'][i:i+chunk_size]}")
            lines.append(f"          {viz['mid'][i:i+chunk_size]}")
            lines.append(f"  Query:  {viz['query'][i:i+chunk_size]}")
            lines.append("")
            
        return "\n".join(lines)

    def print_summary(self):
        """Prints a hierarchical summary of the impact analysis."""
        print("\n" + "="*80)
        print("IMPACT ANALYSIS SUMMARY")
        print("="*80)

        sorted_paths = sorted(self.results.keys())
        if not sorted_paths:
            print("No alignments found to analyze.")
            return

        for path_name in sorted_paths:
            print(f"\n--- Path: {self.colors['BLUE']}{path_name}{self.colors['ENDC']} ---")
            
            root_features = self.gff.roots
            if self.target_genes:
                root_features = [r for r in self.gff.roots if self.gff.get_feature(r)['attributes'].get('gene_name') in self.target_genes]

            sorted_roots = sorted(root_features, key=lambda r: self.gff.get_feature(r)['start'])

            for root_id in sorted_roots:
                if root_id in self.results[path_name]:
                    self._print_feature_summary(root_id, path_name, indent=0)

    def _print_feature_summary(self, feature_id: str, path_name: str, indent: int):
        """Recursively prints the summary for a feature and its children."""
        feature = self.gff.get_feature(feature_id)
        result = self.results[path_name].get(feature_id)
        if not feature or not result:
            return

        prefix = "  " * indent
        
        display_name = feature['attributes'].get('gene_name', feature_id) if feature['type'] == 'gene' else feature_id

        impact_color = {
            LOSS_OF_FUNCTION: self.colors['RED'],
            CHANGE_OF_FUNCTION: self.colors['YELLOW'],
            NO_CHANGE: self.colors['GREEN'],
        }.get(result['impact'], '')

        print(f"{prefix}- {self.colors['BOLD']}{feature['type']}{self.colors['ENDC']} '{display_name}': {impact_color}{result['impact']}{self.colors['ENDC']}")

        identity = result.get('identity', 0.0)
        mismatches = result.get('mismatches', 0)
        insertions = result.get('insertions', 0)
        deletions = result.get('deletions', 0)

        print(f"{prefix}  - Similarity: {identity:.2f}% | Changes: {mismatches} mismatches, {insertions} insertions, {deletions} deletions")

        for detail in result['details']:
            print(f"{prefix}  - {detail}")

        if self.show_alignments and result.get('alignment'):
            print(f"{prefix}  - Alignment:")
            viz_str = self._get_alignment_visualization(feature_id, result['alignment'])
            indented_viz = "\n".join([f"{prefix}  {line}" for line in viz_str.splitlines()])
            print(indented_viz)

        for child_id in sorted(self.gff.get_children(feature_id), key=lambda c: self.gff.get_feature(c)['start']):
            if child_id in self.results[path_name]:
                self._print_feature_summary(child_id, path_name, indent + 1)


def main():
    """Main function to run the impact analysis."""
    parser = argparse.ArgumentParser(description="Analyze functional impact of variants from a PAF file.")
    parser.add_argument("--paf", required=True, type=Path, help="Input PAF file from minimap2.")
    parser.add_argument("--gff", required=True, type=Path, help="GFF3 annotation file used for alignment.")
    parser.add_argument("--reference-fasta", required=True, type=Path, help="Reference genome FASTA file.")
    parser.add_argument("--path-fasta", required=True, type=Path, help="FASTA file for path sequences (haplotypes).")
    parser.add_argument("--genes", type=str, help="Comma-separated list of gene names to analyze.")
    parser.add_argument("--show-alignments", action="store_true", help="Show detailed base-level alignments.")
    parser.add_argument("--no-color", action="store_true", help="Disable color output.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")
    args = parser.parse_args()

    setup_logging(args.verbose)

    try:
        target_genes = [g.strip() for g in args.genes.split(',')] if args.genes else []
        
        gff_parser = GFFParser(args.gff)
        paf_parser = PAFParser(args.paf)

        analyzer = ImpactAnalyzer(
            gff_parser,
            paf_parser,
            args.reference_fasta,
            args.path_fasta,
            target_genes,
            use_color=not args.no_color,
            show_alignments=args.show_alignments
        )
        analyzer.analyze()
        analyzer.print_summary()

    except FileNotFoundError as e:
        logging.error(f"Error: Input file not found - {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}", exc_info=args.verbose)


if __name__ == "__main__":
    main()
