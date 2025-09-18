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

    def __init__(self, gff_parser: GFFParser, paf_parser: PAFParser):
        self.gff = gff_parser
        self.paf = paf_parser
        self.results: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(lambda: defaultdict(dict))

    def analyze(self):
        """Performs impact analysis for all features on all paths."""
        all_paths = self._get_all_paths()
        logging.info(f"Analyzing impacts across {len(all_paths)} paths...")

        for path_name in sorted(list(all_paths)):
            logging.debug(f"Processing path: {path_name}")
            for root_id in self.gff.roots:
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
            impact = self._aggregate_impacts(child_impacts, feature_id)

        self.results[path_name][feature_id] = impact
        return impact

    def _get_direct_impact(self, feature_id: str, path_name: str) -> Dict[str, Any]:
        """Determines the impact on a single leaf feature based on its alignment."""
        feature = self.gff.get_feature(feature_id)
        aln = self.paf.get_alignment(feature_id, path_name)

        if not aln or not aln.get("cigar"):
            return {"impact": LOSS_OF_FUNCTION, "details": ["Feature not aligned or alignment is partial"]}

        # Rule-based impact assessment using CIGAR string
        cigar = aln["cigar"]
        mismatches = len(re.findall(r'\d+X', cigar))
        insertions = sum(int(n) for n in re.findall(r'(\d+)I', cigar))
        deletions = sum(int(n) for n in re.findall(r'(\d+)D', cigar))
        
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

        return {"impact": impact, "details": details}

    def _aggregate_impacts(self, child_impacts: List[Dict[str, Any]], parent_id: str) -> Dict[str, Any]:
        """Aggregates impacts from children to determine parent's impact."""
        if not child_impacts:
            # If a parent has no children in the GFF, analyze it directly
            # This is a fallback, shouldn't happen for genes/transcripts
            return self._get_direct_impact(parent_id, "reference") # HACK: path_name is lost

        most_severe_impact = UNKNOWN
        all_details = []

        for impact_result in child_impacts:
            impact = impact_result["impact"]
            if IMPACT_SEVERITY[impact] > IMPACT_SEVERITY[most_severe_impact]:
                most_severe_impact = impact
            all_details.extend(impact_result["details"])
        
        # Consolidate details
        unique_details = sorted(list(set(all_details)))
        return {"impact": most_severe_impact, "details": unique_details}

    def print_summary(self):
        """Prints a hierarchical summary of the impact analysis."""
        print("\n" + "="*80)
        print("IMPACT ANALYSIS SUMMARY")
        print("="*80)

        for path_name in sorted(self.results.keys()):
            print(f"\n--- Path: {path_name} ---")
            for root_id in sorted(self.gff.roots, key=lambda r: self.gff.get_feature(r)['start']):
                self._print_feature_summary(root_id, path_name, indent=0)

    def _print_feature_summary(self, feature_id: str, path_name: str, indent: int):
        """Recursively prints the summary for a feature and its children."""
        feature = self.gff.get_feature(feature_id)
        result = self.results[path_name].get(feature_id)
        if not feature or not result:
            return

        prefix = "  " * indent
        gene_name = feature['attributes'].get('gene_name', feature_id)
        print(f"{prefix}- {feature['type']} '{gene_name}': {result['impact']}")
        for detail in result['details']:
            print(f"{prefix}  - {detail}")

        for child_id in sorted(self.gff.get_children(feature_id), key=lambda c: self.gff.get_feature(c)['start']):
            self._print_feature_summary(child_id, path_name, indent + 1)


def main():
    """Main function to run the impact analysis."""
    parser = argparse.ArgumentParser(description="Analyze functional impact of variants from a PAF file.")
    parser.add_argument("--paf", required=True, type=Path, help="Input PAF file from minimap2.")
    parser.add_argument("--gff", required=True, type=Path, help="GFF3 annotation file used for alignment.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")
    args = parser.parse_args()

    setup_logging(args.verbose)

    try:
        gff_parser = GFFParser(args.gff)
        paf_parser = PAFParser(args.paf)

        analyzer = ImpactAnalyzer(gff_parser, paf_parser)
        analyzer.analyze()
        analyzer.print_summary()

    except FileNotFoundError as e:
        logging.error(f"Error: Input file not found - {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}", exc_info=args.verbose)


if __name__ == "__main__":
    main()
