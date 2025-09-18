#!/usr/bin/env python3

import argparse
import json
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import gffutils
import pysam
from tqdm import tqdm

# --- Impact Classification ---
FUNCTIONAL = "Functional"
NON_FUNCTIONAL = "Non-functional"
CHANGE_OF_FUNCTION = "Change in Function"
PARTIAL_FUNCTION = "Partial Function"
DUPLICATED = "Duplicated"
TRUNCATED = "Truncated"
FRAMESHIFT = "Frameshift"
NO_CHANGE = "No Change"
NOT_FOUND = "Not Found"
UNKNOWN = "Unknown"

# Severity order for impact aggregation
IMPACT_SEVERITY = {
    NON_FUNCTIONAL: 5,
    FRAMESHIFT: 4,
    TRUNCATED: 4,
    CHANGE_OF_FUNCTION: 3,
    PARTIAL_FUNCTION: 2,
    DUPLICATED: 2,
    FUNCTIONAL: 1,
    NO_CHANGE: 1,
    NOT_FOUND: 0,
    UNKNOWN: 0,
}

# Color codes for terminal output
COLORS = {
    "RED": '\033[91m',
    "GREEN": '\033[92m',
    "YELLOW": '\033[93m',
    "BLUE": '\033[94m',
    "MAGENTA": '\033[95m',
    "CYAN": '\033[96m',
    "BOLD": '\033[1m',
    "UNDERLINE": '\033[4m',
    "ENDC": '\033[0m'
}


def setup_logging(verbose: bool = False):
    """Set up logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level, 
        format="[%(asctime)s] %(levelname)s: %(message)s", 
        datefmt="%Y-%m-%d %H:%M:%S"
    )


class GFFParser:
    """Parses a GFF file using gffutils to build a feature database."""

    def __init__(self, gff_file: Path):
        self.gff_file = gff_file
        self.db_path = gff_file.with_suffix(gff_file.suffix + '.db')
        
        logging.info(f"Creating/loading GFF database from {gff_file}...")
        try:
            # Attempt to load existing DB
            self.db = gffutils.FeatureDB(str(self.db_path))
            logging.info(f"Loaded existing GFF database: {self.db_path}")
        except Exception:
            logging.warning(f"Could not load GFF database '{self.db_path}'. Recreating it.")
            try:
                self.db = gffutils.create_db(
                    str(gff_file),
                    dbfn=str(self.db_path),
                    force=True,
                    keep_order=True,
                    merge_strategy='merge',
                    sort_attributes=True,
                    disable_infer_genes=True,
                    disable_infer_transcripts=True
                )
            except Exception as e:
                logging.error(f"Failed to create gffutils database: {e}")
                # Clean up failed DB file
                if self.db_path.exists():
                    self.db_path.unlink()
                sys.exit(1)

        self.roots = [f for f in self.db.features_of_type('gene')]
        logging.info(f"Loaded {len(list(self.db.all_features()))} features with {len(self.roots)} root gene features")

    def get_feature(self, feature_id: str) -> Optional[gffutils.Feature]:
        """Retrieves a feature by its ID."""
        try:
            return self.db[feature_id]
        except gffutils.exceptions.FeatureNotFoundError:
            return None

    def get_children(self, feature_id: str) -> List[gffutils.Feature]:
        """Retrieves all direct children of a feature."""
        try:
            # Sorting by start position for consistent order
            return list(self.db.children(feature_id, order_by='start'))
        except gffutils.exceptions.FeatureNotFoundError:
            return []


class PAFParser:
    """Parses a PAF file to store alignment information."""

    def __init__(self, paf_file: Path):
        self.paf_file = paf_file
        # {query_name: {target_name: [alignment_dict]}}
        self.alignments: Dict[str, Dict[str, List[Dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
        self.target_paths: Set[str] = set()
        self._parse()

    def _parse(self):
        """Reads the PAF file and stores alignments."""
        logging.info(f"Parsing PAF alignments from {self.paf_file}...")
        
        # Count lines for progress
        total_lines = sum(1 for line in open(self.paf_file, 'r') if not line.startswith('#'))
        
        with open(self.paf_file, 'r') as f:
            with tqdm(total=total_lines, desc="Loading PAF alignments", unit=" alignments") as pbar:
                for line in f:
                    if line.startswith('#'):
                        continue
                    pbar.update(1)
                    
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
                        "coverage": (int(parts[3]) - int(parts[2])) / int(parts[1]) * 100 if int(parts[1]) > 0 else 0
                    }

                    # Find and store CIGAR string and other tags
                    for tag in parts[12:]:
                        if tag.startswith("cg:Z:"):
                            alignment["cigar"] = tag[5:]
                        elif tag.startswith("ft:Z:"):
                            alignment["feature_type"] = tag[5:]
                        elif tag.startswith("gn:Z:"):
                            alignment["gene_name"] = tag[5:]
                    
                    self.alignments[alignment["query_name"]][alignment["target_name"]].append(alignment)
                    self.target_paths.add(alignment["target_name"])
        
        logging.info(f"Loaded alignments for {len(self.alignments)} features across {len(self.target_paths)} paths")

    def get_alignments(self, query_name: str, target_name: str) -> List[Dict[str, Any]]:
        """Returns all alignments for a given query and target."""
        return self.alignments.get(query_name, {}).get(target_name, [])

    def get_best_alignment(self, query_name: str, target_name: str) -> Optional[Dict[str, Any]]:
        """Returns the best alignment for a given query and target."""
        alns = self.get_alignments(query_name, target_name)
        if not alns:
            return None
        # Return the one with the highest mapping quality and coverage
        return max(alns, key=lambda x: (x['mapq'], x['coverage']))

    def count_alignments_per_target(self, query_name: str, target_name: str) -> int:
        """Count how many times a feature aligns to a specific target."""
        return len(self.alignments.get(query_name, {}).get(target_name, []))


class ImpactAnalyzer:
    """Analyzes feature alignments to determine functional impact."""

    def __init__(self, gff_parser: GFFParser, paf_parser: PAFParser, reference_fasta: Path, 
                 path_fasta: Path, target_genes: List[str], use_color: bool, show_alignments: bool):
        self.gff = gff_parser
        self.paf = paf_parser
        self.target_genes = set(target_genes) if target_genes else set()
        self.use_color = use_color
        self.show_alignments = show_alignments
        self.results: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(lambda: defaultdict(dict))

        logging.info("Loading reference genome sequence...")
        self.reference_genome = pysam.FastaFile(str(reference_fasta))
        logging.info("Loading path sequences...")
        self.path_seqs = pysam.FastaFile(str(path_fasta))

        # Set up colors
        self.colors = COLORS if use_color else {k: '' for k in COLORS}

    def analyze(self):
        """Performs impact analysis for all features on all paths."""
        all_paths = sorted(list(self.paf.target_paths))
        logging.info(f"Analyzing impacts across {len(all_paths)} paths...")

        # Filter root features if target genes specified
        root_features = self.gff.roots
        if self.target_genes:
            logging.info(f"Filtering for genes: {', '.join(self.target_genes)}")
            filtered_roots = []
            for feature in self.gff.roots:
                gene_name = feature.attributes.get('gene_name', feature.attributes.get('Name', ['']))[0]
                if gene_name in self.target_genes or feature.id in self.target_genes:
                    filtered_roots.append(feature)
            root_features = filtered_roots
            logging.info(f"Found {len(root_features)} matching root features to analyze")

        # Analyze each path
        for path_name in tqdm(all_paths, desc="Analyzing paths", unit=" paths"):
            for root_feature in root_features:
                self._analyze_feature_recursively(root_feature.id, path_name)
        
        return self.results

    def _analyze_feature_recursively(self, feature_id: str, path_name: str) -> Dict[str, Any]:
        """Analyzes a feature and its children, returning the aggregated impact."""
        if feature_id in self.results[path_name]:
            return self.results[path_name][feature_id]

        feature = self.gff.get_feature(feature_id)
        if not feature:
            return {"impact": UNKNOWN, "details": ["Feature not found in GFF"], "identity": 0.0}

        child_features = self.gff.get_children(feature_id)
        
        if not child_features:
            # This is a leaf node, perform direct analysis
            impact = self._get_direct_impact(feature_id, path_name)
        else:
            # This is a parent node, aggregate from children
            child_impacts = []
            for child_feature in child_features:
                child_impact = self._analyze_feature_recursively(child_feature.id, path_name)
                if child_impact:
                    child_impacts.append(child_impact)
            impact = self._aggregate_impacts(child_impacts, feature_id, path_name)

        self.results[path_name][feature_id] = impact
        return impact

    def _get_direct_impact(self, feature_id: str, path_name: str) -> Dict[str, Any]:
        """Determines the impact on a single leaf feature based on its alignment."""
        feature = self.gff.get_feature(feature_id)
        alignments = self.paf.get_alignments(feature_id, path_name)
        
        # Count duplications for this specific target
        duplication_count = self.paf.count_alignments_per_target(feature_id, path_name)
        
        if not alignments:
            return {
                "impact": NOT_FOUND,
                "details": ["Feature not aligned to this path"],
                "identity": 0.0,
                "coverage": 0.0,
                "duplications": 0,
                "interpretation": "Feature is completely missing from this haplotype"
            }

        # Get best alignment for primary analysis
        best_aln = self.paf.get_best_alignment(feature_id, path_name)
        
        identity = (best_aln['matches'] / best_aln['align_len'] * 100) if best_aln['align_len'] > 0 else 0
        coverage = best_aln['coverage']
        
        # Parse CIGAR for detailed analysis
        insertions = 0
        deletions = 0
        mismatches = 0
        
        if best_aln.get("cigar"):
            cigar_tuples = re.findall(r'(\d+)([MIDNSHPX=])', best_aln["cigar"])
            insertions = sum(int(length) for length, op in cigar_tuples if op == 'I')
            deletions = sum(int(length) for length, op in cigar_tuples if op == 'D')
            
            # Estimate mismatches
            if best_aln['align_len'] > 0:
                mismatches = best_aln['align_len'] - best_aln['matches'] - insertions - deletions
                mismatches = max(0, mismatches)

        details = []
        impact = FUNCTIONAL
        interpretation = ""

        # Special handling for start/stop codons
        if feature.featuretype in ["start_codon", "stop_codon"]:
            if coverage < 100 or identity < 100:
                impact = NON_FUNCTIONAL
                details.append(f"{feature.featuretype} is altered")
                interpretation = f"{feature.featuretype} mutation will affect protein"
            else:
                impact = FUNCTIONAL
                details.append(f"{feature.featuretype} is intact")
                interpretation = f"{feature.featuretype} is preserved"
        # Analyze based on feature type
        elif feature.featuretype == "CDS":
            if coverage < 50:
                impact = TRUNCATED
                details.append(f"Only {coverage:.1f}% of CDS is aligned")
                interpretation = f"CDS is severely truncated, likely non-functional"
            elif coverage < 90:
                impact = PARTIAL_FUNCTION
                details.append(f"CDS coverage is {coverage:.1f}%")
                interpretation = f"CDS is partially present, may have reduced function"
            else:
                # Check for frameshifts
                total_indels = insertions + deletions
                if total_indels > 0 and total_indels % 3 != 0:
                    impact = FRAMESHIFT
                    details.append(f"Frameshift mutation (net indel: {total_indels}bp)")
                    interpretation = f"Frameshift will likely produce non-functional protein"
                elif total_indels > 0:
                    impact = CHANGE_OF_FUNCTION
                    details.append(f"In-frame indel ({total_indels}bp)")
                    interpretation = f"In-frame indel may alter protein function"
                elif mismatches > 0:
                    if identity < 95:
                        impact = CHANGE_OF_FUNCTION
                        details.append(f"{mismatches} substitutions ({identity:.1f}% identity)")
                        interpretation = f"Multiple amino acid changes likely alter protein function"
                    else:
                        impact = FUNCTIONAL
                        details.append(f"{mismatches} substitutions ({identity:.1f}% identity)")
                        interpretation = f"Minor changes, likely retains function"
                else:
                    impact = FUNCTIONAL
                    details.append(f"Perfect alignment ({identity:.1f}% identity)")
                    interpretation = f"CDS is intact and functional"

        elif feature.featuretype in ["five_prime_UTR", "three_prime_UTR"]:
            if coverage < 50:
                impact = PARTIAL_FUNCTION
                details.append(f"UTR coverage is only {coverage:.1f}%")
                interpretation = f"UTR is partially present, may affect regulation"
            elif identity < 90:
                impact = CHANGE_OF_FUNCTION
                details.append(f"UTR has {identity:.1f}% identity")
                interpretation = f"UTR changes may affect gene regulation"
            else:
                impact = FUNCTIONAL
                details.append(f"UTR is well preserved ({identity:.1f}% identity)")
                interpretation = f"UTR is intact, regulation likely preserved"

        elif feature.featuretype == "exon":
            if coverage < 80:
                impact = TRUNCATED
                details.append(f"Exon coverage is only {coverage:.1f}%")
                interpretation = f"Exon is significantly truncated"
            elif identity < 95:
                impact = CHANGE_OF_FUNCTION
                details.append(f"Exon has {identity:.1f}% identity")
                interpretation = f"Exon has significant changes"
            else:
                impact = FUNCTIONAL
                details.append(f"Exon is well preserved")
                interpretation = f"Exon is intact"

        else:  # Other features (intron, gene, transcript, etc.)
            if coverage < 50:
                impact = PARTIAL_FUNCTION
                details.append(f"Feature coverage is {coverage:.1f}%")
                interpretation = f"Feature is partially present"
            elif identity < 85:
                impact = CHANGE_OF_FUNCTION
                details.append(f"Feature has {identity:.1f}% identity")
                interpretation = f"Feature has moderate changes"
            else:
                impact = FUNCTIONAL
                details.append(f"Feature is preserved")
                interpretation = f"Feature is intact"

        # Check for duplications (but not for reference path)
        # Reference path shouldn't show duplications
        is_reference = "grch38" in path_name.lower() or "hg38" in path_name.lower() or "reference" in path_name.lower()
        if duplication_count > 1 and not is_reference:
            if impact == FUNCTIONAL:
                impact = DUPLICATED
            details.append(f"Feature appears {duplication_count} times on this path")
            interpretation += f" (duplicated {duplication_count}x)"

        return {
            "impact": impact,
            "details": details,
            "identity": identity,
            "coverage": coverage,
            "mismatches": mismatches,
            "insertions": insertions,
            "deletions": deletions,
            "duplications": duplication_count if not is_reference else 1,
            "alignment": best_aln,
            "all_alignments": alignments,
            "interpretation": interpretation
        }

    def _aggregate_impacts(self, child_impacts: List[Dict[str, Any]], parent_id: str, path_name: str) -> Dict[str, Any]:
        """Aggregates impacts from children to determine parent's impact."""
        if not child_impacts:
            return self._get_direct_impact(parent_id, path_name)

        feature = self.gff.get_feature(parent_id)
        
        # Count impact types
        impact_counts = defaultdict(int)
        for child in child_impacts:
            impact_counts[child["impact"]] += 1
        
        # Determine overall impact - STRICT aggregation
        # If ANY critical child is non-functional, parent is non-functional
        most_severe_impact = FUNCTIONAL
        
        # Check for critical failures
        critical_child_types = {"CDS", "start_codon", "stop_codon", "exon"}
        critical_failures = 0
        total_critical = 0
        
        for child_impact in child_impacts:
            child_id = child_impact.get("feature_id", "")
            if child_id:
                child_feature = self.gff.get_feature(child_id)
                if child_feature and child_feature.featuretype in critical_child_types:
                    total_critical += 1
                    if child_impact["impact"] in [NON_FUNCTIONAL, NOT_FOUND, FRAMESHIFT, TRUNCATED]:
                        critical_failures += 1
        
        # Strict rules for parent impact based on children
        if feature.featuretype == "transcript":
            # Transcript needs ALL critical components
            if impact_counts[NOT_FOUND] > 0 or impact_counts[NON_FUNCTIONAL] > 0:
                most_severe_impact = NON_FUNCTIONAL
            elif impact_counts[FRAMESHIFT] > 0:
                most_severe_impact = FRAMESHIFT
            elif impact_counts[TRUNCATED] > 0:
                most_severe_impact = TRUNCATED
            elif critical_failures > 0:
                most_severe_impact = NON_FUNCTIONAL
            elif impact_counts[CHANGE_OF_FUNCTION] > 0:
                most_severe_impact = CHANGE_OF_FUNCTION
            elif impact_counts[PARTIAL_FUNCTION] > 0:
                most_severe_impact = PARTIAL_FUNCTION
            else:
                most_severe_impact = FUNCTIONAL
        else:
            # For other parent types, use severity-based aggregation
            for impact in child_impacts:
                if IMPACT_SEVERITY[impact["impact"]] > IMPACT_SEVERITY[most_severe_impact]:
                    most_severe_impact = impact["impact"]
        
        # Calculate aggregate statistics
        avg_identity = sum(c.get('identity', 0) for c in child_impacts) / len(child_impacts) if child_impacts else 0
        avg_coverage = sum(c.get('coverage', 0) for c in child_impacts) / len(child_impacts) if child_impacts else 0
        total_duplications = sum(c.get('duplications', 0) for c in child_impacts)
        
        # Build interpretation
        interpretation = f"{feature.featuretype} has "
        if impact_counts[NOT_FOUND] > 0:
            interpretation += f"{impact_counts[NOT_FOUND]} missing components, "
        if impact_counts[NON_FUNCTIONAL] > 0:
            interpretation += f"{impact_counts[NON_FUNCTIONAL]} non-functional components, "
        if impact_counts[FRAMESHIFT] > 0:
            interpretation += f"{impact_counts[FRAMESHIFT]} frameshift mutations, "
        if impact_counts[TRUNCATED] > 0:
            interpretation += f"{impact_counts[TRUNCATED]} truncated components, "
        if impact_counts[CHANGE_OF_FUNCTION] > 0:
            interpretation += f"{impact_counts[CHANGE_OF_FUNCTION]} altered components, "
        if impact_counts[FUNCTIONAL] > 0:
            interpretation += f"{impact_counts[FUNCTIONAL]} functional components, "
        
        interpretation = interpretation.rstrip(", ")
        
        # Overall assessment
        if most_severe_impact in [NON_FUNCTIONAL, FRAMESHIFT]:
            interpretation += ". Overall: NON-FUNCTIONAL"
        elif most_severe_impact in [TRUNCATED]:
            interpretation += ". Overall: SEVERELY DAMAGED"
        elif most_severe_impact == CHANGE_OF_FUNCTION:
            interpretation += ". Overall: ALTERED FUNCTION"
        elif most_severe_impact == PARTIAL_FUNCTION:
            interpretation += ". Overall: PARTIAL FUNCTION"
        else:
            interpretation += ". Overall: FUNCTIONAL"
        
        details = [f"Aggregated from {len(child_impacts)} sub-features"]
        for impact_type, count in sorted(impact_counts.items(), key=lambda x: -IMPACT_SEVERITY[x[0]]):
            if count > 0:
                details.append(f"  {count} {impact_type}")
        
        return {
            "impact": most_severe_impact,
            "details": details,
            "identity": avg_identity,
            "coverage": avg_coverage,
            "mismatches": sum(c.get('mismatches', 0) for c in child_impacts),
            "insertions": sum(c.get('insertions', 0) for c in child_impacts),
            "deletions": sum(c.get('deletions', 0) for c in child_impacts),
            "duplications": total_duplications,
            "interpretation": interpretation,
            "child_count": len(child_impacts),
            "impact_counts": dict(impact_counts),
            "feature_id": parent_id  # Add feature_id for parent tracking
        }

    def _get_alignment_visualization(self, feature_id: str, aln: Dict[str, Any]) -> str:
        """Generates a colorized string showing the base-level alignment."""
        if not aln or not aln.get('cigar'):
            return "  (No CIGAR string available for visualization)"
        
        feature = self.gff.get_feature(feature_id)
        
        # Get sequences
        try:
            # Get reference sequence
            ref_seq = self.reference_genome.fetch(
                feature.seqid, 
                feature.start - 1, 
                feature.end
            )
            if feature.strand == '-':
                ref_seq = self._reverse_complement(ref_seq)
        except Exception as e:
            return f"  (Could not fetch reference sequence: {e})"
        
        try:
            # Get target sequence region
            target_seq = self.path_seqs.fetch(
                aln['target_name'],
                aln['target_start'],
                aln['target_end']
            )
        except Exception as e:
            return f"  (Could not fetch target sequence: {e})"
        
        # Parse CIGAR and build visualization
        cigar_tuples = re.findall(r'(\d+)([MIDNSHPX=])', aln['cigar'])
        
        query_pos = 0
        target_pos = 0
        
        viz_lines = []
        query_line = ""
        match_line = ""
        target_line = ""
        
        for length_str, op in cigar_tuples:
            length = int(length_str)
            
            if op in ['M', '=', 'X']:  # Match/mismatch
                for i in range(length):
                    if query_pos < len(ref_seq) and target_pos < len(target_seq):
                        q_base = ref_seq[query_pos]
                        t_base = target_seq[target_pos]
                        
                        if q_base.upper() == t_base.upper():
                            query_line += q_base
                            match_line += "|"
                            target_line += t_base
                        else:
                            query_line += self.colors['RED'] + q_base + self.colors['ENDC']
                            match_line += self.colors['YELLOW'] + "X" + self.colors['ENDC']
                            target_line += self.colors['RED'] + t_base + self.colors['ENDC']
                    
                    query_pos += 1
                    target_pos += 1
                    
            elif op == 'I':  # Insertion in query
                for i in range(length):
                    if query_pos < len(ref_seq):
                        query_line += self.colors['GREEN'] + ref_seq[query_pos] + self.colors['ENDC']
                        match_line += self.colors['GREEN'] + "+" + self.colors['ENDC']
                        target_line += self.colors['GREEN'] + "-" + self.colors['ENDC']
                        query_pos += 1
                        
            elif op == 'D':  # Deletion in query
                for i in range(length):
                    if target_pos < len(target_seq):
                        query_line += self.colors['BLUE'] + "-" + self.colors['ENDC']
                        match_line += self.colors['BLUE'] + "-" + self.colors['ENDC']
                        target_line += self.colors['BLUE'] + target_seq[target_pos] + self.colors['ENDC']
                        target_pos += 1
            
            elif op in ['S', 'H']:  # Soft/hard clipping
                if op == 'S':
                    query_pos += length
        
        # Format output in chunks
        chunk_size = 80
        result_lines = []
        
        for i in range(0, len(query_line), chunk_size):
            chunk_end = min(i + chunk_size, len(query_line))
            result_lines.append(f"  Ref:    {query_line[i:chunk_end]}")
            result_lines.append(f"          {match_line[i:chunk_end]}")
            result_lines.append(f"  Target: {target_line[i:chunk_end]}")
            result_lines.append("")
        
        # Add legend
        result_lines.append(f"  Legend: {self.colors['RED']}Mismatch{self.colors['ENDC']} "
                          f"{self.colors['GREEN']}Insertion{self.colors['ENDC']} "
                          f"{self.colors['BLUE']}Deletion{self.colors['ENDC']}")
        
        return "\n".join(result_lines)

    @staticmethod
    def _reverse_complement(dna: str) -> str:
        """Returns the reverse complement of a DNA sequence."""
        complement = str.maketrans('ATCGNRY', 'TAGCNYR')
        return dna.upper().translate(complement)[::-1]

    def print_summary(self):
        """Prints a hierarchical summary of the impact analysis."""
        print("\n" + "="*80)
        print(f"{self.colors['BOLD']}FUNCTIONAL IMPACT ANALYSIS SUMMARY{self.colors['ENDC']}")
        print("="*80)

        sorted_paths = sorted(self.results.keys())
        if not sorted_paths:
            print("No alignments found to analyze.")
            return

        # Summary statistics
        total_features_analyzed = sum(len(path_results) for path_results in self.results.values())
        print(f"\nAnalyzed {total_features_analyzed} feature alignments across {len(sorted_paths)} paths")
        
        for path_name in sorted_paths:
            print(f"\n{self.colors['CYAN']}{'='*80}{self.colors['ENDC']}")
            print(f"{self.colors['BOLD']}{self.colors['CYAN']}Path: {path_name}{self.colors['ENDC']}")
            print(f"{self.colors['CYAN']}{'='*80}{self.colors['ENDC']}")
            
            # Get root features to display
            root_features = self.gff.roots
            if self.target_genes:
                filtered_roots = []
                for feature in self.gff.roots:
                    gene_name = feature.attributes.get('gene_name', feature.attributes.get('Name', ['']))[0]
                    if gene_name in self.target_genes or feature.id in self.target_genes:
                        filtered_roots.append(feature)
                root_features = filtered_roots

            sorted_roots = sorted(root_features, key=lambda f: f.start)

            for root_feature in sorted_roots:
                if root_feature.id in self.results[path_name]:
                    self._print_feature_summary(root_feature.id, path_name, indent=0)

    def _print_feature_summary(self, feature_id: str, path_name: str, indent: int):
        """Recursively prints the summary for a feature and its children."""
        feature = self.gff.get_feature(feature_id)
        result = self.results[path_name].get(feature_id)
        if not feature or not result:
            return

        prefix = "  " * indent
        
        # Get display name
        display_name = feature.id
        if feature.featuretype == 'gene':
            # For genes, prefer a common name
            display_name = feature.attributes.get('gene_name', 
                          feature.attributes.get('Name', [feature.id]))[0]
        
        # Choose color based on impact
        impact_color = {
            NON_FUNCTIONAL: self.colors['RED'],
            FRAMESHIFT: self.colors['RED'],
            TRUNCATED: self.colors['RED'],
            CHANGE_OF_FUNCTION: self.colors['YELLOW'],
            PARTIAL_FUNCTION: self.colors['YELLOW'],
            DUPLICATED: self.colors['MAGENTA'],
            FUNCTIONAL: self.colors['GREEN'],
            NO_CHANGE: self.colors['GREEN'],
            NOT_FOUND: self.colors['BLUE'],
        }.get(result['impact'], '')

        # Print main feature line
        print(f"{prefix}{self.colors['BOLD']}├─ {feature.featuretype}{self.colors['ENDC']} "
              f"'{display_name}' [{feature.seqid}:{feature.start:,}-{feature.end:,}]")
        print(f"{prefix}│  {self.colors['BOLD']}Status:{self.colors['ENDC']} "
              f"{impact_color}{result['impact']}{self.colors['ENDC']}")
        
        # Print statistics
        identity = result.get('identity', 0.0)
        coverage = result.get('coverage', 0.0)
        
        stats_line = f"{prefix}│  {self.colors['BOLD']}Stats:{self.colors['ENDC']} "
        if coverage > 0:
            stats_line += f"Coverage: {coverage:.1f}%, Identity: {identity:.1f}%"
        else:
            stats_line += "Not aligned"
        
        if result.get('mismatches', 0) > 0 or result.get('insertions', 0) > 0 or result.get('deletions', 0) > 0:
            stats_line += f" | Changes: {result.get('mismatches', 0)} subst, "
            stats_line += f"{result.get('insertions', 0)} ins, {result.get('deletions', 0)} del"
        
        print(stats_line)
        
        # Print interpretation
        if result.get('interpretation'):
            print(f"{prefix}│  {self.colors['BOLD']}Interpretation:{self.colors['ENDC']} "
                  f"{result['interpretation']}")
        
        # Print details if any
        if result.get('details') and len(result['details']) > 1:
            print(f"{prefix}│  {self.colors['BOLD']}Details:{self.colors['ENDC']}")
            for detail in result['details']:
                print(f"{prefix}│    • {detail}")
        
        # Print duplications if present
        if result.get('duplications', 0) > 1:
            print(f"{prefix}│  {self.colors['MAGENTA']}⚠ Feature duplicated "
                  f"{result['duplications']} times{self.colors['ENDC']}")
        
        # Show alignment visualization if requested
        if self.show_alignments and result.get('alignment') and feature.featuretype in ['CDS', 'exon']:
            print(f"{prefix}│  {self.colors['BOLD']}Alignment:{self.colors['ENDC']}")
            viz = self._get_alignment_visualization(feature_id, result['alignment'])
            for line in viz.split('\n'):
                print(f"{prefix}│  {line}")
        
        print(f"{prefix}│")
        
        # Process children
        child_features = self.gff.get_children(feature_id)
        if child_features:
            # Children are already sorted by start position from get_children
            for i, child_feature in enumerate(child_features):
                if child_feature.id in self.results[path_name]:
                    self._print_feature_summary(child_feature.id, path_name, indent + 1)


def main():
    """Main function to run the impact analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze functional impact of variants from PAF alignments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Impact Categories:
  - Functional: Feature is intact and likely functional
  - Non-functional: Feature is severely damaged or missing
  - Change in Function: Feature has significant changes that may alter function
  - Partial Function: Feature is partially present/damaged
  - Frameshift: CDS has frameshift mutation
  - Truncated: Feature is significantly shortened
  - Duplicated: Feature appears multiple times
  - Not Found: Feature is not aligned to this path
        """
    )
    parser.add_argument("--paf", required=True, type=Path, help="Input PAF file from hierarchical aligner")
    parser.add_argument("--gff", required=True, type=Path, help="GFF3 annotation file")
    parser.add_argument("--reference-fasta", required=True, type=Path, help="Reference genome FASTA")
    parser.add_argument("--path-fasta", required=True, type=Path, help="Multi-sample FASTA with paths")
    parser.add_argument("--genes", type=str, help="Comma-separated list of gene names to analyze")
    parser.add_argument("--show-alignments", action="store_true", help="Show detailed alignment visualizations")
    parser.add_argument("--no-color", action="store_true", help="Disable colored output")
    parser.add_argument("--output-json", type=Path, help="Save results to JSON file")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    setup_logging(args.verbose)

    try:
        # Parse target genes
        target_genes = []
        if args.genes:
            target_genes = [g.strip() for g in args.genes.split(',')]
            logging.info(f"Targeting genes: {', '.join(target_genes)}")
        
        # Load data
        logging.info("Loading input files...")
        gff_parser = GFFParser(args.gff)
        paf_parser = PAFParser(args.paf)

        # Run analysis
        analyzer = ImpactAnalyzer(
            gff_parser,
            paf_parser,
            args.reference_fasta,
            args.path_fasta,
            target_genes,
            use_color=not args.no_color,
            show_alignments=args.show_alignments
        )
        
        results = analyzer.analyze()
        analyzer.print_summary()
        
        # Save JSON output if requested
        if args.output_json:
            logging.info(f"Saving results to {args.output_json}")
            with open(args.output_json, 'w') as f:
                # Convert results to JSON-serializable format
                json_results = {}
                for path_name, path_results in results.items():
                    json_results[path_name] = {}
                    for feature_id, feature_result in path_results.items():
                        # Remove non-serializable alignment objects
                        clean_result = {k: v for k, v in feature_result.items() 
                                      if k not in ['alignment', 'all_alignments']}
                        json_results[path_name][feature_id] = clean_result
                json.dump(json_results, f, indent=2)
            logging.info(f"Results saved to {args.output_json}")

    except FileNotFoundError as e:
        logging.error(f"Error: Input file not found - {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
