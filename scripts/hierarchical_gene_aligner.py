#!/usr/bin/env python3

import argparse
import json
import logging
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple, Set

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
    def __init__(self, feature_id: str, feature_type: str, mapq: int, target_start: int, target_end: int, identity: float, cigar: str = None):
        self.feature_id = feature_id
        self.feature_type = feature_type
        self.mapq = mapq
        self.target_start = target_start
        self.target_end = target_end
        self.identity = identity
        self.cigar = cigar
        self.children: List[AlignmentResult] = []

    def __repr__(self) -> str:
        return f"AlignmentResult(id={self.feature_id}, type={self.feature_type}, mapq={self.mapq})"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'feature_id': self.feature_id,
            'feature_type': self.feature_type,
            'mapq': self.mapq,
            'target_start': self.target_start,
            'target_end': self.target_end,
            'identity': self.identity,
            'cigar': self.cigar,
            'children': [child.to_dict() for child in self.children]
        }

# --- Core Classes ---

class GFFProcessor:
    """Handles GFF parsing and feature hierarchy for a specific gene."""
    def __init__(self, gff_file: Path, gene_identifier: str):
        self.gff_file = gff_file
        self.gene_identifier = gene_identifier
        self.gene_feature = None
        self.features_by_id = {}
        self.children_map = defaultdict(list)
        
        # Load only the relevant gene and its descendants
        self._load_gene_features_optimized()
    
    def _parse_attributes(self, attr_string: str) -> Dict[str, List[str]]:
        """Parse GFF attribute string into a dictionary."""
        attributes = defaultdict(list)
        for item in attr_string.strip().split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attributes[key].append(value)
        return dict(attributes)
    
    def _load_gene_features_optimized(self):
        """Load only the gene and its related features from the GFF file in a single pass."""
        logging.info(f"Searching for gene '{self.gene_identifier}' and loading its features...")
        
        # Single pass to find gene and collect all potential descendants
        gene_found = False
        gene_chr = None
        gene_start = None
        gene_end = None
        all_features = []  # Store all features temporarily
        parent_to_children = defaultdict(set)
        
        # Count lines for progress bar
        total_lines = sum(1 for line in open(self.gff_file, 'r') if not line.startswith('#'))
        
        with open(self.gff_file, 'r') as f:
            with tqdm(total=total_lines, desc="Processing GFF file", unit=" lines") as pbar:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    pbar.update(1)
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue
                    
                    feature_type = parts[2]
                    attributes = self._parse_attributes(parts[8])
                    feature_id = attributes.get('ID', [''])[0]
                    
                    # Check if this is our gene
                    if not gene_found and feature_type == 'gene':
                        # Check multiple possible attribute names for gene name
                        feature_name = attributes.get('Name', [''])[0]  # Standard GFF3
                        if not feature_name:
                            feature_name = attributes.get('gene_name', [''])[0]  # GENCODE/Ensembl style
                        
                        # Also check gene_id for Ensembl IDs
                        gene_id = attributes.get('gene_id', [''])[0]
                        
                        if (feature_id == self.gene_identifier or 
                            feature_name == self.gene_identifier or
                            gene_id == self.gene_identifier):
                            # Create the gene feature
                            self.gene_feature = gffutils.Feature(
                                seqid=parts[0],
                                source=parts[1],
                                featuretype=parts[2],
                                start=int(parts[3]),
                                end=int(parts[4]),
                                score=parts[5],
                                strand=parts[6],
                                frame=parts[7],
                                attributes=attributes,
                                id=feature_id
                            )
                            self.features_by_id[feature_id] = self.gene_feature
                            gene_found = True
                            gene_chr = parts[0]
                            gene_start = int(parts[3])
                            gene_end = int(parts[4])
                            logging.info(f"Found gene '{self.gene_identifier}' (ID: {feature_id}, Name: {feature_name}) at {gene_chr}:{gene_start}-{gene_end}")
                    
                    # If gene is found, collect features in the same region (optimization)
                    if gene_found:
                        # Only process features on the same chromosome and overlapping the gene region
                        if parts[0] == gene_chr:
                            feat_start = int(parts[3])
                            feat_end = int(parts[4])
                            
                            # Check if feature overlaps with gene region (with some buffer)
                            if feat_start <= gene_end + 10000 and feat_end >= gene_start - 10000:
                                all_features.append((parts, attributes, feature_id))
                                
                                # Build parent-child relationships
                                parent_ids = attributes.get('Parent', [])
                                for parent_id in parent_ids:
                                    parent_to_children[parent_id].add(feature_id)
        
        if not gene_found:
            raise ValueError(f"Gene '{self.gene_identifier}' not found in GFF file")
        
        # Find all descendants of the gene
        relevant_ids = {self.gene_feature.id}
        to_process = [self.gene_feature.id]
        
        while to_process:
            current_id = to_process.pop()
            for child_id in parent_to_children.get(current_id, []):
                if child_id not in relevant_ids:
                    relevant_ids.add(child_id)
                    to_process.append(child_id)
        
        logging.info(f"Found {len(relevant_ids)} features related to gene '{self.gene_identifier}'")
        
        # Now load only the relevant features
        feature_count = 0
        for parts, attributes, feature_id in all_features:
            if feature_id in relevant_ids and feature_id != self.gene_feature.id:
                # Create the feature
                feature = gffutils.Feature(
                    seqid=parts[0],
                    source=parts[1],
                    featuretype=parts[2],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    score=parts[5],
                    strand=parts[6],
                    frame=parts[7],
                    attributes=attributes,
                    id=feature_id
                )
                self.features_by_id[feature_id] = feature
                feature_count += 1
                
                # Build parent-child relationships
                parent_ids = attributes.get('Parent', [])
                for parent_id in parent_ids:
                    if parent_id in relevant_ids:
                        self.children_map[parent_id].append(feature_id)
        
        logging.info(f"Loaded {feature_count + 1} features for gene '{self.gene_identifier}'")
    
    def find_gene(self, gene_identifier: str) -> Optional[gffutils.Feature]:
        """Returns the gene feature if it matches the identifier."""
        if self.gene_feature:
            return self.gene_feature
        return None
    
    def get_children(self, feature: gffutils.Feature) -> Iterator[gffutils.Feature]:
        """Yields all direct children of a feature, sorted by position."""
        child_ids = self.children_map.get(feature.id, [])
        children = []
        for child_id in child_ids:
            if child_id in self.features_by_id:
                children.append(self.features_by_id[child_id])
        
        # Sort by start position
        children.sort(key=lambda f: f.start)
        for child in children:
            yield child


class SequenceExtractor:
    """Extracts feature sequences from a reference FASTA."""
    def __init__(self, reference_fasta: Path):
        logging.info(f"Loading reference FASTA: {reference_fasta}")
        self.fasta = pysam.FastaFile(str(reference_fasta))
        self.available_chromosomes = set(self.fasta.references)
        logging.info(f"Reference FASTA loaded with {len(self.fasta.references)} sequences")

    def get_sequence(self, feature: gffutils.Feature) -> Optional[str]:
        """Extracts a feature's sequence, handling strand."""
        # Check if chromosome exists in reference
        if feature.seqid not in self.available_chromosomes:
            logging.warning(f"Chromosome {feature.seqid} not found in reference FASTA")
            return None
        
        try:
            # GFF coordinates are 1-based, pysam uses 0-based
            seq = self.fasta.fetch(feature.seqid, feature.start - 1, feature.end)
        except (KeyError, ValueError) as e:
            logging.warning(f"Could not fetch sequence for {feature.id} ({feature.seqid}:{feature.start}-{feature.end}): {e}")
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
        self.raw_paf_lines = []  # Store raw PAF output

    def _get_all_features_recursive(self, feature: gffutils.Feature, feature_list: List[gffutils.Feature]):
        """Recursively gather a feature and all its descendants."""
        feature_list.append(feature)
        for child in self.gff.get_children(feature):
            self._get_all_features_recursive(child, feature_list)

    def _try_exact_match(self, feature: gffutils.Feature, feature_seq: str, target_fasta: Path, feature_id: str) -> Optional[List[Dict[str, Any]]]:
        """Try to find exact matches for a feature sequence in the target FASTA."""
        logging.debug(f"Trying exact match for feature {feature_id} ({len(feature_seq)}bp)")
        
        # Load target sequences and search for exact matches
        target_seqs = pysam.FastaFile(str(target_fasta))
            results = []
            
            for ref_name in target_seqs.references:
                # For reference path, we expect to find the feature at its original position
                is_reference = "grch38" in ref_name.lower() or "hg38" in ref_name.lower()
                
                if is_reference:
                    # For reference, check the exact expected position first
                    try:
                        # Get the sequence at the expected position
                        target_len = target_seqs.get_reference_length(ref_name)
                        expected_seq = target_seqs.fetch(ref_name, feature.start - 1, feature.end)
                        if feature.strand == '-':
                            expected_seq = self.seq._reverse_complement(expected_seq)
                        
                        if expected_seq.upper() == feature_seq.upper():
                            results.append({
                                "query_name": feature_id,
                                "query_len": len(feature_seq),
                                "target_name": ref_name,
                                "target_len": target_len,
                                "target_start": feature.start - 1,  # Convert to 0-based
                                "target_end": feature.end,
                                "matches": len(feature_seq),
                                "align_len": len(feature_seq),
                                "mapq": 60,  # High quality for exact match
                                "cigar": f"{len(feature_seq)}M",
                                "is_exact": True
                            })
                            logging.debug(f"Found exact match for {feature_id} in {ref_name} at expected position {feature.start - 1}")
                            # For reference, only take the expected position, so continue to next haplotype
                            continue
                    except Exception as e:
                        logging.debug(f"Could not check expected position for {feature_id}: {e}")
                    # If not found at exact position on reference, let minimap2 handle it.
                else:
                    # For non-reference paths, search the entire sequence
                    ref_seq = target_seqs.fetch(ref_name)
                    target_len = len(ref_seq)
                    pos = ref_seq.upper().find(feature_seq.upper())
                    if pos != -1:
                        results.append({
                            "query_name": feature_id,
                            "query_len": len(feature_seq),
                            "target_name": ref_name,
                            "target_len": target_len,
                            "target_start": pos,
                            "target_end": pos + len(feature_seq),
                            "matches": len(feature_seq),
                            "align_len": len(feature_seq),
                            "mapq": 60,  # High quality for exact match
                            "cigar": f"{len(feature_seq)}M",
                            "is_exact": True
                        })
                        logging.debug(f"Found exact match for {feature_id} in {ref_name} at position {pos}")
            
            target_seqs.close()
            return results if results else None

    def _choose_minimap2_preset(self, feature_length: int) -> Tuple[str, List[str]]:
        """Choose appropriate minimap2 preset based on feature length."""
        if feature_length < 20:
            # For very small features like start/stop codons (3bp), use extremely sensitive settings
            return "sr", ["-k", "3", "-w", "1", "--min-dp-score", "1", "-m", "1", "-n", "1", "-A", "2", "-B", "2", "-O", "4,2", "-E", "2,1", "--score-N", "0"]
        elif feature_length < 50:
            # For small features, use very sensitive settings
            return "sr", ["-k", "7", "-w", "1", "--score-N", "0", "-A", "2", "-B", "4"]
        elif feature_length < 200:
            # For small exons/UTRs
            return "sr", ["-k", "11", "-w", "3"]
        elif feature_length < 1000:
            # For medium features
            return "splice", []
        else:
            # For large features
            return "asm20", []

    def _run_minimap2_batch(self, target_fasta: Path, features_by_size: Dict[str, List[Tuple[gffutils.Feature, str]]]) -> Dict[str, List[Dict[str, Any]]]:
        """Run minimap2 with different settings for different feature sizes."""
        all_alignments = defaultdict(list)
        all_paf_lines = []
        
        # 1. Handle features with exact matching
        if "exact" in features_by_size and features_by_size["exact"]:
            logging.info(f"Running exact match for {len(features_by_size['exact'])} features (<=50bp)...")
            for feature, seq in tqdm(features_by_size["exact"], desc="Exact matching", unit=" features"):
                exact_matches = self._try_exact_match(feature, seq, target_fasta, feature.id)
                if exact_matches:
                    for match in exact_matches:
                        all_alignments[match["target_name"]].append(match)
                        # Create PAF line for exact match
                        paf_line = f"{match['query_name']}\t{match['query_len']}\t0\t{match['query_len']}\t+\t"
                        paf_line += f"{match['target_name']}\t{match['target_len']}\t{match['target_start']}\t{match['target_end']}\t"
                        paf_line += f"{match['matches']}\t{match['align_len']}\t{match['mapq']}\tcg:Z:{match['cigar']}"
                        all_paf_lines.append(paf_line)

        # 2. Run minimap2 on the other categories
        for size_category, features in features_by_size.items():
            if size_category == "exact" or not features:
                continue
                
            # Create temp file for this batch
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa") as query_fasta:
                query_path = Path(query_fasta.name)
                
                # Write sequences
                feature_lengths = {}
                for feature, seq in features:
                    query_fasta.write(f">{feature.id}\n{seq}\n")
                    feature_lengths[feature.id] = len(seq)
                query_fasta.flush()
                
                # Get appropriate preset
                avg_length = sum(len(seq) for _, seq in features) / len(features) if features else 0
                preset, extra_args = self._choose_minimap2_preset(int(avg_length))
                
                logging.info(f"Running minimap2 for {len(features)} {size_category} features (preset: {preset}, avg length: {avg_length:.0f}bp)...")
                
                cmd = [
                    "minimap2",
                    "-x", preset,
                    "-c",  # Output CIGAR
                    "-N", "1", # Report at most 1 secondary alignment to reduce noise
                    "-t", str(self.threads),
                ] + extra_args + [
                    str(target_fasta),
                    str(query_path)
                ]
                
                # For tiny and small features, don't suppress no-hit output
                if size_category not in ["tiny", "small"]:
                    cmd.insert(3, "--paf-no-hit")
                
                try:
                    process = subprocess.run(cmd, capture_output=True, text=True, check=True)
                    
                    lines = process.stdout.strip().split('\n')
                    alignment_count = 0
                    for line in lines:
                        if not line:
                            continue
                        
                        all_paf_lines.append(line)
                        parts = line.split('\t')
                        if len(parts) < 12:
                            continue
                        
                        # Extract CIGAR string
                        cigar = None
                        for tag in parts[12:]:
                            if tag.startswith("cg:Z:"):
                                cigar = tag[5:]
                                break
                        
                        paf_record = {
                            "query_name": parts[0],
                            "query_len": int(parts[1]),
                            "target_name": parts[5],
                            "target_start": int(parts[7]),
                            "target_end": int(parts[8]),
                            "matches": int(parts[9]),
                            "align_len": int(parts[10]),
                            "mapq": int(parts[11]),
                            "cigar": cigar
                        }
                        all_alignments[paf_record["target_name"]].append(paf_record)
                        alignment_count += 1
                    
                    logging.info(f"Found {alignment_count} alignments for {size_category} features")
                        
                except subprocess.CalledProcessError as e:
                    logging.warning(f"minimap2 failed for {size_category} features: {e.stderr}")
                except FileNotFoundError:
                    logging.error("minimap2 not found. Please ensure minimap2 is installed and in your PATH.")
                    sys.exit(1)
                finally:
                    try:
                        query_path.unlink()
                    except Exception:
                        pass
        
        self.raw_paf_lines = all_paf_lines
        return all_alignments

    def align_gene_to_haplotypes(self, gene: gffutils.Feature, haplotypes_fasta: Path) -> Dict[str, AlignmentResult]:
        """Aligns a gene and its sub-features to all sequences in a multi-sample FASTA."""
        # 1. Gather all features and their sequences
        all_features = []
        logging.info("Gathering all features in gene hierarchy...")
        self._get_all_features_recursive(gene, all_features)
        
        logging.info(f"Found {len(all_features)} features (gene and descendants) to align.")
        
        # Categorize features by size for appropriate alignment strategies
        features_by_size = {
            "exact": [],     # <= 50bp
            "small": [],     # 51-200bp
            "medium": [],    # 200-1000bp
            "large": []      # > 1000bp
        }
        
        skipped_count = 0
        for feature in tqdm(all_features, desc="Extracting feature sequences", unit=" features"):
            seq = self.seq.get_sequence(feature)
            if seq and len(seq) > 0:
                length = len(seq)
                if length <= 50:
                    features_by_size["exact"].append((feature, seq))
                elif length < 200:
                    features_by_size["small"].append((feature, seq))
                elif length < 1000:
                    features_by_size["medium"].append((feature, seq))
                else:
                    features_by_size["large"].append((feature, seq))
            else:
                skipped_count += 1
                logging.debug(f"Skipped feature {feature.id} (no sequence)")
        
        if skipped_count > 0:
            logging.warning(f"Skipped {skipped_count} features with no extractable sequence")
        
        total_features = sum(len(f) for f in features_by_size.values())
        if total_features == 0:
            logging.error("No feature sequences could be extracted.")
            return {}
        
        logging.info(f"Feature size distribution: exact={len(features_by_size['exact'])}, "
                    f"small={len(features_by_size['small'])}, medium={len(features_by_size['medium'])}, "
                    f"large={len(features_by_size['large'])}")
        
        # 2. Run minimap2 with appropriate settings for each size category
        all_alignments = self._run_minimap2_batch(haplotypes_fasta, features_by_size)
        
        # 3. Build hierarchical results from flat alignment list
        results = {}
        hap_names = list(all_alignments.keys())
        
        if not hap_names:
            logging.warning("No alignments found for any haplotype.")
            return results
        
        logging.info(f"Building hierarchical alignment results for {len(hap_names)} haplotypes...")
        
        for hap_name in tqdm(hap_names, desc="Building hierarchical results", unit=" haplotypes"):
            is_reference = "grch38" in hap_name.lower() or "hg38" in hap_name.lower()
            
            # Group alignments by query name
            grouped_by_query = defaultdict(list)
            for aln in all_alignments[hap_name]:
                grouped_by_query[aln['query_name']].append(aln)
            
            alignments_for_hap = {}
            for query_name, aln_list in grouped_by_query.items():
                if aln_list:
                    if is_reference:
                        # For reference, be very strict. Prioritize an exact match at the original location.
                        feature = self.gff.features_by_id.get(query_name)
                        if not feature: continue

                        # 1. Look for our own exact match at the precise original location.
                        exact_alns = [a for a in aln_list if a.get("is_exact") and a['target_start'] == feature.start - 1]
                        if exact_alns:
                            alignments_for_hap[query_name] = exact_alns[0]
                            continue

                        # 2. If not found, look for a high-quality minimap2 alignment at the original location.
                        candidate_alns = [
                            aln for aln in aln_list 
                            if aln['target_start'] < feature.end and aln['target_end'] > feature.start - 1
                        ]
                        if candidate_alns:
                            best_aln = max(candidate_alns, key=lambda x: (x['mapq'], x['matches']/x['align_len'] if x['align_len'] > 0 else 0))
                            alignments_for_hap[query_name] = best_aln
                        else:
                            # 3. Fallback: if no alignment is at the original location (e.g., large SV), take the best one anywhere.
                            best_aln = max(aln_list, key=lambda x: (x['mapq'], x['matches']/x['align_len'] if x['align_len'] > 0 else 0))
                            alignments_for_hap[query_name] = best_aln
                    else:
                        # For non-reference, we are more lenient. Just take the best alignment found.
                        # Multiple alignments for the same feature are handled as duplications downstream.
                        best_aln = max(aln_list, key=lambda x: (x['mapq'], x['matches']/x['align_len'] if x['align_len'] > 0 else 0))
                        alignments_for_hap[query_name] = best_aln
            
            # Find the top-level gene alignment
            gene_aln_data = alignments_for_hap.get(gene.id)
            if not gene_aln_data:
                logging.debug(f"No alignment found for gene {gene.id} on haplotype {hap_name}")
                continue
            
            identity = gene_aln_data['matches'] / gene_aln_data['align_len'] if gene_aln_data['align_len'] > 0 else 0
            gene_result = AlignmentResult(
                feature_id=gene.id,
                feature_type=gene.featuretype,
                mapq=gene_aln_data['mapq'],
                target_start=gene_aln_data['target_start'],
                target_end=gene_aln_data['target_end'],
                identity=identity,
                cigar=gene_aln_data.get('cigar')
            )
            
            # Recursively build the tree for children
            # Pass perfect identity info down if parent has 100% identity
            self._build_result_tree(gene, gene_result, alignments_for_hap, 
                                  inherit_perfect=(identity >= 0.999 and is_reference))
            results[hap_name] = gene_result
        
        logging.info(f"Successfully built alignment results for {len(results)} haplotypes")
        return results

    def _build_result_tree(self, parent_feature: gffutils.Feature, parent_result: AlignmentResult, 
                          alignments: Dict[str, Dict], inherit_perfect: bool = False):
        """Recursively constructs the alignment result tree, applying hierarchical constraints."""
        for child_feature in self.gff.get_children(parent_feature):
            # If parent has perfect identity on reference, child inherits it
            if inherit_perfect:
                # Calculate relative position within parent
                parent_length = parent_feature.end - parent_feature.start
                child_rel_start = child_feature.start - parent_feature.start
                child_rel_end = child_feature.end - parent_feature.start
                
                child_result = AlignmentResult(
                    feature_id=child_feature.id,
                    feature_type=child_feature.featuretype,
                    mapq=60,  # High quality inherited from parent
                    target_start=parent_result.target_start + child_rel_start,
                    target_end=parent_result.target_start + child_rel_end,
                    identity=1.0,
                    cigar=None  # Inherited, not directly aligned
                )
                parent_result.children.append(child_result)
                # Recurse with perfect inheritance
                self._build_result_tree(child_feature, child_result, alignments, inherit_perfect=True)
                continue
            
            child_aln_data = alignments.get(child_feature.id)
            if not child_aln_data:
                # If no alignment found but parent has good alignment, mark as likely present
                if parent_result.identity >= 0.95:
                    logging.debug(f"No alignment for {child_feature.id}, but parent has {parent_result.identity:.1%} identity")
                continue
            
            # Hierarchical constraint: child must be within parent's aligned region (with tolerance)
            tolerance = 100  # Increase tolerance for small features
            if not (child_aln_data['target_start'] >= parent_result.target_start - tolerance and 
                    child_aln_data['target_end'] <= parent_result.target_end + tolerance):
                logging.debug(f"Skipping child {child_feature.id} as its alignment [{child_aln_data['target_start']}-{child_aln_data['target_end']}] "
                            f"is outside parent {parent_feature.id}'s region [{parent_result.target_start}-{parent_result.target_end}]")
                continue
            
            identity = child_aln_data['matches'] / child_aln_data['align_len'] if child_aln_data['align_len'] > 0 else 0
            child_result = AlignmentResult(
                feature_id=child_feature.id,
                feature_type=child_feature.featuretype,
                mapq=child_aln_data['mapq'],
                target_start=child_aln_data['target_start'],
                target_end=child_aln_data['target_end'],
                identity=identity,
                cigar=child_aln_data.get('cigar')
            )
            parent_result.children.append(child_result)
            
            # Recurse, potentially with perfect inheritance if this child has perfect identity
            self._build_result_tree(child_feature, child_result, alignments, 
                                  inherit_perfect=(identity >= 0.999))
    
    def get_raw_paf_lines(self) -> List[str]:
        """Returns the raw PAF output lines from the last alignment."""
        return self.raw_paf_lines

# --- Output ---

def print_summary(results: Dict[str, AlignmentResult]):
    """Prints a formatted summary of the alignment results."""
    print("\n" + "="*80)
    print("HIERARCHICAL ALIGNMENT SUMMARY")
    print("="*80)
    
    if not results:
        print("\nNo alignment results found.")
        return

    print(f"\nTotal haplotypes with alignments: {len(results)}")
    
    for hap_name, gene_result in results.items():
        print(f"\n--- Haplotype: {hap_name} ---")
        _print_result_recursively(gene_result, indent=0)

def _print_result_recursively(result: AlignmentResult, indent: int):
    """Helper to print results hierarchically."""
    prefix = "  " * indent
    print(f"{prefix}- {result.feature_type} '{result.feature_id}':")
    print(f"{prefix}  - Aligned to region: {result.target_start:,}-{result.target_end:,}")
    print(f"{prefix}  - MAPQ: {result.mapq}, Identity: {result.identity:.2%}")
    if result.cigar is None and result.identity >= 0.999:
        print(f"{prefix}  - (Inherited from parent's perfect alignment)")

    # Sort children by their start position for logical output
    sorted_children = sorted(result.children, key=lambda r: r.target_start)
    for child in sorted_children:
        _print_result_recursively(child, indent + 1)

def save_results(results: Dict[str, AlignmentResult], output_path: Path):
    """Save alignment results to JSON file."""
    results_dict = {
        hap_name: result.to_dict() 
        for hap_name, result in results.items()
    }
    
    with open(output_path, 'w') as f:
        json.dump(results_dict, f, indent=2)
    
    logging.info(f"Results saved to {output_path}")

def save_paf_alignments(paf_lines: List[str], output_path: Path):
    """Save raw PAF alignment data to file."""
    with open(output_path, 'w') as f:
        for line in paf_lines:
            f.write(line + '\n')
    
    logging.info(f"PAF alignments saved to {output_path}")

def save_alignments_with_metadata(paf_lines: List[str], gff_processor: GFFProcessor, output_path: Path):
    """Save PAF alignments with additional metadata about features."""
    with open(output_path, 'w') as f:
        # Write header
        f.write("# PAF alignments with feature metadata\n")
        f.write("# Additional columns: feature_type, parent_id, gene_name\n")
        
        for line in paf_lines:
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 12:
                f.write(line + '\n')
                continue
            
            query_name = parts[0]
            
            # Get feature metadata
            feature = gff_processor.features_by_id.get(query_name)
            if feature:
                feature_type = feature.featuretype
                parent_ids = feature.attributes.get('Parent', [])
                parent_id = parent_ids[0] if parent_ids else 'NA'
                gene_name = feature.attributes.get('gene_name', ['NA'])[0]
                
                # Add metadata as additional columns
                extended_line = f"{line}\tft:Z:{feature_type}\tpi:Z:{parent_id}\tgn:Z:{gene_name}"
                f.write(extended_line + '\n')
            else:
                f.write(line + '\n')
    
    logging.info(f"Extended PAF alignments saved to {output_path}")

# --- Main Execution ---

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Hierarchically align a gene's features against multiple sequences.")
    parser.add_argument("--multi-fasta", required=True, type=Path, help="Multi-sample FASTA file of haplotypes.")
    parser.add_argument("--gff", required=True, type=Path, help="GFF3 annotation file.")
    parser.add_argument("--reference", required=True, type=Path, help="Reference genome FASTA file.")
    parser.add_argument("--gene", required=True, type=str, help="Name or ID of the gene to align.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for minimap2.")
    parser.add_argument("--output", type=Path, help="Output JSON file for alignment results.")
    parser.add_argument("--paf-output", type=Path, help="Output PAF file for raw alignments.")
    parser.add_argument("--extended-paf", type=Path, help="Output PAF file with feature metadata.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")
    args = parser.parse_args()

    setup_logging(args.verbose)

    try:
        # Validate input files exist
        if not args.multi_fasta.exists():
            logging.error(f"Multi-sample FASTA file not found: {args.multi_fasta}")
            sys.exit(1)
        if not args.gff.exists():
            logging.error(f"GFF file not found: {args.gff}")
            sys.exit(1)
        if not args.reference.exists():
            logging.error(f"Reference FASTA file not found: {args.reference}")
            sys.exit(1)

        # 1. Initialize tools - now passing gene identifier to GFFProcessor
        gff_processor = GFFProcessor(args.gff, args.gene)
        seq_extractor = SequenceExtractor(args.reference)
        aligner = HierarchicalAligner(gff_processor, seq_extractor, threads=args.threads)

        # 2. Get the gene feature
        gene_feature = gff_processor.find_gene(args.gene)
        if not gene_feature:
            logging.error(f"Could not find gene '{args.gene}' in the GFF file.")
            sys.exit(1)
        
        # Get the gene name from attributes (try both Name and gene_name)
        gene_name = gene_feature.attributes.get('Name', gene_feature.attributes.get('gene_name', ['-']))
        if isinstance(gene_name, list):
            gene_name = gene_name[0] if gene_name else '-'
        
        logging.info(f"Found gene '{gene_feature.id}' (Name: {gene_name}) at {gene_feature.seqid}:{gene_feature.start}-{gene_feature.end}")

        # 3. Run alignment
        alignment_results = aligner.align_gene_to_haplotypes(gene_feature, args.multi_fasta)

        # 4. Output results
        print_summary(alignment_results)
        
        if args.output:
            save_results(alignment_results, args.output)
        
        # 5. Save raw PAF alignments if requested
        if args.paf_output:
            paf_lines = aligner.get_raw_paf_lines()
            if paf_lines:
                save_paf_alignments(paf_lines, args.paf_output)
            else:
                logging.warning("No PAF alignments to save.")
        
        # 6. Save extended PAF with metadata if requested
        if args.extended_paf:
            paf_lines = aligner.get_raw_paf_lines()
            if paf_lines:
                save_alignments_with_metadata(paf_lines, gff_processor, args.extended_paf)
            else:
                logging.warning("No PAF alignments to save.")

    except FileNotFoundError as e:
        logging.error(f"Error: Input file not found - {e}")
        sys.exit(1)
    except ValueError as e:
        logging.error(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
