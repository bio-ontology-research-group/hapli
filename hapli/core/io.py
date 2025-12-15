import logging
import gffutils
import pysam
from pathlib import Path
from typing import Optional, Iterator, Dict, List, Set, Tuple
from collections import defaultdict
from tqdm import tqdm

class GFFProcessor:
    """Handles GFF parsing and feature hierarchy."""
    def __init__(self, gff_file: Path, target_region: Optional[Tuple[str, int, int]] = None, target_gene: Optional[str] = None):
        self.gff_file = gff_file
        self.features_by_id = {}
        self.children_map = defaultdict(list)
        self.target_gene = target_gene
        self.target_region = target_region
        
        # Initialize database or load features
        self._load_features()

    def _parse_attributes(self, attr_string: str) -> Dict[str, List[str]]:
        """Parse GFF attribute string into a dictionary."""
        attributes = defaultdict(list)
        for item in attr_string.strip().split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attributes[key].append(value)
        return dict(attributes)
    
    def _load_features(self):
        """Load features from the GFF file."""
        logging.info(f"Loading features from {self.gff_file}")
        
        # If we have a specific gene target, optimize
        if self.target_gene:
            self._load_gene_features_optimized(self.target_gene)
            return

        # Otherwise load everything (could be slow for whole genome)
        # For this prototype, we'll assume we might want to load a region if specified
        # or we rely on the caller to provide a filtered GFF or use GFFUtils DB which is slower to create but faster to query.
        
        # Given the "rethink" scope, let's stick to the efficient single-pass parser 
        # but generalize it slightly or keep the gene-centric optimization as it's a key requirement.
        
        # Fallback: simple full load if no target specified (careful with large files)
        self._load_all_features()

    def _load_gene_features_optimized(self, gene_identifier: str):
        """Load only the gene and its related features from the GFF file in a single pass."""
        logging.info(f"Searching for gene '{gene_identifier}' and loading its features...")
        
        gene_found = False
        gene_chr = None
        gene_start = None
        gene_end = None
        all_features = []
        parent_to_children = defaultdict(set)
        
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
                        feature_name = attributes.get('Name', [''])[0]
                        if not feature_name:
                            feature_name = attributes.get('gene_name', [''])[0]
                        gene_id_attr = attributes.get('gene_id', [''])[0]
                        
                        if (feature_id == gene_identifier or 
                            feature_name == gene_identifier or
                            gene_id_attr == gene_identifier):
                            
                            self._add_feature(parts, attributes, feature_id)
                            gene_found = True
                            gene_chr = parts[0]
                            gene_start = int(parts[3])
                            gene_end = int(parts[4])
                            logging.info(f"Found gene '{gene_identifier}' at {gene_chr}:{gene_start}-{gene_end}")
                    
                    # If gene is found, collect features in the same region (optimization)
                    if gene_found:
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
                            
                            # Optimization: Stop scanning if we are far past the gene
                            elif feat_start > gene_end + 50000:
                                logging.info(f"Passed gene region ({feat_start} > {gene_end} + buffer). Stopping scan.")
                                break
                        elif parts[0] != gene_chr:
                             # We moved to another chromosome after finding the gene
                             logging.info(f"Moved from {gene_chr} to {parts[0]}. Stopping scan.")
                             break
        
        if not gene_found:
            raise ValueError(f"Gene '{gene_identifier}' not found in GFF file")
        
        # Traverse hierarchy
        relevant_ids = {self.target_gene} # This might be wrong if target_gene is Name not ID.
        # But wait, self.features_by_id is populated in _add_feature. 
        # We need to find the ID of the gene we found.
        # Let's fix the logic: we store the gene feature in _add_feature, so we can find it.
        
        # Actually, let's keep it simple: filter all_features based on the gene we found.
        # We need the ID of the gene.
        gene_obj = None
        for fid, feat in self.features_by_id.items():
            if feat.featuretype == 'gene':
                gene_obj = feat
                break
        
        if not gene_obj:
             # Should not happen if gene_found is True and _add_feature was called
             pass
        else:
             relevant_ids = {gene_obj.id}
             
        to_process = list(relevant_ids)
        processed = set(relevant_ids)

        while to_process:
            current_id = to_process.pop()
            for child_id in parent_to_children.get(current_id, []):
                if child_id not in processed:
                    processed.add(child_id)
                    relevant_ids.add(child_id)
                    to_process.append(child_id)
        
        # Now add relevant features to self.features_by_id
        self.features_by_id = {} # Clear initial
        self.children_map = defaultdict(list)
        
        # Re-add gene
        if gene_obj:
             self.features_by_id[gene_obj.id] = gene_obj
             
        for parts, attributes, feature_id in all_features:
            if feature_id in relevant_ids and feature_id not in self.features_by_id:
                self._add_feature(parts, attributes, feature_id)

    def _load_all_features(self):
        """Load all features (fallback)."""
        logging.info("Loading all features...")
        with open(self.gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                attributes = self._parse_attributes(parts[8])
                feature_id = attributes.get('ID', [''])[0]
                if feature_id:
                    self._add_feature(parts, attributes, feature_id)

    def _add_feature(self, parts, attributes, feature_id):
        feature = gffutils.Feature(
            seqid=parts[0], source=parts[1], featuretype=parts[2],
            start=int(parts[3]), end=int(parts[4]), score=parts[5],
            strand=parts[6], frame=parts[7], attributes=attributes, id=feature_id
        )
        self.features_by_id[feature_id] = feature
        for parent_id in attributes.get('Parent', []):
            self.children_map[parent_id].append(feature_id)

    def get_children(self, feature: gffutils.Feature) -> Iterator[gffutils.Feature]:
        child_ids = self.children_map.get(feature.id, [])
        children = []
        for child_id in child_ids:
            if child_id in self.features_by_id:
                children.append(self.features_by_id[child_id])
        children.sort(key=lambda f: f.start)
        for child in children:
            yield child

class SequenceExtractor:
    """Extracts feature sequences from a reference FASTA."""
    def __init__(self, reference_fasta: Path):
        self.fasta = pysam.FastaFile(str(reference_fasta))
        self.available_chromosomes = set(self.fasta.references)

    def get_sequence(self, feature: gffutils.Feature) -> Optional[str]:
        if feature.seqid not in self.available_chromosomes:
            return None
        try:
            seq = self.fasta.fetch(feature.seqid, feature.start - 1, feature.end)
        except Exception:
            return None

        if feature.strand == '-':
            return self._reverse_complement(seq)
        return seq

    @staticmethod
    def _reverse_complement(dna: str) -> str:
        complement = str.maketrans('ATCGNRY', 'TAGCNYR')
        return dna.upper().translate(complement)[::-1]
