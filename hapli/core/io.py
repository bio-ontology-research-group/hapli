import logging
import pysam
from pathlib import Path
from typing import Optional, Iterator, Dict, List, Set, Tuple, Any
from collections import defaultdict
from tqdm import tqdm

class GFFProcessor:
    """Handles GFF parsing and feature hierarchy efficiently for specific targets."""
    def __init__(self, gff_file: Path, target_gene: Optional[str] = None):
        self.gff_file = gff_file
        self.features_by_id = {}
        self.children_map = defaultdict(list)
        self.target_gene = target_gene
        self.logger = logging.getLogger(__name__)
        
        if self.target_gene:
            self._load_gene_features_optimized(self.target_gene)
        else:
            self.logger.warning("No target gene specified. Loading ALL features (this might be slow/memory intensive).")
            self._load_all_features()

    def _parse_attributes(self, attr_string: str) -> Dict[str, List[str]]:
        """Parse GFF attribute string into a dictionary."""
        attributes = defaultdict(list)
        for item in attr_string.strip().split(';'):
            if not item: continue
            if '=' in item:
                key, value = item.split('=', 1)
                # Handle comma-separated values in attributes
                for v in value.split(','):
                    attributes[key].append(v)
        return dict(attributes)
    
    def _load_gene_features_optimized(self, gene_identifier: str):
        """Load only the gene and its related features from the GFF file in a single pass."""
        self.logger.info(f"Searching for gene '{gene_identifier}'...")
        
        gene_found = False
        gene_chr = None
        gene_start = None
        gene_end = None
        
        # Buffer to hold features that might be children but appear before the gene line (unlikely in GFF3 but possible)
        # Actually, simpler: 1st pass find gene coords. 2nd pass load region. 
        # But single pass is better.
        # Strategy:
        # Scan until gene found. Record coords.
        # Scan region around gene.
        # But children usually follow parent or are grouped.
        # Let's assume standard sorted GFF3 (grouped by chrom).
        
        potential_features = []
        
        with open(self.gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                
                attributes = self._parse_attributes(parts[8])
                feature_id = attributes.get('ID', [''])[0]
                feature_type = parts[2]
                
                # Check for gene match
                if not gene_found and feature_type == 'gene':
                    names = [feature_id] + attributes.get('Name', []) + attributes.get('gene_name', [])
                    if gene_identifier in names:
                        gene_found = True
                        gene_chr = parts[0]
                        gene_start = int(parts[3])
                        gene_end = int(parts[4])
                        self._add_feature(parts, attributes, feature_id)
                        self.logger.info(f"Found gene {feature_id} at {gene_chr}:{gene_start}-{gene_end}")
                        continue

                if gene_found:
                    # If we are on same chromosome
                    if parts[0] == gene_chr:
                        start = int(parts[3])
                        end = int(parts[4])
                        
                        # Check overlap with generous buffer (transcripts can be long)
                        # We also check Parent attribute later to confirm relationship
                        if start <= gene_end + 100000 and end >= gene_start - 100000:
                            self._add_feature(parts, attributes, feature_id)
                        elif start > gene_end + 100000:
                            # Assume sorted, we passed it
                            break
                    else:
                        # different chromosome, assume sorted so we are done
                        break
        
        if not gene_found:
            self.logger.error(f"Gene {gene_identifier} not found.")
            return

        # Prune features: Only keep those that are children/descendants of the gene
        # Root is the gene(s) matching identifier
        roots = [f for f in self.features_by_id.values() if f.featuretype == 'gene']
        # Filter strictly by ID match if multiple genes loaded (unlikely with break logic but safe)
        
        valid_ids = set()
        stack = [r.id for r in roots]
        while stack:
            pid = stack.pop()
            valid_ids.add(pid)
            children = self.children_map.get(pid, [])
            stack.extend(children)
            
        # Remove orphans
        all_ids = list(self.features_by_id.keys())
        for fid in all_ids:
            if fid not in valid_ids:
                del self.features_by_id[fid]
                
    def _load_all_features(self):
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
        # We use a simple object instead of gffutils.Feature to avoid dependency complexity if we wanted
        # But sticking to gffutils structure is fine if we mimic it or use it.
        # Since I imported gffutils, let's use a SimpleFeature compatible class or just a namedtuple/dict
        # to avoid strict gffutils dependency if possible? No, gffutils is in requirements.
        # But 'gffutils.Feature' constructor might validate heavily.
        # Let's use a lightweight named tuple or simple class.
        
        feature = Feature(
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
        for parent_id in attributes.get('Parent', []):
            self.children_map[parent_id].append(feature_id)

    def get_children(self, feature_id: str) -> Iterator['Feature']:
        child_ids = self.children_map.get(feature_id, [])
        children = []
        for child_id in child_ids:
            if child_id in self.features_by_id:
                children.append(self.features_by_id[child_id])
        children.sort(key=lambda f: f.start)
        for child in children:
            yield child

class Feature:
    """Simple feature class mimicking gffutils.Feature."""
    def __init__(self, seqid, source, featuretype, start, end, score, strand, frame, attributes, id):
        self.seqid = seqid
        self.source = source
        self.featuretype = featuretype
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attributes = attributes
        self.id = id

class SequenceExtractor:
    """Extracts feature sequences from a reference FASTA."""
    def __init__(self, reference_fasta: Path):
        self.fasta = pysam.FastaFile(str(reference_fasta))
        self.available_chromosomes = set(self.fasta.references)

    def get_sequence(self, feature: Feature) -> Optional[str]:
        if feature.seqid not in self.available_chromosomes:
            return None
        try:
            # pysam fetch is 0-based, [start, end)
            # GFF is 1-based inclusive.
            # So start-1, end
            seq = self.fasta.fetch(feature.seqid, feature.start - 1, feature.end)
        except Exception as e:
            logging.error(f"Error extracting sequence: {e}")
            return None

        if feature.strand == '-':
            return self._reverse_complement(seq)
        return seq

    @staticmethod
    def _reverse_complement(dna: str) -> str:
        complement = str.maketrans('ATCGNRYatcgnry', 'TAGCNRYtagcnry')
        return dna.translate(complement)[::-1]
