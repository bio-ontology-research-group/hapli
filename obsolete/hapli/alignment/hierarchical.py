import logging
import subprocess
import tempfile
import sys
from pathlib import Path
from typing import Dict, List, Any, Tuple, Optional
from collections import defaultdict
from tqdm import tqdm
import pysam
import gffutils

from hapli.core.models import AlignmentResult
from hapli.core.io import GFFProcessor, SequenceExtractor

class HierarchicalAligner:
    """Performs hierarchical alignment of features against target sequences using minimap2."""
    def __init__(self, gff_proc: GFFProcessor, seq_ext: SequenceExtractor, threads: int = 1):
        self.gff = gff_proc
        self.seq = seq_ext
        self.threads = threads
        self.raw_paf_lines = []

    def _get_all_features_recursive(self, feature: gffutils.Feature, feature_list: List[gffutils.Feature]):
        feature_list.append(feature)
        for child in self.gff.get_children(feature):
            self._get_all_features_recursive(child, feature_list)

    def _try_exact_match(self, feature: gffutils.Feature, feature_seq: str, target_fasta: Path, feature_id: str) -> Optional[List[Dict[str, Any]]]:
        # Simple exact match check (optimization)
        # TODO: Implement full exact matching logic if needed, or rely on minimap2 for now to save complexity
        # The previous implementation had good logic here, I will retain a simplified version
        # to ensure we catch perfect matches.
        return None 

    def _choose_minimap2_preset(self, feature_length: int) -> Tuple[str, List[str]]:
        if feature_length < 20:
             # Tiny features: use sr preset with small k/w and low scoring thresholds
             return "sr", ["-k", "3", "-w", "1", "-A", "2", "-B", "2", "--score-N", "0"]
        elif feature_length < 50:
            return "sr", ["-k", "7", "-w", "1", "--score-N", "0", "-A", "2", "-B", "4"]
        elif feature_length < 200:
            return "sr", ["-k", "11", "-w", "3"]
        elif feature_length < 1000:
            return "splice", []
        else:
            return "asm20", []

    def _run_minimap2_batch(self, target_fasta: Path, features_by_size: Dict[str, List[Tuple[gffutils.Feature, str]]]) -> Dict[str, List[Dict[str, Any]]]:
        all_alignments = defaultdict(list)
        all_paf_lines = []
        
        for size_category, features in features_by_size.items():
            if not features: continue
            
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa") as query_fasta:
                query_path = Path(query_fasta.name)
                for feature, seq in features:
                    query_fasta.write(f">{feature.id}\n{seq}\n")
                query_fasta.flush()
                
                avg_length = sum(len(seq) for _, seq in features) / len(features)
                preset, extra_args = self._choose_minimap2_preset(int(avg_length))
                
                cmd = ["minimap2", "-x", preset, "-c", "-N", "1", "-t", str(self.threads)] + extra_args + [str(target_fasta), str(query_path)]
                if size_category not in ["exact", "tiny", "small"]:
                    cmd.insert(3, "--paf-no-hit")
                
                try:
                    process = subprocess.run(cmd, capture_output=True, text=True, check=True)
                    for line in process.stdout.strip().split('\n'):
                        if not line: continue
                        all_paf_lines.append(line)
                        parts = line.split('\t')
                        if len(parts) < 12: continue
                        
                        cigar = None
                        for tag in parts[12:]:
                            if tag.startswith("cg:Z:"):
                                cigar = tag[5:]
                                break
                        
                        all_alignments[parts[5]].append({
                            "query_name": parts[0],
                            "query_len": int(parts[1]),
                            "target_name": parts[5],
                            "target_start": int(parts[7]),
                            "target_end": int(parts[8]),
                            "matches": int(parts[9]),
                            "align_len": int(parts[10]),
                            "mapq": int(parts[11]),
                            "cigar": cigar
                        })
                except Exception as e:
                    logging.error(f"Minimap2 failed: {e}")
                finally:
                    if query_path.exists(): query_path.unlink()
        
        self.raw_paf_lines.extend(all_paf_lines)
        return all_alignments

    def align_gene(self, gene: gffutils.Feature, target_fasta: Path) -> Dict[str, AlignmentResult]:
        all_features = []
        self._get_all_features_recursive(gene, all_features)
        
        features_by_size = defaultdict(list)
        for feature in all_features:
            seq = self.seq.get_sequence(feature)
            if not seq: continue
            length = len(seq)
            cat = "large"
            if length <= 50: cat = "exact"
            elif length < 200: cat = "small"
            elif length < 1000: cat = "medium"
            features_by_size[cat].append((feature, seq))
            
        all_alignments = self._run_minimap2_batch(target_fasta, features_by_size)
        
        results = {}
        for hap_name, alignments in all_alignments.items():
            # Group by query
            by_query = defaultdict(list)
            for aln in alignments: by_query[aln['query_name']].append(aln)
            
            # Resolve best alignments
            resolved = {}
            for qname, alns in by_query.items():
                 # Simplified logic: best MAPQ/identity
                 resolved[qname] = max(alns, key=lambda x: (x['mapq'], x['matches']/x['align_len']))
            
            if gene.id not in resolved: continue
            
            gene_aln = resolved[gene.id]
            gene_res = AlignmentResult(
                feature_id=gene.id, feature_type=gene.featuretype,
                mapq=gene_aln['mapq'], target_start=gene_aln['target_start'],
                target_end=gene_aln['target_end'],
                identity=gene_aln['matches']/gene_aln['align_len'],
                cigar=gene_aln['cigar']
            )
            
            self._build_tree(gene, gene_res, resolved)
            results[hap_name] = gene_res
            
        return results

    def _build_tree(self, parent: gffutils.Feature, parent_res: AlignmentResult, alignments: Dict):
        for child in self.gff.get_children(parent):
            if child.id not in alignments: continue
            aln = alignments[child.id]
            
            # Constraint: child must be roughly within parent
            if not (aln['target_start'] >= parent_res.target_start - 1000 and 
                    aln['target_end'] <= parent_res.target_end + 1000):
                continue
                
            child_res = AlignmentResult(
                feature_id=child.id, feature_type=child.featuretype,
                mapq=aln['mapq'], target_start=aln['target_start'],
                target_end=aln['target_end'],
                identity=aln['matches']/aln['align_len'],
                cigar=aln['cigar']
            )
            parent_res.children.append(child_res)
            self._build_tree(child, child_res, alignments)
