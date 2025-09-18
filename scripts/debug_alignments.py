#!/usr/bin/env python3

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pysam

# Color codes for terminal output
COLORS = {
    "RED": '\033[91m',
    "GREEN": '\033[92m',
    "YELLOW": '\033[93m',
    "BLUE": '\033[94m',
    "ENDC": '\033[0m'
}

def _reverse_complement(dna: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    complement = str.maketrans('ATCGNRY', 'TAGCNYR')
    return dna.upper().translate(complement)[::-1]

def get_alignment_visualization(feature: dict, aln: dict, ref_genome: pysam.FastaFile, path_seqs: pysam.FastaFile, use_color: bool) -> str:
    """Generates a colorized string showing the base-level alignment."""
    colors = COLORS if use_color else {k: '' for k in COLORS}
    if not aln or not aln.get('cigar'):
        return "  (No CIGAR string available for visualization)"

    # Get sequences
    try:
        ref_seq = ref_genome.fetch(feature['chrom'], feature['start'] - 1, feature['end'])
        if feature['strand'] == '-':
            ref_seq = _reverse_complement(ref_seq)
    except Exception as e:
        return f"  (Could not fetch reference sequence: {e})"

    try:
        target_seq = path_seqs.fetch(aln['target_name'], aln['target_start'], aln['target_end'])
    except Exception as e:
        return f"  (Could not fetch target sequence: {e})"

    # Parse CIGAR and build visualization
    cigar_tuples = re.findall(r'(\d+)([MIDNSHPX=])', aln['cigar'])
    query_pos, target_pos = 0, 0
    query_line, match_line, target_line = "", "", ""

    for length_str, op in cigar_tuples:
        length = int(length_str)
        if op in ['M', '=', 'X']:
            for _ in range(length):
                if query_pos < len(ref_seq) and target_pos < len(target_seq):
                    q_base, t_base = ref_seq[query_pos], target_seq[target_pos]
                    if q_base.upper() == t_base.upper():
                        query_line += q_base
                        match_line += "|"
                        target_line += t_base
                    else:
                        query_line += f"{colors['RED']}{q_base}{colors['ENDC']}"
                        match_line += f"{colors['YELLOW']}X{colors['ENDC']}"
                        target_line += f"{colors['RED']}{t_base}{colors['ENDC']}"
                query_pos += 1
                target_pos += 1
        elif op == 'I':
            for _ in range(length):
                if query_pos < len(ref_seq):
                    query_line += f"{colors['GREEN']}{ref_seq[query_pos]}{colors['ENDC']}"
                    match_line += f"{colors['GREEN']}+{colors['ENDC']}"
                    target_line += f"{colors['GREEN']}-{colors['ENDC']}"
                    query_pos += 1
        elif op == 'D':
            for _ in range(length):
                if target_pos < len(target_seq):
                    query_line += f"{colors['BLUE']}-{colors['ENDC']}"
                    match_line += f"{colors['BLUE']}-{colors['ENDC']}"
                    target_line += f"{colors['BLUE']}{target_seq[target_pos]}{colors['ENDC']}"
                    target_pos += 1
        elif op == 'S':
            query_pos += length

    # Format output in chunks
    chunk_size = 80
    result_lines = []
    for i in range(0, len(query_line), chunk_size):
        chunk_end = min(i + chunk_size, len(query_line))
        result_lines.append(f"        Ref:    {query_line[i:chunk_end]}")
        result_lines.append(f"                {' ' * len('Ref:    ')}{match_line[i:chunk_end]}")
        result_lines.append(f"        Target: {target_line[i:chunk_end]}")
        result_lines.append("")
    
    return "\n".join(result_lines)


def analyze_paf(paf_file: Path, gff_file: Path, reference_fasta: Path, path_fasta: Path):
    """Analyze PAF file to understand alignment patterns."""
    
    # Load GFF features for visualization
    gff_features = {}
    if gff_file and gff_file.exists():
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                attributes = {k: v for k, v in (item.split('=', 1) for item in parts[8].split(';') if '=' in item)}
                feature_id = attributes.get("ID")
                if feature_id:
                    gff_features[feature_id] = {
                        "id": feature_id, "type": parts[2], "chrom": parts[0],
                        "start": int(parts[3]), "end": int(parts[4]), "strand": parts[6]
                    }

    # Load FASTA files for visualization
    ref_genome, path_seqs = None, None
    can_visualize = gff_file and reference_fasta and reference_fasta.exists() and path_fasta and path_fasta.exists()
    if can_visualize:
        ref_genome = pysam.FastaFile(str(reference_fasta))
        path_seqs = pysam.FastaFile(str(path_fasta))

    # Track alignments per feature per target
    alignments = defaultdict(lambda: defaultdict(list))
    features_seen = set()
    targets_seen = set()
    
    print(f"\nAnalyzing PAF file: {paf_file}")
    print("="*80)
    
    with open(paf_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            
            query_name = parts[0]
            query_len = int(parts[1])
            target_name = parts[5]
            target_start = int(parts[7])
            target_end = int(parts[8])
            matches = int(parts[9])
            align_len = int(parts[10])
            mapq = int(parts[11])
            
            cigar = None
            for tag in parts[12:]:
                if tag.startswith("cg:Z:"):
                    cigar = tag[5:]
                    break

            features_seen.add(query_name)
            targets_seen.add(target_name)
            
            identity = (matches / align_len * 100) if align_len > 0 else 0
            coverage = ((int(parts[3]) - int(parts[2])) / query_len * 100) if query_len > 0 else 0
            
            alignments[target_name][query_name].append({
                'line': line_num,
                'target_name': target_name,
                'target_start': target_start,
                'target_end': target_end,
                'identity': identity,
                'coverage': coverage,
                'mapq': mapq,
                'cigar': cigar
            })
    
    print(f"Total features: {len(features_seen)}")
    print(f"Total targets: {len(targets_seen)}")
    print()
    
    # Analyze each target
    for target in sorted(targets_seen):
        is_reference = "grch38" in target.lower() or "hg38" in target.lower()
        
        print(f"\nTarget: {target}")
        if is_reference:
            print("  ** THIS IS A REFERENCE PATH **")
        
        print(f"  Features aligned: {len(alignments[target])}")
        
        # Check for duplications
        duplicated_features = []
        missing_critical = []
        
        for feature, alns in alignments[target].items():
            if len(alns) > 1:
                duplicated_features.append((feature, len(alns)))
        
        if duplicated_features:
            print(f"  DUPLICATED features ({len(duplicated_features)}):")
            for feat, count in sorted(duplicated_features)[:10]:  # Show first 10
                print(f"    - {feat}: {count} alignments")
                if is_reference:
                    print(f"      WARNING: Reference should not have duplications!")
                
                # Visualize the different alignments for duplicated features
                if can_visualize and feat in gff_features:
                    feature_data = gff_features[feat]
                    for i, aln in enumerate(alignments[target][feat]):
                        print(f"      Alignment #{i+1} (MAPQ: {aln['mapq']}, Target: {aln['target_start']}-{aln['target_end']}):")
                        viz = get_alignment_visualization(feature_data, aln, ref_genome, path_seqs, True)
                        print(viz)
        
        # Check for critical missing features
        critical_patterns = ['start_codon', 'stop_codon', 'UTR', 'exon']
        aligned_features = set(alignments[target].keys())
        
        for feature in features_seen:
            if any(pattern in feature for pattern in critical_patterns):
                if feature not in aligned_features:
                    missing_critical.append(feature)
        
        if missing_critical:
            print(f"  MISSING critical features ({len(missing_critical)}):")
            for feat in sorted(missing_critical)[:10]:  # Show first 10
                print(f"    - {feat}")
                if is_reference:
                    print(f"      ERROR: Reference should have ALL features!")
        
        # Show alignment quality distribution
        quality_stats = defaultdict(int)
        for feature, alns in alignments[target].items():
            for aln in alns:
                if aln['identity'] == 100 and aln['coverage'] == 100:
                    quality_stats['perfect'] += 1
                elif aln['identity'] >= 95:
                    quality_stats['high'] += 1
                elif aln['identity'] >= 80:
                    quality_stats['medium'] += 1
                else:
                    quality_stats['low'] += 1
        
        print(f"  Alignment quality:")
        print(f"    Perfect (100%): {quality_stats['perfect']}")
        print(f"    High (95-99%): {quality_stats['high']}")
        print(f"    Medium (80-94%): {quality_stats['medium']}")
        print(f"    Low (<80%): {quality_stats['low']}")

def analyze_json(json_file: Path):
    """Analyze JSON output from hierarchical aligner."""
    
    print(f"\nAnalyzing JSON file: {json_file}")
    print("="*80)
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    for target, gene_data in data.items():
        is_reference = "grch38" in target.lower() or "hg38" in target.lower()
        
        print(f"\nTarget: {target}")
        if is_reference:
            print("  ** THIS IS A REFERENCE PATH **")
        
        # Recursively count features
        def count_features(node, depth=0):
            counts = {'total': 0, 'perfect': 0, 'missing': 0}
            
            if node.get('identity', 0) >= 0.999:
                counts['perfect'] += 1
            elif node.get('identity', 0) == 0:
                counts['missing'] += 1
            counts['total'] += 1
            
            for child in node.get('children', []):
                child_counts = count_features(child, depth + 1)
                for key in counts:
                    counts[key] += child_counts[key]
            
            return counts
        
        counts = count_features(gene_data)
        print(f"  Total features: {counts['total']}")
        print(f"  Perfect alignments: {counts['perfect']}")
        print(f"  Missing features: {counts['missing']}")
        
        if is_reference and counts['missing'] > 0:
            print(f"  ERROR: Reference has {counts['missing']} missing features!")
        
        # Show tree structure
        def print_tree(node, indent=0):
            prefix = "  " * (indent + 1)
            identity = node.get('identity', 0) * 100
            print(f"{prefix}- {node['feature_type']} {node['feature_id']}: {identity:.1f}% identity")
            
            if is_reference and identity < 100 and indent > 0:
                print(f"{prefix}  WARNING: Reference feature not 100%!")
            
            for child in node.get('children', []):
                print_tree(child, indent + 1)
        
        print("  Feature tree:")
        print_tree(gene_data)

def main():
    parser = argparse.ArgumentParser(description="Debug alignment issues")
    parser.add_argument("--paf", type=Path, help="PAF file to analyze")
    parser.add_argument("--json", type=Path, help="JSON file from hierarchical aligner")
    parser.add_argument("--gff", type=Path, help="GFF3 file (for visualization)")
    parser.add_argument("--reference-fasta", type=Path, help="Reference FASTA (for visualization)")
    parser.add_argument("--path-fasta", type=Path, help="Path FASTA (for visualization)")
    args = parser.parse_args()
    
    if not args.paf and not args.json:
        print("Error: Provide at least one of --paf or --json")
        sys.exit(1)
    
    if args.paf and args.paf.exists():
        analyze_paf(args.paf, args.gff, args.reference_fasta, args.path_fasta)
    
    if args.json and args.json.exists():
        analyze_json(args.json)

if __name__ == "__main__":
    main()
