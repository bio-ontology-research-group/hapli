#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from collections import defaultdict
import sys

def analyze_paf(paf_file: Path):
    """Analyze PAF file to understand alignment patterns."""
    
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
            
            features_seen.add(query_name)
            targets_seen.add(target_name)
            
            identity = (matches / align_len * 100) if align_len > 0 else 0
            coverage = ((int(parts[3]) - int(parts[2])) / query_len * 100) if query_len > 0 else 0
            
            alignments[target_name][query_name].append({
                'line': line_num,
                'target_start': target_start,
                'target_end': target_end,
                'identity': identity,
                'coverage': coverage,
                'mapq': mapq
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
    args = parser.parse_args()
    
    if not args.paf and not args.json:
        print("Error: Provide at least one of --paf or --json")
        sys.exit(1)
    
    if args.paf and args.paf.exists():
        analyze_paf(args.paf)
    
    if args.json and args.json.exists():
        analyze_json(args.json)

if __name__ == "__main__":
    main()
