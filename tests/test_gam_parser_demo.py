#!/usr/bin/env python3
"""
Demo/test script for GAM parser functionality.

This script demonstrates how to use the GAM parser to load and analyze
GAM alignment data, showing the hierarchical structure and relationships.
"""

import json
import logging
import tempfile
from pathlib import Path
from hapli.gam_parser import GAMParser


def create_test_gam_file() -> Path:
    """
    Create a small test GAM file in JSON format for demonstration.
    
    Returns:
        Path to the temporary GAM file
    """
    # Sample GAM alignments in JSON format (what vg view -a would output)
    test_alignments = [
        {
            "name": "gene_001",
            "sequence": "ATCGATCGATCGATCG",
            "quality": "IIIIIIIIIIIIIIII",
            "path": {
                "name": "sample1#0#chr1",
                "mapping_quality": 60
            },
            "score": 95,
            "identity": 0.98
        },
        {
            "name": "exon_001_01",
            "sequence": "ATCGATCG",
            "quality": "IIIIIIII",
            "path": {
                "name": "sample1#0#chr1",
                "mapping_quality": 58
            },
            "score": 88,
            "identity": 0.96
        },
        {
            "name": "cds_001_01",
            "sequence": "ATCGAT",
            "quality": "IIIIII",
            "path": {
                "name": "sample1#0#chr1",
                "mapping_quality": 55
            },
            "score": 82,
            "identity": 0.94
        },
        {
            "name": "gene_002",
            "sequence": "GCTAGCTAGCTAGCTA",
            "quality": "IIIIIIIIIIIIIIII",
            "path": {
                "name": "sample1#1#chr1",
                "mapping_quality": 62
            },
            "score": 92,
            "identity": 0.97
        },
        {
            "name": "transcript_002",
            "sequence": "GCTAGCTAGCTA",
            "quality": "IIIIIIIIIIII",
            "path": {
                "name": "sample1#1#chr1",
                "mapping_quality": 60
            },
            "score": 89,
            "identity": 0.96
        },
        {
            "name": "gene_003",
            "sequence": "TTAAGGCCTTAAGGCC",
            "quality": "IIIIIIIIIIIIIIII",
            "path": {
                "name": "sample2.0.chr2",
                "mapping_quality": 61
            },
            "score": 94,
            "identity": 0.99
        },
        {
            "name": "exon_003_01",
            "sequence": "TTAAGGCC",
            "quality": "IIIIIIII",
            "path": {
                "name": "sample2.0.chr2",
                "mapping_quality": 59
            },
            "score": 87,
            "identity": 0.95
        },
        {
            "name": "utr_003_5prime",
            "sequence": "TTAA",
            "quality": "IIII",
            "path": {
                "name": "sample2.0.chr2",
                "mapping_quality": 57
            },
            "score": 78,
            "identity": 0.92
        }
    ]
    
    # Create temporary file
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False)
    
    # Write alignments as JSON lines (GAM format when converted to JSON)
    for alignment in test_alignments:
        json.dump(alignment, temp_file)
        temp_file.write('\n')
    
    temp_file.close()
    return Path(temp_file.name)


def analyze_sample_haplotype_data(sample: str, haplotype: str, data: dict) -> dict:
    """
    Analyze data for a specific sample/haplotype combination.
    
    Args:
        sample: Sample name
        haplotype: Haplotype identifier
        data: Feature data for this sample/haplotype
        
    Returns:
        Analysis summary
    """
    analysis = {
        'sample': sample,
        'haplotype': haplotype,
        'feature_types': list(data.keys()),
        'feature_counts': {},
        'top_features': {},
        'example_relationships': []
    }
    
    # Count features by type and find top features
    for feature_type, features in data.items():
        analysis['feature_counts'][feature_type] = len(features)
        
        # Sort by score and get top feature
        sorted_features = sorted(features, key=lambda x: x['score'], reverse=True)
        if sorted_features:
            top_feature = sorted_features[0]
            analysis['top_features'][feature_type] = {
                'name': top_feature['read_name'],
                'score': top_feature['score'],
                'identity': top_feature['identity']
            }
    
    # Find example parent-child relationships
    # Look for features that could be hierarchically related based on naming
    for feature_type, features in data.items():
        for feature in features:
            read_name = feature['read_name']
            
            # Look for potential children (features with similar base names)
            base_name = read_name.split('_')[0] if '_' in read_name else read_name
            
            for other_type, other_features in data.items():
                if other_type != feature_type:
                    for other_feature in other_features:
                        other_name = other_feature['read_name']
                        
                        # Check if names suggest parent-child relationship
                        if (base_name in other_name and 
                            read_name != other_name and
                            len(analysis['example_relationships']) < 2):
                            
                            relationship = {
                                'potential_parent': {
                                    'name': read_name,
                                    'type': feature_type,
                                    'score': feature['score']
                                },
                                'potential_child': {
                                    'name': other_name,
                                    'type': other_type,
                                    'score': other_feature['score']
                                }
                            }
                            analysis['example_relationships'].append(relationship)
    
    return analysis


def main():
    """Main demonstration function."""
    print("=" * 60)
    print("GAM Parser Demo")
    print("=" * 60)
    
    # Setup logging
    logging.basicConfig(level=logging.INFO)
    
    try:
        # Create test GAM file
        print("\n1. Creating test GAM data...")
        test_gam_path = create_test_gam_file()
        print(f"   Created test file: {test_gam_path}")
        
        # Initialize parser (mock vg command since we have JSON directly)
        print("\n2. Initializing GAM parser...")
        
        # For this demo, we'll directly load the JSON we created
        # In real usage, the parser would call vg view -a on a GAM file
        parser = GAMParser.__new__(GAMParser)  # Create without __init__
        parser.gam_file = test_gam_path
        parser.vg_executable = "vg"
        parser.alignments = []
        
        # Load the JSON alignments directly
        with open(test_gam_path, 'r') as f:
            for line in f:
                if line.strip():
                    alignment = json.loads(line)
                    parser.alignments.append(alignment)
        
        print(f"   Loaded {len(parser.alignments)} alignments")
        
        # Group alignments by sample/haplotype
        print("\n3. Grouping alignments by sample/haplotype...")
        grouped_data = parser.group_alignments_by_sample_haplotype()
        
        # Pretty-print the full nested structure
        print("\n4. Complete nested structure:")
        print("-" * 40)
        print(json.dumps(grouped_data, indent=2, default=str))
        
        # Analyze each sample/haplotype combination
        print("\n5. Sample/Haplotype Analysis:")
        print("-" * 40)
        
        for sample in sorted(grouped_data.keys()):
            for haplotype in sorted(grouped_data[sample].keys()):
                analysis = analyze_sample_haplotype_data(
                    sample, haplotype, grouped_data[sample][haplotype]
                )
                
                print(f"\nSample: {sample}, Haplotype: {haplotype}")
                print(f"  Feature types: {', '.join(analysis['feature_types'])}")
                print(f"  Feature counts: {analysis['feature_counts']}")
                
                print("  Top features by type:")
                for ftype, top_feature in analysis['top_features'].items():
                    print(f"    {ftype}: {top_feature['name']} "
                          f"(score: {top_feature['score']}, "
                          f"identity: {top_feature['identity']:.3f})")
                
                if analysis['example_relationships']:
                    print("  Example parent-child relationships:")
                    for rel in analysis['example_relationships']:
                        parent = rel['potential_parent']
                        child = rel['potential_child']
                        print(f"    {parent['name']} ({parent['type']}) -> "
                              f"{child['name']} ({child['type']})")
                else:
                    print("  No clear parent-child relationships detected")
        
        # Get and display statistics
        print("\n6. Alignment Statistics:")
        print("-" * 40)
        stats = parser.get_alignment_statistics()
        print(json.dumps(stats, indent=2, default=str))
        
        print("\n7. Summary:")
        print("-" * 40)
        print(f"Total alignments: {stats['total_alignments']}")
        print(f"Samples found: {', '.join(stats['samples'])}")
        print(f"Haplotypes found: {', '.join(stats['haplotypes'])}")
        print(f"Feature types: {', '.join(stats['feature_types'].keys())}")
        print(f"Average score: {stats['avg_score']:.1f}")
        print(f"Average identity: {stats['avg_identity']:.3f}")
        
    finally:
        # Clean up temporary file
        if 'test_gam_path' in locals():
            test_gam_path.unlink(missing_ok=True)
            print(f"\nCleaned up temporary file: {test_gam_path}")
    
    print("\n" + "=" * 60)
    print("Demo completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
