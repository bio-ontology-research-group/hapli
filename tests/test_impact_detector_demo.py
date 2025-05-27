#!/usr/bin/env python3
"""
Demo/test script for ImpactDetector functionality.

This script demonstrates how to use the ImpactDetector to analyze structural
impacts on genomic features, showing comprehensive impact analysis with
colored output and detailed explanations.
"""

import json
import logging
import tempfile
from pathlib import Path
from hapli.impact_detector import ImpactDetector

# ANSI color codes for terminal output
class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Impact type color mapping
IMPACT_COLORS = {
    'INTACT': Colors.OKGREEN,
    'TRUNCATED': Colors.WARNING,
    'SPLIT': Colors.FAIL,
    'MISSING': Colors.FAIL
}

STRUCTURAL_COLORS = {
    'INVERSION': Colors.FAIL,
    'DUPLICATION': Colors.WARNING,
    'DELETION': Colors.FAIL,
    'TRANSLOCATION': Colors.FAIL,
    'COPY_NUMBER_CHANGE': Colors.WARNING,
    'COMPLEX_REARRANGEMENT': Colors.FAIL
}

def create_test_gfa_file() -> Path:
    """
    Create a small test GFA file for demonstration.
    
    Returns:
        Path to the temporary GFA file
    """
    gfa_content = """H	VN:Z:1.0
S	1	ATCGATCGATCGATCG
S	2	GCTAGCTAGCTAGCTA
S	3	TTAAGGCCTTAAGGCC
S	4	CCGGAATTCCGGAATT
S	5	AAAATTTTAAAATTTT
S	6	GGGGCCCCGGGGCCCC
S	7	TACGTACGTACGTACG
S	8	CGTACGTACGTACGTA
L	1	+	2	+	4M
L	2	+	3	+	4M
L	3	+	4	+	4M
L	4	+	5	+	4M
L	1	+	6	+	4M
L	6	+	7	+	4M
L	7	+	8	+	4M
P	reference	1+,2+,3+,4+,5+	*
P	sample1#0#chr1	1+,2+,3+,4+,5+	*
P	sample1#1#chr1	1+,6+,7+,8+	*
P	sample2#0#chr1	1+,2-,3+,4+,5+	*
"""
    
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.gfa', delete=False)
    temp_file.write(gfa_content)
    temp_file.close()
    return Path(temp_file.name)

def create_test_gam_data() -> dict:
    """
    Create test GAM data that corresponds to the GFA structure.
    
    Returns:
        GAM data in the format expected by ImpactDetector
    """
    return {
        "sample1": {
            "0": {
                "gene": [
                    {
                        "read_name": "gene_001",
                        "sequence": "ATCGATCGATCGATCGGCTAGCTAGCTAGCTA",
                        "quality": "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                        "feature_type": "gene",
                        "score": 95,
                        "identity": 0.98,
                        "path_positions": [
                            {"node_id": "1", "offset": 0, "is_reverse": False},
                            {"node_id": "2", "offset": 0, "is_reverse": False}
                        ]
                    }
                ],
                "CDS": [
                    {
                        "read_name": "cds_001_01",
                        "sequence": "ATGGATCGATCGATAG",
                        "quality": "IIIIIIIIIIIIIII",
                        "feature_type": "CDS",
                        "score": 88,
                        "identity": 0.96,
                        "path_positions": [
                            {"node_id": "1", "offset": 0, "is_reverse": False}
                        ]
                    }
                ],
                "exon": [
                    {
                        "read_name": "exon_001_01",
                        "sequence": "ATCGATCG",
                        "quality": "IIIIIIII",
                        "feature_type": "exon",
                        "score": 82,
                        "identity": 0.94,
                        "path_positions": [
                            {"node_id": "1", "offset": 0, "is_reverse": False}
                        ]
                    }
                ]
            },
            "1": {
                "gene": [
                    {
                        "read_name": "gene_002",
                        "sequence": "ATCGATCGATCGATCGGGGGCCCCGGGGCCCC",
                        "quality": "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                        "feature_type": "gene",
                        "score": 92,
                        "identity": 0.97,
                        "path_positions": [
                            {"node_id": "1", "offset": 0, "is_reverse": False},
                            {"node_id": "6", "offset": 0, "is_reverse": False}
                        ]
                    }
                ],
                "promoter": [
                    {
                        "read_name": "promoter_002",
                        "sequence": "GGGGCCCC",
                        "quality": "IIIIIIII",
                        "feature_type": "promoter",
                        "score": 75,
                        "identity": 0.85,
                        "path_positions": [
                            {"node_id": "6", "offset": 0, "is_reverse": False}
                        ]
                    }
                ]
            }
        },
        "sample2": {
            "0": {
                "gene": [
                    {
                        "read_name": "gene_003",
                        "sequence": "ATCGATCGATCGATCGTAGCTAGCTAGCTATT",
                        "quality": "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                        "feature_type": "gene",
                        "score": 70,
                        "identity": 0.89,
                        "path_positions": [
                            {"node_id": "1", "offset": 0, "is_reverse": False},
                            {"node_id": "2", "offset": 0, "is_reverse": True}  # Inverted
                        ]
                    }
                ],
                "CDS": [
                    {
                        "read_name": "cds_003_01",
                        "sequence": "ATGGATCGATCGTAG",
                        "quality": "IIIIIIIIIIIIIII",
                        "feature_type": "CDS",
                        "score": 65,
                        "identity": 0.87,
                        "path_positions": [
                            {"node_id": "1", "offset": 0, "is_reverse": False},
                            {"node_id": "2", "offset": 0, "is_reverse": True}  # Inverted
                        ]
                    }
                ],
                "UTR": [
                    {
                        "read_name": "utr_003_5prime",
                        "sequence": "TTAA",
                        "quality": "IIII",
                        "feature_type": "UTR",
                        "score": 45,
                        "identity": 0.70,
                        "path_positions": [
                            {"node_id": "3", "offset": 0, "is_reverse": False}
                        ]
                    }
                ]
            }
        }
    }

def get_impact_explanation(feature_name: str, feature_type: str, impact_type: str, 
                         consequence: str, structural_impacts: list, 
                         details: dict) -> str:
    """
    Generate a human-readable explanation for the impact.
    
    Args:
        feature_name: Name of the feature
        feature_type: Type of the feature
        impact_type: General impact type
        consequence: Specific consequence
        structural_impacts: List of structural impacts
        details: Impact analysis details
        
    Returns:
        Human-readable explanation string
    """
    explanations = []
    
    # Coverage information
    coverage = details.get('coverage', 0.0)
    coverage_pct = coverage * 100
    
    # Basic impact explanation
    if impact_type == 'INTACT':
        explanations.append(f"fully preserved ({coverage_pct:.1f}% coverage)")
    elif impact_type == 'TRUNCATED':
        explanations.append(f"truncated ({coverage_pct:.1f}% coverage)")
    elif impact_type == 'SPLIT':
        num_fragments = len(details.get('fragments', []))
        explanations.append(f"split into {num_fragments} fragments")
    elif impact_type == 'MISSING':
        explanations.append("completely missing or unmappable")
    
    # Structural impact explanations
    if 'INVERSION' in structural_impacts:
        explanations.append("sequence inverted")
    if 'DUPLICATION' in structural_impacts:
        copy_number = details.get('copy_number', 1)
        explanations.append(f"duplicated ({copy_number} copies)")
    if 'DELETION' in structural_impacts:
        explanations.append("deleted")
    if 'TRANSLOCATION' in structural_impacts:
        explanations.append("translocated across non-adjacent nodes")
    if 'COPY_NUMBER_CHANGE' in structural_impacts:
        copy_number = details.get('copy_number', 1)
        expected = details.get('expected_copy_number', 1)
        # Only show copy number change if there's actually a difference
        if copy_number != expected:
            explanations.append(f"copy number changed ({copy_number} vs {expected} expected)")
    if 'COMPLEX_REARRANGEMENT' in structural_impacts:
        explanations.append("complex rearrangement detected")
    
    # Orientation changes
    orientation_changes = details.get('orientation_changes', [])
    if orientation_changes:
        explanations.append(f"{len(orientation_changes)} orientation changes")
    
    # Path discontinuities
    discontinuities = details.get('path_discontinuities', [])
    if discontinuities:
        explanations.append(f"{len(discontinuities)} path discontinuities")
    
    # Component spanning
    components = details.get('components_spanned', [])
    if len(components) > 1:
        explanations.append(f"spans {len(components)} graph components")
    
    # Feature-specific explanations
    if feature_type.lower() == 'cds':
        if consequence == 'FRAMESHIFT':
            explanations.append("reading frame disrupted")
        elif consequence == 'START_LOST':
            explanations.append("start codon lost")
        elif consequence == 'STOP_LOST':
            explanations.append("stop codon lost")
        elif consequence == 'INFRAME_INDEL':
            explanations.append("in-frame indel")
    elif feature_type.lower() == 'promoter':
        if consequence == 'CORE_DISRUPTED':
            explanations.append("core promoter disrupted")
        elif consequence == 'TSS_SHIFTED':
            explanations.append("transcription start site shifted")
    
    return "; ".join(explanations)

def print_colored(text: str, color: str = Colors.ENDC, bold: bool = False) -> None:
    """Print text with color formatting."""
    format_code = color
    if bold:
        format_code += Colors.BOLD
    print(f"{format_code}{text}{Colors.ENDC}")

def print_feature_analysis(sample: str, haplotype: str, feature_name: str, 
                         analysis: dict) -> None:
    """
    Print detailed analysis for a single feature.
    
    Args:
        sample: Sample name
        haplotype: Haplotype identifier
        feature_name: Feature name
        analysis: Impact analysis results
    """
    impact_type = analysis['type']
    consequence = analysis['consequence']
    structural_impacts = analysis.get('structural_impacts', [])
    details = analysis['details']
    feature_type = details.get('feature_type', 'unknown')
    
    # Choose color based on impact severity
    impact_color = IMPACT_COLORS.get(impact_type, Colors.ENDC)
    
    # Feature header
    print_colored(f"  📍 {feature_name} ({feature_type})", Colors.BOLD)
    
    # Impact type and consequence
    print(f"     Impact: ", end="")
    print_colored(impact_type, impact_color, bold=True)
    print(f"     Consequence: {consequence}")
    
    # Structural impacts
    if structural_impacts:
        print("     Structural impacts:")
        for struct_impact in structural_impacts:
            color = STRUCTURAL_COLORS.get(struct_impact, Colors.WARNING)
            print(f"       • ", end="")
            print_colored(struct_impact, color, bold=True)
    
    # Explanation
    explanation = get_impact_explanation(feature_name, feature_type, impact_type, 
                                       consequence, structural_impacts, details)
    print(f"     📝 {explanation}")
    
    # Technical details
    coverage = details.get('coverage', 0.0)
    identity = details.get('identity', 0.0)
    score = details.get('score', 0)
    print(f"     📊 Coverage: {coverage:.1%}, Identity: {identity:.1%}, Score: {score}")
    
    print()  # Blank line for spacing

def analyze_and_display_results(gam_data: dict, gfa_file: Path) -> None:
    """
    Run impact analysis and display results with colored output.
    
    Args:
        gam_data: GAM data dictionary
        gfa_file: Path to GFA file
    """
    print_colored("🔬 Initializing ImpactDetector...", Colors.HEADER, bold=True)
    
    # Initialize impact detector
    detector = ImpactDetector(
        gfa_file,
        min_alignment_coverage=0.8,
        min_identity_threshold=0.85  # Lower threshold for demo data
    )
    
    print_colored("🧬 Analyzing feature impacts...", Colors.HEADER, bold=True)
    
    # Run impact analysis
    impact_results = detector.analyze_impacts(gam_data)
    
    # Display results
    print_colored("\n📋 IMPACT ANALYSIS RESULTS", Colors.HEADER, bold=True)
    print_colored("=" * 50, Colors.HEADER)
    
    total_features = 0
    structural_variants_found = 0
    
    for sample in sorted(impact_results.keys()):
        print_colored(f"\n🧬 Sample: {sample}", Colors.OKBLUE, bold=True)
        
        for haplotype in sorted(impact_results[sample].keys()):
            print_colored(f"  🔗 Haplotype: {haplotype}", Colors.OKCYAN, bold=True)
            
            sample_features = impact_results[sample][haplotype]
            if not sample_features:
                print("    (No features found)")
                continue
            
            for feature_name in sorted(sample_features.keys()):
                analysis = sample_features[feature_name]
                total_features += 1
                
                # Check for structural variants
                structural_impacts = analysis.get('structural_impacts', [])
                if structural_impacts:
                    structural_variants_found += 1
                
                print_feature_analysis(sample, haplotype, feature_name, analysis)
    
    # Summary statistics
    summary = detector.get_impact_summary(impact_results)
    
    print_colored("📊 SUMMARY STATISTICS", Colors.HEADER, bold=True)
    print_colored("=" * 30, Colors.HEADER)
    
    print(f"Total features analyzed: {summary['total_features']}")
    print(f"Structural variants detected: {structural_variants_found}")
    
    print("\n🎯 Impact Type Distribution:")
    for impact_type, count in summary['impact_type_counts'].items():
        if count > 0:
            percentage = (count / summary['total_features'] * 100) if summary['total_features'] > 0 else 0
            color = IMPACT_COLORS.get(impact_type, Colors.ENDC)
            print(f"  ", end="")
            print_colored(f"{impact_type}: {count} ({percentage:.1f}%)", color)
    
    print("\n🔄 Structural Impact Distribution:")
    for impact_type, count in summary['structural_impact_counts'].items():
        if count > 0:
            percentage = (count / summary['total_features'] * 100) if summary['total_features'] > 0 else 0
            color = STRUCTURAL_COLORS.get(impact_type, Colors.WARNING)
            print(f"  ", end="")
            print_colored(f"{impact_type}: {count} ({percentage:.1f}%)", color)
    
    print(f"\n📈 Analysis Metrics:")
    print(f"  Average coverage: {summary['avg_coverage']:.1%}")
    print(f"  Average identity: {summary['avg_identity']:.1%}")
    print(f"  Features spanning multiple components: {summary['features_spanning_multiple_components']}")
    print(f"  Orientation changes detected: {summary['orientation_changes_detected']}")
    print(f"  Path discontinuities detected: {summary['path_discontinuities_detected']}")
    
    # Feature type breakdown
    print(f"\n🧮 Feature Type Analysis:")
    for feature_type, data in summary['feature_type_analysis'].items():
        if data['total'] > 0:
            print(f"  {feature_type}: {data['total']} features")
            for impact_type, count in data['impact_types'].items():
                if count > 0:
                    color = IMPACT_COLORS.get(impact_type, Colors.ENDC)
                    print(f"    ", end="")
                    print_colored(f"{impact_type}: {count}", color)

def main():
    """Main demonstration function."""
    print_colored("=" * 60, Colors.HEADER)
    print_colored("🧬 IMPACT DETECTOR DEMONSTRATION", Colors.HEADER, bold=True)
    print_colored("=" * 60, Colors.HEADER)
    
    # Setup logging
    logging.basicConfig(level=logging.INFO)
    
    try:
        # Create test data
        print_colored("\n1. Creating test GFA graph...", Colors.OKBLUE, bold=True)
        gfa_file = create_test_gfa_file()
        print(f"   Created test GFA: {gfa_file}")
        
        print_colored("\n2. Creating test GAM data...", Colors.OKBLUE, bold=True)
        gam_data = create_test_gam_data()
        print(f"   Created GAM data for {len(gam_data)} samples")
        
        # Analyze impacts
        print_colored("\n3. Running impact analysis...", Colors.OKBLUE, bold=True)
        analyze_and_display_results(gam_data, gfa_file)
        
        # Legend
        print_colored("\n🎨 COLOR LEGEND", Colors.HEADER, bold=True)
        print_colored("-" * 20, Colors.HEADER)
        print("Impact Types:")
        for impact_type, color in IMPACT_COLORS.items():
            print(f"  ", end="")
            print_colored(f"■ {impact_type}", color, bold=True)
        
        print("\nStructural Variants:")
        for struct_type, color in STRUCTURAL_COLORS.items():
            print(f"  ", end="")
            print_colored(f"▲ {struct_type}", color, bold=True)
        
    except Exception as e:
        print_colored(f"\n❌ Error during analysis: {e}", Colors.FAIL, bold=True)
        raise
    
    finally:
        # Clean up temporary files
        if 'gfa_file' in locals():
            gfa_file.unlink(missing_ok=True)
            print(f"\n🧹 Cleaned up temporary file: {gfa_file}")
    
    print_colored("\n" + "=" * 60, Colors.HEADER)
    print_colored("✅ IMPACT DETECTOR DEMO COMPLETED!", Colors.OKGREEN, bold=True)
    print_colored("=" * 60, Colors.HEADER)

if __name__ == "__main__":
    main()
