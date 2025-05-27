#!/usr/bin/env python3
"""
Demo/test script for DiploidAnalyzer functionality.

This script demonstrates how to use the DiploidAnalyzer to analyze diploid
impacts across haplotypes, showing zygosity patterns and loss-of-function
analysis with rich formatting.
"""

import json
import logging
import tempfile
from pathlib import Path
from hapli.impact_detector import ImpactDetector
from hapli.diploid_analyzer import DiploidAnalyzer

try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich.text import Text
    from rich.layout import Layout
    from rich.columns import Columns
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False
    print("Rich not available, using basic formatting")

try:
    from tabulate import tabulate
    TABULATE_AVAILABLE = True
except ImportError:
    TABULATE_AVAILABLE = False
    print("Tabulate not available, using basic formatting")


def create_test_impact_data() -> dict:
    """
    Create comprehensive test impact data for diploid analysis.
    
    Returns:
        Impact data in the format expected by DiploidAnalyzer
    """
    return {
        "sample1": {
            "0": {  # Haplotype 1
                "gene_001": {
                    "type": "INTACT",
                    "consequence": "SYNONYMOUS",
                    "structural_impacts": [],
                    "details": {
                        "feature_type": "gene",
                        "coverage": 0.98,
                        "identity": 0.96,
                        "score": 95
                    }
                },
                "cds_001_01": {
                    "type": "INTACT",
                    "consequence": "SYNONYMOUS",
                    "structural_impacts": [],
                    "details": {
                        "feature_type": "CDS",
                        "coverage": 0.98,
                        "identity": 0.96,
                        "score": 88
                    }
                },
                "gene_002": {
                    "type": "TRUNCATED",
                    "consequence": "PARTIAL",
                    "structural_impacts": ["DELETION"],
                    "details": {
                        "feature_type": "gene",
                        "coverage": 0.65,
                        "identity": 0.92,
                        "score": 70
                    }
                },
                "cds_002_01": {
                    "type": "SPLIT",
                    "consequence": "FRAMESHIFT",
                    "structural_impacts": ["DELETION", "TRANSLOCATION"],
                    "details": {
                        "feature_type": "CDS",
                        "coverage": 0.45,
                        "identity": 0.89,
                        "score": 55
                    }
                },
                "gene_003": {
                    "type": "MISSING",
                    "consequence": "MISSING",
                    "structural_impacts": ["DELETION"],
                    "details": {
                        "feature_type": "gene",
                        "coverage": 0.1,
                        "identity": 0.0,
                        "score": 0
                    }
                }
            },
            "1": {  # Haplotype 2
                "gene_001": {
                    "type": "INTACT",
                    "consequence": "SYNONYMOUS",
                    "structural_impacts": [],
                    "details": {
                        "feature_type": "gene",
                        "coverage": 0.97,
                        "identity": 0.95,
                        "score": 93
                    }
                },
                "cds_001_01": {
                    "type": "INTACT",
                    "consequence": "SYNONYMOUS",
                    "structural_impacts": [],
                    "details": {
                        "feature_type": "CDS",
                        "coverage": 0.97,
                        "identity": 0.95,
                        "score": 86
                    }
                },
                "gene_002": {
                    "type": "INTACT",
                    "consequence": "SYNONYMOUS",
                    "structural_impacts": [],
                    "details": {
                        "feature_type": "gene",
                        "coverage": 0.99,
                        "identity": 0.98,
                        "score": 98
                    }
                },
                "cds_002_01": {
                    "type": "INTACT",
                    "consequence": "SYNONYMOUS",
                    "structural_impacts": [],
                    "details": {
                        "feature_type": "CDS",
                        "coverage": 0.99,
                        "identity": 0.98,
                        "score": 95
                    }
                },
                "gene_003": {
                    "type": "MISSING",
                    "consequence": "MISSING",
                    "structural_impacts": ["DELETION"],
                    "details": {
                        "feature_type": "gene",
                        "coverage": 0.05,
                        "identity": 0.0,
                        "score": 0
                    }
                }
            }
        },
        "sample2": {
            "0": {  # Haplotype 1
                "gene_001": {
                    "type": "TRUNCATED",
                    "consequence": "INFRAME_INDEL",
                    "structural_impacts": ["INVERSION"],
                    "details": {
                        "feature_type": "gene",
                        "coverage": 0.75,
                        "identity": 0.88,
                        "score": 75
                    }
                },
                "cds_001_01": {
                    "type": "SPLIT",
                    "consequence": "FRAMESHIFT",
                    "structural_impacts": ["INVERSION"],
                    "details": {
                        "feature_type": "CDS",
                        "coverage": 0.60,
                        "identity": 0.85,
                        "score": 65
                    }
                },
                "gene_004": {
                    "type": "SPLIT",
                    "consequence": "START_LOST",
                    "structural_impacts": ["COMPLEX_REARRANGEMENT"],
                    "details": {
                        "feature_type": "gene",
                        "coverage": 0.40,
                        "identity": 0.80,
                        "score": 50
                    }
                },
                "cds_004_01": {
                    "type": "MISSING",
                    "consequence": "START_LOST",
                    "structural_impacts": ["DELETION"],
                    "details": {
                        "feature_type": "CDS",
                        "coverage": 0.15,
                        "identity": 0.70,
                        "score": 25
                    }
                }
            },
            "1": {  # Haplotype 2
                "gene_001": {
                    "type": "SPLIT",
                    "consequence": "STOP_LOST",
                    "structural_impacts": ["TRANSLOCATION"],
                    "details": {
                        "feature_type": "gene",
                        "coverage": 0.55,
                        "identity": 0.82,
                        "score": 60
                    }
                },
                "cds_001_01": {
                    "type": "TRUNCATED",
                    "consequence": "STOP_LOST",
                    "structural_impacts": ["TRANSLOCATION"],
                    "details": {
                        "feature_type": "CDS",
                        "coverage": 0.70,
                        "identity": 0.85,
                        "score": 70
                    }
                },
                "gene_004": {
                    "type": "SPLIT",
                    "consequence": "FRAMESHIFT",
                    "structural_impacts": ["DUPLICATION"],
                    "details": {
                        "feature_type": "gene",
                        "coverage": 0.35,
                        "identity": 0.75,
                        "score": 45
                    }
                },
                "cds_004_01": {
                    "type": "SPLIT",
                    "consequence": "FRAMESHIFT",
                    "structural_impacts": ["DUPLICATION"],
                    "details": {
                        "feature_type": "CDS",
                        "coverage": 0.30,
                        "identity": 0.75,
                        "score": 40
                    }
                }
            }
        }
    }


def get_zygosity_color(zygosity: str) -> str:
    """Get color code for zygosity type."""
    if not RICH_AVAILABLE:
        return ""
    
    colors = {
        'HOMOZYGOUS_REF': 'green',
        'HETEROZYGOUS': 'yellow',
        'HOMOZYGOUS_ALT': 'red',
        'COMPOUND_HET': 'magenta'
    }
    return colors.get(zygosity, 'white')


def get_lof_symbol(is_lof: bool) -> str:
    """Get symbol for loss-of-function status."""
    return "⚠️" if is_lof else "✅"


def print_rich_summary(diploid_results: dict, summary: dict) -> None:
    """Print summary using Rich formatting."""
    if not RICH_AVAILABLE:
        print_basic_summary(diploid_results, summary)
        return
    
    console = Console()
    
    # Title
    console.print("\n🧬 DIPLOID ANALYSIS RESULTS", style="bold magenta", justify="center")
    console.print("=" * 60, style="magenta")
    
    # Overall summary
    summary_table = Table(title="📊 Overall Summary", show_header=True, header_style="bold blue")
    summary_table.add_column("Metric", style="cyan")
    summary_table.add_column("Count", justify="right", style="white")
    summary_table.add_column("Percentage", justify="right", style="yellow")
    
    total_genes = summary['total_genes']
    summary_table.add_row("Total Samples", str(summary['total_samples']), "-")
    summary_table.add_row("Total Genes", str(total_genes), "-")
    summary_table.add_row("Biallelic LOF", str(summary['biallelic_lof_genes']), 
                         f"{(summary['biallelic_lof_genes']/total_genes*100):.1f}%" if total_genes > 0 else "0%")
    summary_table.add_row("Compound Het", str(summary['compound_het_genes']), 
                         f"{(summary['compound_het_genes']/total_genes*100):.1f}%" if total_genes > 0 else "0%")
    
    console.print(summary_table)
    
    # Zygosity distribution
    zyg_table = Table(title="🎯 Zygosity Distribution", show_header=True, header_style="bold blue")
    zyg_table.add_column("Zygosity", style="cyan")
    zyg_table.add_column("Count", justify="right", style="white")
    zyg_table.add_column("Percentage", justify="right", style="yellow")
    
    for zyg_type, count in summary['zygosity_counts'].items():
        percentage = f"{(count/total_genes*100):.1f}%" if total_genes > 0 else "0%"
        color = get_zygosity_color(zyg_type)
        zyg_table.add_row(zyg_type, str(count), percentage, style=color)
    
    console.print(zyg_table)


def print_rich_detailed_results(diploid_results: dict) -> None:
    """Print detailed results using Rich formatting."""
    if not RICH_AVAILABLE:
        print_basic_detailed_results(diploid_results)
        return
    
    console = Console()
    
    for sample_name, genes in diploid_results.items():
        console.print(f"\n🧬 Sample: [bold cyan]{sample_name}[/bold cyan]")
        
        # Create table for this sample
        gene_table = Table(show_header=True, header_style="bold blue")
        gene_table.add_column("Gene", style="cyan")
        gene_table.add_column("Zygosity", style="white")
        gene_table.add_column("Hap1 Impact", style="white")
        gene_table.add_column("Hap2 Impact", style="white")
        gene_table.add_column("Hap1 LOF", justify="center", style="white")
        gene_table.add_column("Hap2 LOF", justify="center", style="white")
        gene_table.add_column("Gene Impact", style="white")
        
        # Sort genes for consistent display
        for gene_name in sorted(genes.keys()):
            analysis = genes[gene_name]
            zygosity = analysis['zygosity']
            details = analysis['details']
            
            # Get impact summaries
            hap1_impact = analysis['hap1_impact']['overall_impact']
            hap2_impact = analysis['hap2_impact']['overall_impact']
            
            # LOF status
            hap1_lof = get_lof_symbol(details['hap1_lof'])
            hap2_lof = get_lof_symbol(details['hap2_lof'])
            
            # Color coding
            zyg_color = get_zygosity_color(zygosity)
            
            gene_table.add_row(
                gene_name,
                f"[{zyg_color}]{zygosity}[/{zyg_color}]",
                hap1_impact,
                hap2_impact,
                hap1_lof,
                hap2_lof,
                details['gene_impact']
            )
        
        console.print(gene_table)
        
        # Highlight compound heterozygous cases
        compound_hets = [gene for gene, analysis in genes.items() 
                        if analysis['zygosity'] == 'COMPOUND_HET']
        
        if compound_hets:
            console.print(f"\n🔍 [bold magenta]Compound Heterozygous Cases in {sample_name}:[/bold magenta]")
            for gene in compound_hets:
                analysis = genes[gene]
                details = analysis['details']['compound_het_details']
                if details:
                    console.print(f"  • [cyan]{gene}[/cyan]: {details['impact_difference']}")


def create_rich_summary_matrix(diploid_results: dict) -> None:
    """Create a gene × sample matrix showing zygosity patterns."""
    if not RICH_AVAILABLE:
        create_basic_summary_matrix(diploid_results)
        return
    
    console = Console()
    
    # Collect all genes and samples
    all_genes = set()
    all_samples = list(diploid_results.keys())
    
    for sample_genes in diploid_results.values():
        all_genes.update(sample_genes.keys())
    
    all_genes = sorted(list(all_genes))
    
    # Create matrix table
    matrix_table = Table(title="📊 Gene × Sample Zygosity Matrix", show_header=True, header_style="bold blue")
    matrix_table.add_column("Gene", style="cyan")
    
    for sample in all_samples:
        matrix_table.add_column(sample, justify="center", style="white")
    
    # Add data rows
    for gene in all_genes:
        row = [gene]
        for sample in all_samples:
            if gene in diploid_results[sample]:
                zygosity = diploid_results[sample][gene]['zygosity']
                # Abbreviate zygosity for matrix display
                abbrev = {
                    'HOMOZYGOUS_REF': 'REF',
                    'HETEROZYGOUS': 'HET',
                    'HOMOZYGOUS_ALT': 'ALT',
                    'COMPOUND_HET': 'C.HET'
                }
                color = get_zygosity_color(zygosity)
                row.append(f"[{color}]{abbrev.get(zygosity, zygosity)}[/{color}]")
            else:
                row.append("-")
        
        matrix_table.add_row(*row)
    
    console.print(f"\n{matrix_table}")


def print_basic_summary(diploid_results: dict, summary: dict) -> None:
    """Print summary using basic formatting."""
    print("\n" + "=" * 60)
    print("DIPLOID ANALYSIS RESULTS")
    print("=" * 60)
    
    total_genes = summary['total_genes']
    print(f"\nOverall Summary:")
    print(f"  Total Samples: {summary['total_samples']}")
    print(f"  Total Genes: {total_genes}")
    print(f"  Biallelic LOF: {summary['biallelic_lof_genes']} ({(summary['biallelic_lof_genes']/total_genes*100):.1f}%)" if total_genes > 0 else "  Biallelic LOF: 0 (0%)")
    print(f"  Compound Het: {summary['compound_het_genes']} ({(summary['compound_het_genes']/total_genes*100):.1f}%)" if total_genes > 0 else "  Compound Het: 0 (0%)")
    
    print(f"\nZygosity Distribution:")
    for zyg_type, count in summary['zygosity_counts'].items():
        percentage = f"{(count/total_genes*100):.1f}%" if total_genes > 0 else "0%"
        print(f"  {zyg_type}: {count} ({percentage})")


def print_basic_detailed_results(diploid_results: dict) -> None:
    """Print detailed results using basic formatting."""
    for sample_name, genes in diploid_results.items():
        print(f"\nSample: {sample_name}")
        print("-" * 40)
        
        if TABULATE_AVAILABLE:
            headers = ["Gene", "Zygosity", "Hap1 Impact", "Hap2 Impact", "Hap1 LOF", "Hap2 LOF", "Gene Impact"]
            rows = []
            
            for gene_name in sorted(genes.keys()):
                analysis = genes[gene_name]
                details = analysis['details']
                
                rows.append([
                    gene_name,
                    analysis['zygosity'],
                    analysis['hap1_impact']['overall_impact'],
                    analysis['hap2_impact']['overall_impact'],
                    "YES" if details['hap1_lof'] else "NO",
                    "YES" if details['hap2_lof'] else "NO",
                    details['gene_impact']
                ])
            
            print(tabulate(rows, headers=headers, tablefmt="grid"))
        else:
            for gene_name in sorted(genes.keys()):
                analysis = genes[gene_name]
                details = analysis['details']
                print(f"  {gene_name}:")
                print(f"    Zygosity: {analysis['zygosity']}")
                print(f"    Hap1 Impact: {analysis['hap1_impact']['overall_impact']}")
                print(f"    Hap2 Impact: {analysis['hap2_impact']['overall_impact']}")
                print(f"    LOF Status: Hap1={'YES' if details['hap1_lof'] else 'NO'}, Hap2={'YES' if details['hap2_lof'] else 'NO'}")
                print(f"    Gene Impact: {details['gene_impact']}")
        
        # Highlight compound heterozygous cases
        compound_hets = [gene for gene, analysis in genes.items() 
                        if analysis['zygosity'] == 'COMPOUND_HET']
        
        if compound_hets:
            print(f"\nCompound Heterozygous Cases in {sample_name}:")
            for gene in compound_hets:
                analysis = genes[gene]
                details = analysis['details']['compound_het_details']
                if details:
                    print(f"  • {gene}: {details['impact_difference']}")


def create_basic_summary_matrix(diploid_results: dict) -> None:
    """Create a gene × sample matrix using basic formatting."""
    # Collect all genes and samples
    all_genes = set()
    all_samples = list(diploid_results.keys())
    
    for sample_genes in diploid_results.values():
        all_genes.update(sample_genes.keys())
    
    all_genes = sorted(list(all_genes))
    
    if TABULATE_AVAILABLE:
        headers = ["Gene"] + all_samples
        rows = []
        
        for gene in all_genes:
            row = [gene]
            for sample in all_samples:
                if gene in diploid_results[sample]:
                    zygosity = diploid_results[sample][gene]['zygosity']
                    # Abbreviate zygosity for matrix display
                    abbrev = {
                        'HOMOZYGOUS_REF': 'REF',
                        'HETEROZYGOUS': 'HET',
                        'HOMOZYGOUS_ALT': 'ALT',
                        'COMPOUND_HET': 'C.HET'
                    }
                    row.append(abbrev.get(zygosity, zygosity))
                else:
                    row.append("-")
            rows.append(row)
        
        print(f"\nGene × Sample Zygosity Matrix:")
        print(tabulate(rows, headers=headers, tablefmt="grid"))
    else:
        print(f"\nGene × Sample Zygosity Matrix:")
        print("Gene\t" + "\t".join(all_samples))
        for gene in all_genes:
            row = [gene]
            for sample in all_samples:
                if gene in diploid_results[sample]:
                    zygosity = diploid_results[sample][gene]['zygosity']
                    abbrev = {
                        'HOMOZYGOUS_REF': 'REF',
                        'HETEROZYGOUS': 'HET',
                        'HOMOZYGOUS_ALT': 'ALT',
                        'COMPOUND_HET': 'C.HET'
                    }
                    row.append(abbrev.get(zygosity, zygosity))
                else:
                    row.append("-")
            print("\t".join(row))


def main():
    """Main demonstration function."""
    if RICH_AVAILABLE:
        console = Console()
        console.print("🧬 DIPLOID ANALYZER DEMONSTRATION", style="bold magenta", justify="center")
        console.print("=" * 60, style="magenta")
    else:
        print("=" * 60)
        print("DIPLOID ANALYZER DEMONSTRATION")
        print("=" * 60)
    
    # Setup logging
    logging.basicConfig(level=logging.INFO)
    
    try:
        # Create test impact data
        if RICH_AVAILABLE:
            console = Console()
            console.print("\n1. Creating test impact data...", style="bold blue")
        else:
            print("\n1. Creating test impact data...")
        
        impact_data = create_test_impact_data()
        print(f"   Created impact data for {len(impact_data)} samples")
        
        # Initialize diploid analyzer
        if RICH_AVAILABLE:
            console.print("\n2. Initializing DiploidAnalyzer...", style="bold blue")
        else:
            print("\n2. Initializing DiploidAnalyzer...")
        
        analyzer = DiploidAnalyzer(lof_coverage_threshold=0.5)
        
        # Analyze diploid impacts
        if RICH_AVAILABLE:
            console.print("\n3. Analyzing diploid impacts...", style="bold blue")
        else:
            print("\n3. Analyzing diploid impacts...")
        
        diploid_results = analyzer.analyze_diploid_impacts(impact_data)
        
        # Get summary
        summary = analyzer.get_diploid_summary(diploid_results)
        
        # Display results
        if RICH_AVAILABLE:
            console.print("\n4. Displaying results...", style="bold blue")
            print_rich_summary(diploid_results, summary)
            print_rich_detailed_results(diploid_results)
            create_rich_summary_matrix(diploid_results)
        else:
            print("\n4. Displaying results...")
            print_basic_summary(diploid_results, summary)
            print_basic_detailed_results(diploid_results)
            create_basic_summary_matrix(diploid_results)
        
        # Key findings
        if RICH_AVAILABLE:
            console.print("\n🔍 KEY FINDINGS", style="bold green")
        else:
            print("\nKEY FINDINGS")
            print("-" * 20)
        
        key_findings = []
        
        # Find biallelic LOF genes
        biallelic_lof = []
        compound_hets = []
        
        for sample, genes in diploid_results.items():
            for gene, analysis in genes.items():
                if analysis['details']['biallelic_lof']:
                    biallelic_lof.append(f"{sample}:{gene}")
                if analysis['zygosity'] == 'COMPOUND_HET':
                    compound_hets.append(f"{sample}:{gene}")
        
        if biallelic_lof:
            finding = f"⚠️  Biallelic LOF detected in: {', '.join(biallelic_lof)}"
            if RICH_AVAILABLE:
                console.print(finding, style="bold red")
            else:
                print(finding)
        
        if compound_hets:
            finding = f"🔀 Compound heterozygous variants in: {', '.join(compound_hets)}"
            if RICH_AVAILABLE:
                console.print(finding, style="bold magenta")
            else:
                print(finding)
        
        # Sample-specific findings
        for sample in diploid_results:
            sample_summary = summary['sample_summaries'][sample]
            if sample_summary['biallelic_lof_count'] > 0:
                finding = f"📊 {sample}: {sample_summary['biallelic_lof_count']} genes with biallelic LOF"
                if RICH_AVAILABLE:
                    console.print(finding, style="yellow")
                else:
                    print(finding)
        
    except Exception as e:
        if RICH_AVAILABLE:
            console = Console()
            console.print(f"\n❌ Error during analysis: {e}", style="bold red")
        else:
            print(f"\nError during analysis: {e}")
        raise
    
    if RICH_AVAILABLE:
        console = Console()
        console.print("\n" + "=" * 60, style="magenta")
        console.print("✅ DIPLOID ANALYZER DEMO COMPLETED!", style="bold green", justify="center")
        console.print("=" * 60, style="magenta")
    else:
        print("\n" + "=" * 60)
        print("DIPLOID ANALYZER DEMO COMPLETED!")
        print("=" * 60)


if __name__ == "__main__":
    main()
