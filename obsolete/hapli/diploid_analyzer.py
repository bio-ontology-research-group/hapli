#!/usr/bin/env python3
"""
Diploid analyzer for comparing impacts across haplotypes.

This module provides functionality to analyze impact data from both haplotypes
of diploid samples and classify zygosity patterns.
"""

import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any, Union
import json

logger = logging.getLogger(__name__)


class DiploidAnalyzer:
    """
    Analyze impact data across haplotypes to determine zygosity patterns.
    
    Compares impact analysis results between haplotypes to classify features
    as homozygous reference, heterozygous, homozygous alternate, or compound
    heterozygous variants.
    """
    
    # Zygosity classifications
    ZYGOSITY_TYPES = {
        'HOMOZYGOUS_REF': 'Both haplotypes are intact/reference-like',
        'HETEROZYGOUS': 'One haplotype affected, other intact',
        'HOMOZYGOUS_ALT': 'Both haplotypes have similar impacts',
        'COMPOUND_HET': 'Both haplotypes affected but with different impacts'
    }
    
    # Loss-of-function impact types and consequences
    LOF_IMPACT_TYPES = {'MISSING', 'SPLIT'}
    LOF_CONSEQUENCES = {
        'MISSING', 'FRAMESHIFT', 'START_LOST', 'STOP_LOST', 
        'NONSENSE', 'LOST', 'CORE_DISRUPTED'
    }
    
    # Structural variants that may cause LOF
    LOF_STRUCTURAL_IMPACTS = {
        'DELETION', 'COMPLEX_REARRANGEMENT'
    }
    
    def __init__(self, lof_coverage_threshold: float = 0.5):
        """
        Initialize diploid analyzer.
        
        Args:
            lof_coverage_threshold: Coverage threshold below which to consider LOF
        """
        self.lof_coverage_threshold = lof_coverage_threshold
    
    def analyze_diploid_impacts(self, impact_results: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]]) -> Dict[str, Dict[str, Dict[str, Any]]]:
        """
        Analyze impact data across haplotypes for diploid samples.
        
        Args:
            impact_results: Impact analysis results from ImpactDetector
            
        Returns:
            Diploid analysis: {sample: {gene: {'zygosity': type, 'hap1_impact': ..., 'hap2_impact': ..., 'details': {...}}}}
        """
        logger.info("Analyzing diploid impacts across haplotypes...")
        
        diploid_results = {}
        
        for sample in impact_results:
            logger.debug(f"Processing sample: {sample}")
            
            # Get haplotype data
            haplotypes = impact_results[sample]
            
            # We expect haplotypes '0' and '1' for diploid analysis
            if '0' not in haplotypes or '1' not in haplotypes:
                logger.warning(f"Sample {sample} missing expected haplotypes (0, 1). "
                             f"Found: {list(haplotypes.keys())}")
                continue
            
            hap1_data = haplotypes['0']
            hap2_data = haplotypes['1']
            
            # Group features by gene
            sample_results = self._analyze_sample_diploid_impacts(sample, hap1_data, hap2_data)
            diploid_results[sample] = sample_results
        
        return diploid_results
    
    def _analyze_sample_diploid_impacts(self, sample: str, 
                                      hap1_data: Dict[str, Dict[str, Any]], 
                                      hap2_data: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """
        Analyze diploid impacts for a single sample.
        
        Args:
            sample: Sample name
            hap1_data: Haplotype 1 impact data
            hap2_data: Haplotype 2 impact data
            
        Returns:
            Gene-level diploid analysis
        """
        # Group features by gene
        hap1_genes = self._group_features_by_gene(hap1_data)
        hap2_genes = self._group_features_by_gene(hap2_data)
        
        # Get all genes present in either haplotype
        all_genes = set(hap1_genes.keys()) | set(hap2_genes.keys())
        
        sample_results = {}
        
        for gene in all_genes:
            gene_analysis = self._analyze_gene_diploid_impact(
                gene, 
                hap1_genes.get(gene, {}), 
                hap2_genes.get(gene, {})
            )
            sample_results[gene] = gene_analysis
        
        return sample_results
    
    def _group_features_by_gene(self, haplotype_data: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Dict[str, Any]]]:
        """
        Group features by their parent gene.
        
        Args:
            haplotype_data: Feature impact data for one haplotype
            
        Returns:
            Dictionary grouped by gene: {gene_name: {feature_name: impact_data}}
        """
        genes = defaultdict(dict)
        
        for feature_name, impact_data in haplotype_data.items():
            gene_name = self._extract_gene_name(feature_name)
            genes[gene_name][feature_name] = impact_data
        
        return dict(genes)
    
    def _extract_gene_name(self, feature_name: str) -> str:
        """
        Extract gene name from feature name.
        
        Handles patterns like:
        - gene_001 -> gene_001
        - cds_001_01 -> gene_001 (inferred)
        - exon_001_01 -> gene_001 (inferred)
        - promoter_002 -> gene_002 (inferred)
        
        Args:
            feature_name: Feature name
            
        Returns:
            Extracted or inferred gene name
        """
        # Direct gene names
        if feature_name.startswith('gene_'):
            return feature_name
        
        # Extract from feature patterns
        patterns = [
            r'^(cds|exon|utr|promoter)_(\d+)(?:_\d+)?',  # cds_001_01 -> 001
            r'^(\w+)_(\d+)',  # any_prefix_number
        ]
        
        for pattern in patterns:
            match = re.match(pattern, feature_name)
            if match:
                number = match.group(2)
                return f"gene_{number}"
        
        # Fallback: use the feature name as gene name
        return feature_name
    
    def _analyze_gene_diploid_impact(self, gene_name: str, 
                                   hap1_features: Dict[str, Dict[str, Any]], 
                                   hap2_features: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """
        Analyze diploid impact for a single gene.
        
        Args:
            gene_name: Gene name
            hap1_features: Haplotype 1 features for this gene
            hap2_features: Haplotype 2 features for this gene
            
        Returns:
            Gene-level diploid analysis
        """
        # Aggregate impact information for each haplotype
        hap1_summary = self._summarize_gene_impact(hap1_features)
        hap2_summary = self._summarize_gene_impact(hap2_features)
        
        # Determine zygosity
        zygosity = self._determine_zygosity(hap1_summary, hap2_summary)
        
        # Check for loss-of-function
        hap1_lof = self._is_loss_of_function(hap1_summary)
        hap2_lof = self._is_loss_of_function(hap2_summary)
        
        # Determine overall gene impact
        gene_impact = self._determine_gene_impact(hap1_summary, hap2_summary, zygosity)
        
        return {
            'zygosity': zygosity,
            'hap1_impact': hap1_summary,
            'hap2_impact': hap2_summary,
            'details': {
                'gene_name': gene_name,
                'hap1_lof': hap1_lof,
                'hap2_lof': hap2_lof,
                'biallelic_lof': hap1_lof and hap2_lof,
                'gene_impact': gene_impact,
                'hap1_feature_count': len(hap1_features),
                'hap2_feature_count': len(hap2_features),
                'compound_het_details': self._analyze_compound_het(hap1_summary, hap2_summary) if zygosity == 'COMPOUND_HET' else None
            }
        }
    
    def _summarize_gene_impact(self, features: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """
        Summarize impact across all features of a gene.
        
        Args:
            features: Feature impact data for the gene
            
        Returns:
            Summarized gene impact
        """
        if not features:
            return {
                'overall_impact': 'MISSING',
                'worst_consequence': 'MISSING',
                'structural_impacts': [],
                'features': {},
                'avg_coverage': 0.0,
                'avg_identity': 0.0,
                'has_intact_features': False,
                'has_lof_features': False
            }
        
        # Collect all impacts and consequences
        impact_types = []
        consequences = []
        structural_impacts = []
        total_coverage = 0.0
        total_identity = 0.0
        
        feature_summaries = {}
        
        for feature_name, impact_data in features.items():
            impact_type = impact_data.get('type', 'MISSING')
            consequence = impact_data.get('consequence', 'MISSING')
            struct_impacts = impact_data.get('structural_impacts', [])
            details = impact_data.get('details', {})
            
            impact_types.append(impact_type)
            consequences.append(consequence)
            structural_impacts.extend(struct_impacts)
            
            total_coverage += details.get('coverage', 0.0)
            total_identity += details.get('identity', 0.0)
            
            feature_summaries[feature_name] = {
                'type': impact_type,
                'consequence': consequence,
                'structural_impacts': struct_impacts,
                'coverage': details.get('coverage', 0.0),
                'identity': details.get('identity', 0.0)
            }
        
        # Determine overall impact (worst case)
        overall_impact = self._determine_worst_impact(impact_types)
        worst_consequence = self._determine_worst_consequence(consequences)
        
        # Calculate averages
        num_features = len(features)
        avg_coverage = total_coverage / num_features if num_features > 0 else 0.0
        avg_identity = total_identity / num_features if num_features > 0 else 0.0
        
        # Check for intact and LOF features
        has_intact = any(impact == 'INTACT' for impact in impact_types)
        has_lof = any(self._is_lof_consequence(cons) for cons in consequences)
        
        return {
            'overall_impact': overall_impact,
            'worst_consequence': worst_consequence,
            'structural_impacts': list(set(structural_impacts)),
            'features': feature_summaries,
            'avg_coverage': avg_coverage,
            'avg_identity': avg_identity,
            'has_intact_features': has_intact,
            'has_lof_features': has_lof
        }
    
    def _determine_worst_impact(self, impact_types: List[str]) -> str:
        """Determine the worst (most severe) impact type."""
        if not impact_types:
            return 'MISSING'
        
        # Severity order (most to least severe)
        severity_order = ['MISSING', 'SPLIT', 'TRUNCATED', 'INTACT']
        
        for impact in severity_order:
            if impact in impact_types:
                return impact
        
        return impact_types[0]  # Fallback
    
    def _determine_worst_consequence(self, consequences: List[str]) -> str:
        """Determine the worst (most severe) consequence."""
        if not consequences:
            return 'MISSING'
        
        # Severity order for consequences
        lof_consequences = ['MISSING', 'FRAMESHIFT', 'START_LOST', 'STOP_LOST', 'NONSENSE', 'LOST']
        moderate_consequences = ['TRUNCATED', 'INFRAME_INDEL', 'PARTIAL', 'CORE_DISRUPTED']
        
        # Check for LOF consequences first
        for cons in lof_consequences:
            if cons in consequences:
                return cons
        
        # Then moderate consequences
        for cons in moderate_consequences:
            if cons in consequences:
                return cons
        
        return consequences[0]  # Fallback
    
    def _determine_zygosity(self, hap1_summary: Dict[str, Any], hap2_summary: Dict[str, Any]) -> str:
        """
        Determine zygosity based on haplotype impact summaries.
        
        Args:
            hap1_summary: Haplotype 1 impact summary
            hap2_summary: Haplotype 2 impact summary
            
        Returns:
            Zygosity classification
        """
        hap1_impact = hap1_summary['overall_impact']
        hap2_impact = hap2_summary['overall_impact']
        
        hap1_intact = hap1_impact == 'INTACT'
        hap2_intact = hap2_impact == 'INTACT'
        
        # Both intact = homozygous reference
        if hap1_intact and hap2_intact:
            return 'HOMOZYGOUS_REF'
        
        # One intact, one affected = heterozygous
        if hap1_intact and not hap2_intact:
            return 'HETEROZYGOUS'
        if hap2_intact and not hap1_intact:
            return 'HETEROZYGOUS'
        
        # Both affected - check if same or different
        if hap1_impact == hap2_impact:
            # Same impact type - check consequences
            hap1_cons = hap1_summary['worst_consequence']
            hap2_cons = hap2_summary['worst_consequence']
            
            if hap1_cons == hap2_cons:
                return 'HOMOZYGOUS_ALT'
            else:
                return 'COMPOUND_HET'
        else:
            # Different impact types
            return 'COMPOUND_HET'
    
    def _is_loss_of_function(self, gene_summary: Dict[str, Any]) -> bool:
        """
        Determine if a gene summary represents loss-of-function.
        
        Args:
            gene_summary: Gene impact summary
            
        Returns:
            True if gene is likely loss-of-function
        """
        overall_impact = gene_summary['overall_impact']
        worst_consequence = gene_summary['worst_consequence']
        structural_impacts = gene_summary['structural_impacts']
        avg_coverage = gene_summary['avg_coverage']
        
        # Check impact type
        if overall_impact in self.LOF_IMPACT_TYPES:
            return True
        
        # Check consequence
        if worst_consequence in self.LOF_CONSEQUENCES:
            return True
        
        # Check structural impacts
        if any(si in self.LOF_STRUCTURAL_IMPACTS for si in structural_impacts):
            return True
        
        # Check low coverage
        if avg_coverage < self.lof_coverage_threshold:
            return True
        
        return False
    
    def _is_lof_consequence(self, consequence: str) -> bool:
        """Check if a consequence is loss-of-function."""
        return consequence in self.LOF_CONSEQUENCES
    
    def _determine_gene_impact(self, hap1_summary: Dict[str, Any], 
                             hap2_summary: Dict[str, Any], 
                             zygosity: str) -> str:
        """
        Determine overall gene impact based on both haplotypes.
        
        Args:
            hap1_summary: Haplotype 1 summary
            hap2_summary: Haplotype 2 summary
            zygosity: Zygosity classification
            
        Returns:
            Overall gene impact
        """
        hap1_lof = self._is_loss_of_function(hap1_summary)
        hap2_lof = self._is_loss_of_function(hap2_summary)
        
        if hap1_lof and hap2_lof:
            return 'BIALLELIC_LOF'
        elif hap1_lof or hap2_lof:
            return 'MONOALLELIC_LOF'
        elif zygosity == 'HOMOZYGOUS_REF':
            return 'NO_IMPACT'
        elif zygosity in ['HETEROZYGOUS', 'HOMOZYGOUS_ALT', 'COMPOUND_HET']:
            return 'VARIANT_IMPACT'
        else:
            return 'UNKNOWN'
    
    def _analyze_compound_het(self, hap1_summary: Dict[str, Any], 
                            hap2_summary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Analyze compound heterozygous impacts in detail.
        
        Args:
            hap1_summary: Haplotype 1 summary
            hap2_summary: Haplotype 2 summary
            
        Returns:
            Detailed compound heterozygous analysis
        """
        return {
            'hap1_impact_type': hap1_summary['overall_impact'],
            'hap2_impact_type': hap2_summary['overall_impact'],
            'hap1_consequence': hap1_summary['worst_consequence'],
            'hap2_consequence': hap2_summary['worst_consequence'],
            'hap1_structural': hap1_summary['structural_impacts'],
            'hap2_structural': hap2_summary['structural_impacts'],
            'impact_difference': self._describe_impact_difference(hap1_summary, hap2_summary)
        }
    
    def _describe_impact_difference(self, hap1_summary: Dict[str, Any], 
                                  hap2_summary: Dict[str, Any]) -> str:
        """
        Describe the difference between compound heterozygous impacts.
        
        Args:
            hap1_summary: Haplotype 1 summary
            hap2_summary: Haplotype 2 summary
            
        Returns:
            Description of the difference
        """
        hap1_impact = hap1_summary['overall_impact']
        hap2_impact = hap2_summary['overall_impact']
        hap1_cons = hap1_summary['worst_consequence']
        hap2_cons = hap2_summary['worst_consequence']
        
        if hap1_impact != hap2_impact:
            return f"Different impact types: {hap1_impact} vs {hap2_impact}"
        elif hap1_cons != hap2_cons:
            return f"Same impact type ({hap1_impact}) but different consequences: {hap1_cons} vs {hap2_cons}"
        else:
            return f"Same impact and consequence but different details"
    
    def get_diploid_summary(self, diploid_results: Dict[str, Dict[str, Dict[str, Any]]]) -> Dict[str, Any]:
        """
        Generate summary statistics for diploid analysis.
        
        Args:
            diploid_results: Results from analyze_diploid_impacts()
            
        Returns:
            Summary statistics
        """
        summary = {
            'total_samples': len(diploid_results),
            'total_genes': 0,
            'zygosity_counts': {zyg_type: 0 for zyg_type in self.ZYGOSITY_TYPES},
            'gene_impact_counts': defaultdict(int),
            'biallelic_lof_genes': 0,
            'monoallelic_lof_genes': 0,
            'compound_het_genes': 0,
            'sample_summaries': {}
        }
        
        for sample, genes in diploid_results.items():
            sample_summary = {
                'total_genes': len(genes),
                'zygosity_counts': {zyg_type: 0 for zyg_type in self.ZYGOSITY_TYPES},
                'biallelic_lof_count': 0,
                'compound_het_count': 0
            }
            
            for gene_name, analysis in genes.items():
                summary['total_genes'] += 1
                
                zygosity = analysis['zygosity']
                gene_impact = analysis['details']['gene_impact']
                
                # Count zygosity types
                summary['zygosity_counts'][zygosity] += 1
                sample_summary['zygosity_counts'][zygosity] += 1
                
                # Count gene impacts
                summary['gene_impact_counts'][gene_impact] += 1
                
                # Count special cases
                if analysis['details']['biallelic_lof']:
                    summary['biallelic_lof_genes'] += 1
                    sample_summary['biallelic_lof_count'] += 1
                
                if analysis['details']['hap1_lof'] or analysis['details']['hap2_lof']:
                    summary['monoallelic_lof_genes'] += 1
                
                if zygosity == 'COMPOUND_HET':
                    summary['compound_het_genes'] += 1
                    sample_summary['compound_het_count'] += 1
            
            summary['sample_summaries'][sample] = sample_summary
        
        # Convert defaultdict to regular dict
        summary['gene_impact_counts'] = dict(summary['gene_impact_counts'])
        
        return summary
    
    def save_results(self, diploid_results: Dict[str, Dict[str, Dict[str, Any]]], 
                    output_file: Path) -> None:
        """
        Save diploid analysis results to JSON file.
        
        Args:
            diploid_results: Results from analyze_diploid_impacts()
            output_file: Output file path
        """
        output_data = {
            'diploid_analysis': diploid_results,
            'summary': self.get_diploid_summary(diploid_results),
            'parameters': {
                'lof_coverage_threshold': self.lof_coverage_threshold
            },
            'zygosity_type_descriptions': self.ZYGOSITY_TYPES
        }
        
        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        logger.info(f"Saved diploid analysis results to {output_file}")


def main():
    """Command-line interface for diploid analyzer."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze diploid impacts across haplotypes")
    parser.add_argument("impact_json", help="Input impact data (JSON from ImpactDetector)")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file for diploid analysis")
    parser.add_argument("--lof-coverage-threshold", type=float, default=0.5,
                       help="Coverage threshold for loss-of-function classification")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Load impact data
    logger.info(f"Loading impact data from {args.impact_json}")
    with open(args.impact_json, 'r') as f:
        impact_data = json.load(f)
    
    # Extract impact analysis from the input (might be nested)
    if 'impact_analysis' in impact_data:
        impact_results = impact_data['impact_analysis']
    else:
        impact_results = impact_data
    
    # Initialize diploid analyzer
    analyzer = DiploidAnalyzer(
        lof_coverage_threshold=args.lof_coverage_threshold
    )
    
    # Analyze diploid impacts
    diploid_results = analyzer.analyze_diploid_impacts(impact_results)
    
    # Save results
    analyzer.save_results(diploid_results, Path(args.output))
    
    # Print summary
    summary = analyzer.get_diploid_summary(diploid_results)
    print(f"\nDiploid Analysis Summary:")
    print(f"Total samples analyzed: {summary['total_samples']}")
    print(f"Total genes analyzed: {summary['total_genes']}")
    
    print(f"\nZygosity distribution:")
    for zyg_type, count in summary['zygosity_counts'].items():
        percentage = (count / summary['total_genes'] * 100) if summary['total_genes'] > 0 else 0
        print(f"  {zyg_type}: {count} ({percentage:.1f}%)")
    
    print(f"\nGene impact distribution:")
    for impact_type, count in summary['gene_impact_counts'].items():
        percentage = (count / summary['total_genes'] * 100) if summary['total_genes'] > 0 else 0
        print(f"  {impact_type}: {count} ({percentage:.1f}%)")
    
    print(f"\nSpecial cases:")
    print(f"  Biallelic LOF genes: {summary['biallelic_lof_genes']}")
    print(f"  Monoallelic LOF genes: {summary['monoallelic_lof_genes']}")
    print(f"  Compound heterozygous genes: {summary['compound_het_genes']}")


if __name__ == "__main__":
    main()
