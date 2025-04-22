"""
Unit tests for the reporting module.
"""

import os
import unittest
import tempfile
import json
import shutil
from datetime import datetime
from typing import Dict, List, Optional, Tuple

from Bio.SeqFeature import SeqFeature, FeatureLocation
from rdflib import Graph

from src.analysis.impact_classifier import ImpactType
from src.analysis.summary_generator import AnalysisSummary, FeatureSummary
from src.analysis.variant_detector import Variant, VariantType
from src.reporting.report_generator import ReportGenerator
from src.reporting.comparative_report import ComparativeReportGenerator

class TestReporting(unittest.TestCase):
    """Tests for the reporting module."""
    
    def setUp(self):
        """Set up test data."""
        self.test_dir = tempfile.mkdtemp()
        
        # Create data directories
        self.output_dir = os.path.join(self.test_dir, "output")
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Create test data
        self.single_summary = self._create_test_summary("path1", "ref1")
        self.multiple_summaries = {
            "path1": self._create_test_summary("path1", "ref1"),
            "path2": self._create_test_summary("path2", "ref1", variant=True)
        }
        
        # Create report generators
        self.report_generator = ReportGenerator()
        self.comparative_generator = ComparativeReportGenerator()
        
        # Paths to ShEx schema
        self.shex_schema = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 
                                      "src", "reporting", "schemas", "shex_schema.shex")
    
    def tearDown(self):
        """Clean up temporary files."""
        shutil.rmtree(self.test_dir)
    
    def _create_test_summary(self, path_id: str, ref_id: str, variant: bool = False) -> AnalysisSummary:
        """Create a test analysis summary."""
        # Create features
        feature_summaries = {}
        
        # Gene feature
        gene_feature = self._create_mock_feature(100, 500, 1, "gene")
        gene_summary = FeatureSummary(
            feature_id="gene1",
            feature_type="gene",
            impact_type=ImpactType.PRESENT,
            variants=[
                Variant(VariantType.SNP, 150, "A", "G", 1, 40.0),
                Variant(VariantType.SNP, 250, "C", "T", 1, 40.0),
                Variant(VariantType.INSERTION, 350, "", "ATG", 3, 35.0)
            ] if variant else None,
            parent_features=[],
            sequence_identity=0.95,
            coverage=0.98,
            path_id=path_id,
            location=(100, 500, 1)
        )
        feature_summaries["gene1"] = gene_summary
        
        # Exon feature
        exon_feature = self._create_mock_feature(200, 300, 1, "exon")
        exon_summary = FeatureSummary(
            feature_id="exon1",
            feature_type="exon",
            impact_type=ImpactType.PRESENT,
            variants=[
                Variant(VariantType.SNP, 250, "C", "T", 1, 40.0)
            ] if variant else None,
            parent_features=["gene1"],
            sequence_identity=0.98,
            coverage=1.0,
            path_id=path_id,
            location=(200, 300, 1)
        )
        feature_summaries["exon1"] = exon_summary
        
        # Modified feature
        mod_feature = self._create_mock_feature(600, 800, 1, "gene")
        mod_summary = FeatureSummary(
            feature_id="gene2",
            feature_type="gene",
            impact_type=ImpactType.MODIFIED,
            variants=[
                Variant(VariantType.SNP, 650, "A", "G", 1, 40.0),
                Variant(VariantType.DELETION, 700, "ACGT", "", 4, 35.0)
            ] if variant else None,
            parent_features=[],
            sequence_identity=0.85,
            coverage=0.95,
            path_id=path_id,
            location=(600, 800, 1)
        )
        feature_summaries["gene2"] = mod_summary
        
        # Absent feature
        abs_feature = self._create_mock_feature(900, 1000, 1, "gene")
        abs_summary = FeatureSummary(
            feature_id="gene3",
            feature_type="gene",
            impact_type=ImpactType.ABSENT,
            parent_features=[],
            path_id=path_id,
            location=(900, 1000, 1)
        )
        feature_summaries["gene3"] = abs_summary
        
        # Create impact counts
        impact_counts = {
            "present": 2,
            "modified": 1,
            "absent": 1
        }
        
        # Create variant counts
        variant_counts = {
            "SNP": 3,
            "INSERTION": 1,
            "DELETION": 1
        } if variant else {}
        
        # Create analysis summary
        return AnalysisSummary(
            path_id=path_id,
            feature_count=len(feature_summaries),
            feature_by_impact=impact_counts,
            variant_counts=variant_counts,
            reconciliation_counts={},
            feature_summaries=feature_summaries
        )
    
    def _create_mock_feature(self, start: int, end: int, strand: int = 1, feature_type: str = "gene") -> SeqFeature:
        """Create a mock SeqFeature object."""
        location = FeatureLocation(start, end, strand)
        qualifiers = {"ID": [f"{feature_type}1"], "Name": [f"Test {feature_type}"]}
        return SeqFeature(location, type=feature_type, qualifiers=qualifiers)
    
    def test_single_haplotype_tsv_report(self):
        """Test generating a TSV report for a single haplotype."""
        output_file = os.path.join(self.output_dir, "single_report.tsv")
        
        # Generate the report
        result_file = self.report_generator.generate_report(
            self.single_summary, output_file, format_type='tsv')
            
        # Check that the file exists
        self.assertTrue(os.path.exists(result_file))
        
        # Check file content
        with open(result_file, 'r') as f:
            content = f.read()
            
        # Verify key content is present
        self.assertIn("feature_id", content)  # Header
        self.assertIn("gene1", content)  # Gene feature
        self.assertIn("exon1", content)  # Exon feature
        self.assertIn("PRESENT", content)  # Impact type
        self.assertIn("MODIFIED", content)  # Another impact type
        self.assertIn("ABSENT", content)  # Another impact type
    
    def test_single_haplotype_json_report(self):
        """Test generating a JSON report for a single haplotype."""
        output_file = os.path.join(self.output_dir, "single_report.json")
        
        # Generate the report
        result_file = self.report_generator.generate_report(
            self.single_summary, output_file, format_type='json')
            
        # Check that the file exists
        self.assertTrue(os.path.exists(result_file))
        
        # Load and validate JSON
        with open(result_file, 'r') as f:
            json_data = json.load(f)
            
        # Check structure
        self.assertIn('metadata', json_data)
        self.assertIn('features', json_data)
        self.assertEqual(json_data['metadata']['path_id'], "path1")
        self.assertEqual(len(json_data['features']), 4)
        self.assertIn('gene1', json_data['features'])
        self.assertEqual(json_data['features']['gene1']['impact_type'], 'present')
    
    def test_single_haplotype_rdf_report(self):
        """Test generating an RDF report for a single haplotype."""
        output_file = os.path.join(self.output_dir, "single_report.ttl")
        
        # Generate the report
        result_file = self.report_generator.generate_report(
            self.single_summary, output_file, format_type='rdf', rdf_format='turtle')
            
        # Check that the file exists
        self.assertTrue(os.path.exists(result_file))
        
        # Validate RDF structure
        g = Graph()
        g.parse(result_file, format='turtle')
        
        # Check basic structure
        self.assertGreater(len(g), 0, "Graph should not be empty")
        
        # Check the ShEx validation
        if os.path.exists(self.shex_schema):
            is_valid = self.report_generator.validate_rdf_report(result_file, self.shex_schema)
            self.assertTrue(is_valid, "RDF report should validate against the ShEx schema")
    
    def test_comparative_report_generation(self):
        """Test generating a comparative report."""
        output_file = os.path.join(self.output_dir, "comparative_report.json")
        
        # Generate the report
        result_file = self.comparative_generator.generate_comparative_report(
            self.multiple_summaries, output_file, format_type='json')
            
        # Check that the file exists
        self.assertTrue(os.path.exists(result_file))
        
        # Load and validate JSON
        with open(result_file, 'r') as f:
            json_data = json.load(f)
            
        # Check structure
        self.assertIn('metadata', json_data)
        self.assertIn('features', json_data)
        self.assertIn('paths', json_data['metadata'])
        self.assertIn('path1', json_data['metadata']['paths'])
        self.assertIn('path2', json_data['metadata']['paths'])
        self.assertEqual(len(json_data['features']), 4)
    
    def test_impact_type_filtering(self):
        """Test filtering reports by impact type."""
        output_file = os.path.join(self.output_dir, "filtered_report.json")
        
        # Generate the report with filtering
        result_file = self.report_generator.generate_report(
            self.single_summary, output_file, format_type='json',
            impact_types=[ImpactType.PRESENT])
            
        # Check that the file exists
        self.assertTrue(os.path.exists(result_file))
        
        # Load and validate JSON
        with open(result_file, 'r') as f:
            json_data = json.load(f)
            
        # Check filtering worked
        self.assertEqual(len(json_data['features']), 2)  # Only PRESENT features
        self.assertIn('gene1', json_data['features'])
        self.assertIn('exon1', json_data['features'])
        self.assertNotIn('gene2', json_data['features'])  # MODIFIED, should be filtered out
        self.assertNotIn('gene3', json_data['features'])  # ABSENT, should be filtered out
    
    def test_region_filtering(self):
        """Test filtering reports by genomic region."""
        output_file = os.path.join(self.output_dir, "region_filtered_report.json")
        
        # Generate the report with region filtering
        result_file = self.report_generator.generate_report(
            self.single_summary, output_file, format_type='json',
            region={'start': 50, 'end': 550})
            
        # Check that the file exists
        self.assertTrue(os.path.exists(result_file))
        
        # Load and validate JSON
        with open(result_file, 'r') as f:
            json_data = json.load(f)
            
        # Check filtering worked
        self.assertEqual(len(json_data['features']), 2)  # Only features in region
        self.assertIn('gene1', json_data['features'])
        self.assertIn('exon1', json_data['features'])
        self.assertNotIn('gene2', json_data['features'])  # Outside region (600-800)
        self.assertNotIn('gene3', json_data['features'])  # Outside region (900-1000)
    
    def test_rdf_different_formats(self):
        """Test generating RDF reports in different formats."""
        # Test formats: turtle, xml, json-ld, ntriples
        formats = [
            ('turtle', '.ttl'),
            ('xml', '.rdf'),
            ('json-ld', '.jsonld'),
            ('ntriples', '.nt')
        ]
        
        for rdf_format, extension in formats:
            output_file = os.path.join(self.output_dir, f"report{extension}")
            
            # Generate the report
            result_file = self.report_generator.generate_report(
                self.single_summary, output_file, format_type='rdf', rdf_format=rdf_format)
                
            # Check that the file exists
            self.assertTrue(os.path.exists(result_file), f"Failed to create {rdf_format} format")
            
            # Validate it's a parseable RDF file
            g = Graph()
            g.parse(result_file, format=rdf_format)
            self.assertGreater(len(g), 0, f"Graph should not be empty for {rdf_format}")
    
    def test_consensus_feature_identification(self):
        """Test identifying consensus features across multiple paths."""
        # Find consensus features
        consensus = self.comparative_generator.identify_consensus_features(
            self.multiple_summaries, threshold=0.9)
            
        # Both summaries have the same features, so all should be consensus except ABSENT ones
        self.assertIn('gene1', consensus)
        self.assertIn('exon1', consensus)
        self.assertIn('gene2', consensus)
        self.assertNotIn('gene3', consensus)  # ABSENT in all paths
    
    def test_discriminating_feature_identification(self):
        """Test identifying discriminating features across multiple paths."""
        # Create a modified set of summaries with a discriminating feature
        modified_summaries = self.multiple_summaries.copy()
        
        # Add a feature that only exists in path2
        path2_summary = modified_summaries['path2']
        
        # Create a new feature with discriminating characteristics
        new_feature_summary = FeatureSummary(
            feature_id="gene4",
            feature_type="gene",
            impact_type=ImpactType.PRESENT,
            variants=None,
            parent_features=[],
            sequence_identity=0.98,
            coverage=1.0,
            path_id="path2",
            location=(1100, 1200, 1)
        )
        
        # Add it to the path2 summary's feature_summaries
        path2_summary.feature_summaries['gene4'] = new_feature_summary
        path2_summary.feature_count += 1
        path2_summary.feature_by_impact['present'] = path2_summary.feature_by_impact.get('present', 0) + 1
        
        # Find discriminating features
        discriminating = self.comparative_generator.identify_discriminating_features(
            modified_summaries)
            
        # Check that gene4 is identified as discriminating
        self.assertIn('gene4', discriminating)
        self.assertEqual(discriminating['gene4'], ['path2'])

if __name__ == '__main__':
    unittest.main()
