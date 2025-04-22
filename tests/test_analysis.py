"""
Unit tests for the annotation analysis module.
"""

import os
import unittest
from unittest.mock import MagicMock, patch

import tempfile
from io import StringIO
from typing import Dict, List

import networkx as nx
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from src.analysis.impact_classifier import ImpactClassifier, ImpactType
from src.analysis.variant_detector import VariantDetector, Variant, VariantType
from src.analysis.reconciliation import FeatureReconciler, ReconciliationStrategy, ReconciliationResult
from src.analysis.summary_generator import SummaryGenerator, FeatureSummary, AnalysisSummary

# Test data directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")

class TestImpactClassifier(unittest.TestCase):
    """Tests for the ImpactClassifier module."""
    
    def setUp(self):
        self.classifier = ImpactClassifier()
    
    def _create_mock_feature(self, start, end, strand=1, qualifiers=None):
        """Create a mock SeqFeature for testing."""
        feature = MagicMock(spec=SeqFeature)
        feature.location = FeatureLocation(start, end, strand=strand)
        feature.qualifiers = qualifiers or {}
        return feature
    
    def test_classify_present_feature(self):
        """Test classification of a present feature with high identity and coverage."""
        ref_feature = self._create_mock_feature(100, 200)
        aligned_feature = self._create_mock_feature(300, 400, qualifiers={
            "coverage": ["0.95"],
            "identity": ["0.98"]
        })
        
        impact_type, metadata = self.classifier.classify_feature_impact(
            ref_feature, aligned_feature
        )
        
        self.assertEqual(impact_type, ImpactType.PRESENT)
        self.assertAlmostEqual(metadata["coverage"], 0.95)
        self.assertAlmostEqual(metadata["identity"], 0.98)
    
    def test_classify_absent_feature(self):
        """Test classification of an absent feature."""
        ref_feature = self._create_mock_feature(100, 200)
        
        impact_type, metadata = self.classifier.classify_feature_impact(
            ref_feature, None
        )
        
        self.assertEqual(impact_type, ImpactType.ABSENT)
        self.assertEqual(metadata["reason"], "No alignment found")
    
    def test_classify_modified_feature(self):
        """Test classification of a modified feature with medium identity."""
        ref_feature = self._create_mock_feature(100, 200)
        aligned_feature = self._create_mock_feature(300, 400, qualifiers={
            "coverage": ["0.9"],
            "identity": ["0.7"]  # Below threshold for "present"
        })
        
        impact_type, metadata = self.classifier.classify_feature_impact(
            ref_feature, aligned_feature
        )
        
        self.assertEqual(impact_type, ImpactType.MODIFIED)
        self.assertAlmostEqual(metadata["coverage"], 0.9)
        self.assertAlmostEqual(metadata["identity"], 0.7)
    
    def test_classify_truncated_feature(self):
        """Test classification of a truncated feature."""
        ref_feature = self._create_mock_feature(100, 200)  # Length 100
        aligned_feature = self._create_mock_feature(300, 350, qualifiers={
            "coverage": ["0.5"],
            "identity": ["0.95"]
        })  # Length 50
        
        impact_type, metadata = self.classifier.classify_feature_impact(
            ref_feature, aligned_feature
        )
        
        self.assertEqual(impact_type, ImpactType.TRUNCATED)
        self.assertAlmostEqual(metadata["length_ratio"], 0.5)
    
    def test_classify_expanded_feature(self):
        """Test classification of an expanded feature."""
        ref_feature = self._create_mock_feature(100, 200)  # Length 100
        aligned_feature = self._create_mock_feature(300, 450, qualifiers={
            "coverage": ["1.0"],
            "identity": ["0.95"]
        })  # Length 150
        
        impact_type, metadata = self.classifier.classify_feature_impact(
            ref_feature, aligned_feature
        )
        
        self.assertEqual(impact_type, ImpactType.EXPANDED)
        self.assertAlmostEqual(metadata["length_ratio"], 1.5)
    
    def test_classify_feature_set(self):
        """Test classification of a set of features."""
        ref_features = {
            "gene1": self._create_mock_feature(100, 200),
            "gene2": self._create_mock_feature(300, 400),
            "gene3": self._create_mock_feature(500, 600)
        }
        
        aligned_features = {
            "gene1": [self._create_mock_feature(1100, 1200, qualifiers={
                "coverage": ["0.95"],
                "identity": ["0.98"]
            })],
            "gene3": [self._create_mock_feature(1500, 1550, qualifiers={
                "coverage": ["0.5"],
                "identity": ["0.95"]
            })]
            # gene2 is missing
        }
        
        results = self.classifier.classify_feature_set(ref_features, aligned_features)
        
        self.assertEqual(len(results), 3)
        self.assertEqual(results["gene1"][0], ImpactType.PRESENT)
        self.assertEqual(results["gene2"][0], ImpactType.ABSENT)
        self.assertEqual(results["gene3"][0], ImpactType.TRUNCATED)


class TestVariantDetector(unittest.TestCase):
    """Tests for the VariantDetector module."""
    
    def setUp(self):
        self.detector = VariantDetector()
    
    def _create_mock_feature(self, start, end, strand=1, qualifiers=None):
        """Create a mock SeqFeature for testing."""
        feature = MagicMock(spec=SeqFeature)
        feature.location = FeatureLocation(start, end, strand=strand)
        feature.qualifiers = qualifiers or {}
        return feature
    
    def test_detect_snps(self):
        """Test detection of SNPs."""
        ref_feature = self._create_mock_feature(2, 12)
        aln_feature = self._create_mock_feature(102, 112, qualifiers={"cigar": ["10M"]})
        
        # Create sequences with 2 SNPs at positions 5 and 8
        ref_seq = "ACGTACGTACGTACGT"
        aln_seq = "ACGTGCGTTCGTACGT"
        
        variants = self.detector.detect_variants(
            ref_feature, ref_seq, aln_feature, aln_seq
        )
        
        # Should detect 2 SNPs
        self.assertEqual(len(variants), 2)
        self.assertEqual(variants[0].variant_type, VariantType.SNP)
        self.assertEqual(variants[0].position, 7)  # 2 (start) + 5 (offset)
        self.assertEqual(variants[0].reference, "A")
        self.assertEqual(variants[0].alternate, "G")
        
        self.assertEqual(variants[1].variant_type, VariantType.SNP)
        self.assertEqual(variants[1].position, 10)  # 2 (start) + 8 (offset)
        self.assertEqual(variants[1].reference, "A")
        self.assertEqual(variants[1].alternate, "T")
    
    def test_detect_insertion(self):
        """Test detection of an insertion."""
        ref_feature = self._create_mock_feature(2, 7)
        aln_feature = self._create_mock_feature(102, 110, qualifiers={"cigar": ["5M3I2M"]})
        
        # Create sequences with an insertion of "AAA" after the 5th base
        ref_seq = "ACGTACGT"
        aln_seq = "ACGTAAAAACGT"
        
        variants = self.detector.detect_variants(
            ref_feature, ref_seq, aln_feature, aln_seq
        )
        
        # Should detect 1 insertion
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0].variant_type, VariantType.INSERTION)
        self.assertEqual(variants[0].position, 7)  # 2 (start) + 5 (offset)
        self.assertEqual(variants[0].reference, "")
        self.assertEqual(variants[0].alternate, "AAA")
        self.assertEqual(variants[0].length, 3)
    
    def test_detect_deletion(self):
        """Test detection of a deletion."""
        ref_feature = self._create_mock_feature(2, 12)
        aln_feature = self._create_mock_feature(102, 107, qualifiers={"cigar": ["5M5D"]})
        
        # Create sequences with a deletion of "ACGTA" after the 5th base
        ref_seq = "ACGTAACGTACGT"
        aln_seq = "ACGTACGT"
        
        variants = self.detector.detect_variants(
            ref_feature, ref_seq, aln_feature, aln_seq
        )
        
        # Should detect 1 deletion
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0].variant_type, VariantType.DELETION)
        self.assertEqual(variants[0].position, 7)  # 2 (start) + 5 (offset)
        self.assertEqual(variants[0].reference, "ACGTA")
        self.assertEqual(variants[0].alternate, "")
        self.assertEqual(variants[0].length, 5)
    
    def test_filter_variants(self):
        """Test filtering variants by quality."""
        variants = [
            Variant(VariantType.SNP, 100, "A", "G", 1, 60.0),
            Variant(VariantType.SNP, 200, "C", "T", 1, 20.0),
            Variant(VariantType.INSERTION, 300, "", "AGT", 3, 40.0)
        ]
        
        # Set min_quality to 30
        filtered = self.detector.filter_variants(variants, 30.0)
        
        # Should keep 2 variants
        self.assertEqual(len(filtered), 2)
        self.assertEqual(filtered[0].position, 100)  # High quality SNP
        self.assertEqual(filtered[1].position, 300)  # Medium quality insertion
    
    def test_summarize_variants(self):
        """Test summarizing variants by type."""
        variants = [
            Variant(VariantType.SNP, 100, "A", "G", 1, 60.0),
            Variant(VariantType.SNP, 200, "C", "T", 1, 20.0),
            Variant(VariantType.INSERTION, 300, "", "AGT", 3, 40.0),
            Variant(VariantType.DELETION, 400, "TGA", "", 3, 50.0),
            Variant(VariantType.SNP, 500, "G", "A", 1, 70.0)
        ]
        
        summary = self.detector.summarize_variants(variants)
        
        self.assertEqual(summary["snp"], 3)
        self.assertEqual(summary["insertion"], 1)
        self.assertEqual(summary["deletion"], 1)
        self.assertEqual(summary["complex"], 0)


class TestFeatureReconciler(unittest.TestCase):
    """Tests for the FeatureReconciler module."""
    
    def setUp(self):
        self.reconciler = FeatureReconciler(tolerance=5)
    
    def _create_mock_feature(self, start, end, strand=1, qualifiers=None):
        """Create a mock SeqFeature for testing."""
        feature = MagicMock(spec=SeqFeature)
        feature.location = FeatureLocation(start, end, strand=strand)
        feature.qualifiers = qualifiers or {}
        return feature
    
    def test_needs_reconciliation_child_within_parent(self):
        """Test that a child within parent doesn't need reconciliation."""
        parent = self._create_mock_feature(100, 200)
        child = self._create_mock_feature(120, 180)
        
        needs_reconciliation, reason = self.reconciler._needs_reconciliation(parent, child)
        
        self.assertFalse(needs_reconciliation)
        self.assertEqual(reason, "No reconciliation needed")
    
    def test_needs_reconciliation_child_outside_parent(self):
        """Test that a child outside parent needs reconciliation."""
        parent = self._create_mock_feature(100, 200)
        child = self._create_mock_feature(250, 300)
        
        needs_reconciliation, reason = self.reconciler._needs_reconciliation(parent, child)
        
        self.assertTrue(needs_reconciliation)
        self.assertEqual(reason, "Child completely outside parent")
    
    def test_needs_reconciliation_child_overlaps_start(self):
        """Test that a child overlapping parent's start needs reconciliation."""
        parent = self._create_mock_feature(100, 200)
        child = self._create_mock_feature(90, 150)
        
        needs_reconciliation, reason = self.reconciler._needs_reconciliation(parent, child)
        
        self.assertTrue(needs_reconciliation)
        self.assertEqual(reason, "Child starts before parent")
    
    def test_within_tolerance(self):
        """Test that a child just slightly outside parent is within tolerance."""
        parent = self._create_mock_feature(100, 200)
        # Just 3bp outside, within tolerance (5bp)
        child = self._create_mock_feature(97, 150)
        
        needs_reconciliation, reason = self.reconciler._needs_reconciliation(parent, child)
        
        self.assertFalse(needs_reconciliation)
        self.assertEqual(reason, "Within tolerance")
    
    def test_reconcile_adjust_child(self):
        """Test reconciliation by adjusting the child."""
        parent = self._create_mock_feature(100, 200)
        child = self._create_mock_feature(90, 150)
        
        result = self.reconciler._reconcile_feature(
            "parent1", parent, "child1", child, "Child starts before parent"
        )
        
        self.assertEqual(result.strategy, ReconciliationStrategy.ADJUST_CHILD)
        self.assertEqual(result.reconciled_feature.location.start, 100)
        self.assertEqual(result.reconciled_feature.location.end, 150)
        self.assertGreater(result.confidence, 0.5)  # Should be reasonably confident
    
    def test_reconcile_feature_hierarchy(self):
        """Test reconciliation of a feature hierarchy."""
        # Create a simple feature graph
        feature_graph = {
            "gene1": ["exon1", "exon2"],
            "gene2": ["exon3"]
        }
        
        # Create aligned features, with exon2 extending beyond gene1
        aligned_features = {
            "gene1": self._create_mock_feature(100, 200),
            "exon1": self._create_mock_feature(120, 150),
            "exon2": self._create_mock_feature(160, 220),  # Extends beyond gene1
            "gene2": self._create_mock_feature(300, 400),
            "exon3": self._create_mock_feature(320, 380)
        }
        
        results = self.reconciler.reconcile_feature_hierarchy(feature_graph, aligned_features)
        
        # Check reconciliation of exon2
        self.assertEqual(len(results["exon2"]), 1)
        self.assertEqual(results["exon2"][0].strategy, ReconciliationStrategy.ADJUST_CHILD)
        self.assertEqual(results["exon2"][0].reconciled_feature.location.end, 200)
        
        # exon3 should not need reconciliation
        self.assertEqual(len(results["exon3"]), 0)


class TestSummaryGenerator(unittest.TestCase):
    """Tests for the SummaryGenerator module."""
    
    def setUp(self):
        self.generator = SummaryGenerator()
    
    def _create_mock_feature(self, start, end, strand=1, type="gene", qualifiers=None):
        """Create a mock SeqFeature for testing."""
        feature = MagicMock(spec=SeqFeature)
        feature.location = FeatureLocation(start, end, strand=strand)
        feature.qualifiers = qualifiers or {}
        feature.type = type
        return feature
    
    def test_generate_feature_summary(self):
        """Test generating a summary for a single feature."""
        feature = self._create_mock_feature(100, 200, type="gene")
        
        impact_result = (ImpactType.MODIFIED, {"identity": 0.75, "coverage": 0.9})
        
        variants = [
            Variant(VariantType.SNP, 100, "A", "G", 1, 60.0),
            Variant(VariantType.INSERTION, 150, "", "AAA", 3, 40.0)
        ]
        
        reconciliation = ReconciliationResult(
            strategy=ReconciliationStrategy.ADJUST_CHILD,
            original_feature=feature,
            reconciled_feature=feature,
            parent_id="parent1",
            description="Adjusted to fit parent",
            confidence=0.8
        )
        
        summary = self.generator.generate_feature_summary(
            "gene1", feature, "path1",
            impact_result=impact_result,
            variants=variants,
            reconciliation=reconciliation,
            parent_features=["parent1"],
            child_features=["child1", "child2"]
        )
        
        self.assertEqual(summary.feature_id, "gene1")
        self.assertEqual(summary.feature_type, "gene")
        self.assertEqual(summary.impact_type, ImpactType.MODIFIED)
        self.assertEqual(len(summary.variants), 2)
        self.assertEqual(summary.reconciliation, reconciliation)
        self.assertEqual(summary.parent_features, ["parent1"])
        self.assertEqual(summary.child_features, ["child1", "child2"])
        self.assertAlmostEqual(summary.sequence_identity, 0.75)
        self.assertAlmostEqual(summary.coverage, 0.9)
        self.assertEqual(summary.path_id, "path1")
        self.assertEqual(summary.location, (100, 200, 1))
    
    def test_generate_path_summary(self):
        """Test generating a summary for all features in a path."""
        # Create several feature summaries
        gene1_summary = FeatureSummary(
            feature_id="gene1",
            feature_type="gene",
            impact_type=ImpactType.PRESENT,
            variants=[Variant(VariantType.SNP, 100, "A", "G", 1, 60.0)],
            sequence_identity=0.95,
            coverage=0.98,
            path_id="path1",
            location=(100, 200, 1)
        )
        
        gene2_summary = FeatureSummary(
            feature_id="gene2",
            feature_type="gene",
            impact_type=ImpactType.ABSENT,
            variants=[],
            sequence_identity=0,
            coverage=0,
            path_id="path1",
            location=(300, 400, 1)
        )
        
        gene3_summary = FeatureSummary(
            feature_id="gene3",
            feature_type="gene",
            impact_type=ImpactType.MODIFIED,
            variants=[
                Variant(VariantType.SNP, 500, "C", "T", 1, 70.0),
                Variant(VariantType.DELETION, 550, "AAA", "", 3, 50.0)
            ],
            reconciliation=ReconciliationResult(
                strategy=ReconciliationStrategy.ADJUST_PARENT,
                original_feature=MagicMock(),
                reconciled_feature=MagicMock(),
                parent_id="parent1",
                description="Adjust parent",
                confidence=0.7
            ),
            sequence_identity=0.85,
            coverage=0.9,
            path_id="path1",
            location=(500, 600, 1)
        )
        
        feature_summaries = {
            "gene1": gene1_summary,
            "gene2": gene2_summary,
            "gene3": gene3_summary
        }
        
        summary = self.generator.generate_path_summary("path1", feature_summaries)
        
        # Check the summary statistics
        self.assertEqual(summary.feature_count, 3)
        self.assertEqual(summary.feature_by_impact["present"], 1)
        self.assertEqual(summary.feature_by_impact["absent"], 1)
        self.assertEqual(summary.feature_by_impact["modified"], 1)
        
        self.assertEqual(summary.variant_counts["snp"], 2)
        self.assertEqual(summary.variant_counts["deletion"], 1)
        
        self.assertEqual(summary.reconciliation_counts["adjust_parent"], 1)
        
        # Check that all features are included
        self.assertEqual(len(summary.feature_summaries), 3)
    
    def test_export_summary_tsv(self):
        """Test exporting a summary to TSV."""
        # Create a simple analysis summary
        feature_summaries = {
            "gene1": FeatureSummary(
                feature_id="gene1",
                feature_type="gene",
                impact_type=ImpactType.PRESENT,
                variants=[Variant(VariantType.SNP, 100, "A", "G", 1, 60.0)],
                sequence_identity=0.95,
                coverage=0.98,
                path_id="path1",
                location=(100, 200, 1)
            )
        }
        
        summary = AnalysisSummary(
            path_id="path1",
            feature_count=1,
            feature_by_impact={"present": 1},
            variant_counts={"snp": 1},
            reconciliation_counts={},
            feature_summaries=feature_summaries
        )
        
        # Use a temporary file
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp:
            tmp_path = tmp.name
            
            # Export the summary
            self.generator.export_summary_tsv(summary, tmp_path)
            
            # Read the file and check the content
            with open(tmp_path, 'r') as f:
                content = f.read()
                
                self.assertIn("feature_id\tfeature_type\tpath_id\timpact\tidentity\tcoverage", content)
                self.assertIn("gene1\tgene\tpath1\tpresent\t0.95\t0.98", content)
                self.assertIn("snp:100", content)
        
        # Clean up
        os.unlink(tmp_path)


if __name__ == '__main__':
    unittest.main()
