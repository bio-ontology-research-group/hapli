import unittest
from Bio.SeqFeature import SeqFeature, FeatureLocation
from src.analysis.impact_classifier import ImpactClassifier, ImpactType
from src.analysis.variant_detector import Variant, VariantType

class TestImpactClassifier(unittest.TestCase):
    def setUp(self):
        self.classifier = ImpactClassifier()
        self.ref_feature = SeqFeature(
            FeatureLocation(1000, 2000),
            type="gene",
            qualifiers={"ID": ["gene1"]}
        )

    def create_aligned_feature(self, coverage=1.0, identity=1.0, length=1000, variants=None, sub_features=None):
        feature = SeqFeature(
            FeatureLocation(0, length),
            type="gene",
            qualifiers={
                "coverage": [str(coverage)],
                "identity": [str(identity)],
                "alignment_score": ["60"]
            }
        )
        if variants:
            feature.qualifiers["variants"] = variants
        if sub_features:
            feature.sub_features = sub_features
        return feature

    def test_absent_feature(self):
        impact_type, metadata = self.classifier.classify_feature_impact(
            self.ref_feature, None
        )
        self.assertEqual(impact_type, ImpactType.ABSENT)
        self.assertEqual(metadata["reason"], "No alignment found")

    def test_present_feature(self):
        aligned = self.create_aligned_feature(coverage=0.95, identity=0.98)
        impact_type, metadata = self.classifier.classify_feature_impact(
            self.ref_feature, aligned
        )
        self.assertEqual(impact_type, ImpactType.PRESENT)
        self.assertAlmostEqual(metadata["coverage"], 0.95)
        self.assertAlmostEqual(metadata["identity"], 0.98)

    def test_modified_feature(self):
        aligned = self.create_aligned_feature(coverage=0.90, identity=0.85)
        impact_type, _ = self.classifier.classify_feature_impact(
            self.ref_feature, aligned
        )
        self.assertEqual(impact_type, ImpactType.MODIFIED)

    def test_truncated_feature(self):
        aligned = self.create_aligned_feature(coverage=0.95, identity=0.98, length=800)
        impact_type, metadata = self.classifier.classify_feature_impact(
            self.ref_feature, aligned
        )
        self.assertEqual(impact_type, ImpactType.TRUNCATED)
        self.assertAlmostEqual(metadata["length_ratio"], 0.8)

    def test_expanded_feature(self):
        aligned = self.create_aligned_feature(coverage=0.95, identity=0.98, length=1200)
        impact_type, metadata = self.classifier.classify_feature_impact(
            self.ref_feature, aligned
        )
        self.assertEqual(impact_type, ImpactType.EXPANDED)
        self.assertAlmostEqual(metadata["length_ratio"], 1.2)

    def test_fragmented_feature(self):
        sub_features = [
            SeqFeature(FeatureLocation(0, 500)),
            SeqFeature(FeatureLocation(600, 1000))
        ]
        aligned = self.create_aligned_feature(sub_features=sub_features)
        impact_type, metadata = self.classifier.classify_feature_impact(
            self.ref_feature, aligned
        )
        self.assertEqual(impact_type, ImpactType.FRAGMENTED)
        self.assertEqual(metadata["fragments"], 2)

    def test_inverted_feature(self):
        variant = Variant(
            variant_type=VariantType.INVERSION,
            position=1500,
            reference="N",
            alternate="N",
            length=500,
            metadata={"new_position": 2000}
        )
        aligned = self.create_aligned_feature(variants=[variant])
        impact_type, metadata = self.classifier.classify_feature_impact(
            self.ref_feature, aligned
        )
        self.assertEqual(impact_type, ImpactType.INVERTED)
        self.assertEqual(metadata["inversion_length"], 500)

    def test_feature_set_classification(self):
        ref_features = {
            "gene1": self.ref_feature,
            "gene2": SeqFeature(FeatureLocation(3000, 4000), type="gene")
        }
        
        aligned_features = {
            "gene1": [self.create_aligned_feature(coverage=0.95, identity=0.98)],
            "gene2": []
        }
        
        results = self.classifier.classify_feature_set(ref_features, aligned_features)
        
        self.assertEqual(results["gene1"][0], ImpactType.PRESENT)
        self.assertEqual(results["gene2"][0], ImpactType.ABSENT)

if __name__ == '__main__':
    unittest.main()
