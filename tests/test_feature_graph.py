import unittest
import tempfile
import os
from src.parsers.gff_parser import GFF3Parser
from src.parsers.feature_graph import FeatureGraph
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

class TestFeatureGraph(unittest.TestCase):
    """Tests for the Feature Relationship Graph."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.graph = FeatureGraph()
        self.data_dir = "data"
        self.example_gff = os.path.join(self.data_dir, "example.gff3")
        
        # Create a test GFF3 file with parent-child relationships
        self.temp_dir = tempfile.TemporaryDirectory()
        self.test_gff = os.path.join(self.temp_dir.name, "test.gff3")
        with open(self.test_gff, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("##sequence-region ref_chr1 1 1000\n")
            f.write("ref_chr1\t.\tgene\t1\t500\t.\t+\t.\tID=gene1\n")
            f.write("ref_chr1\t.\tmRNA\t1\t500\t.\t+\t.\tID=mRNA1;Parent=gene1\n")
            f.write("ref_chr1\t.\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=mRNA1\n")
            f.write("ref_chr1\t.\texon\t201\t300\t.\t+\t.\tID=exon2;Parent=mRNA1\n")
            f.write("ref_chr1\t.\tCDS\t1\t100\t.\t+\t0\tID=cds1;Parent=mRNA1\n")
            f.write("ref_chr1\t.\tCDS\t201\t300\t.\t+\t0\tID=cds2;Parent=mRNA1\n")
            # Add a gene without parents
            f.write("ref_chr1\t.\tgene\t600\t800\t.\t+\t.\tID=gene2\n")
            # Add an orphaned feature (no parent)
            f.write("ref_chr1\t.\tmisc_feature\t900\t950\t.\t+\t.\tID=misc1\n")
            # Add a feature with non-existent parent
            f.write("ref_chr1\t.\texon\t600\t650\t.\t+\t.\tID=exon3;Parent=nonexistent\n")
        
    def tearDown(self):
        """Clean up test fixtures."""
        self.temp_dir.cleanup()
        
    def test_build_from_features(self):
        """Test building a graph from GFF3 features."""
        # Parse the GFF3 file
        parser = GFF3Parser()
        parser.parse(self.test_gff)
        
        # Build the graph
        graph = self.graph.build_from_features(parser.features_by_id)
        
        # Check the graph
        self.assertGreater(graph.number_of_nodes(), 0)
        self.assertGreater(graph.number_of_edges(), 0)
        
        # Check orphans
        orphans = self.graph.get_orphans()
        self.assertGreaterEqual(len(orphans), 2)  # gene1, gene2, misc1
        self.assertIn('gene1', orphans)
        self.assertIn('gene2', orphans)
        self.assertIn('misc1', orphans)
        
    def test_get_children(self):
        """Test getting children of a feature."""
        # Parse the GFF3 file
        parser = GFF3Parser()
        parser.parse(self.test_gff)
        
        # Build the graph
        self.graph.build_from_features(parser.features_by_id)
        
        # Check children
        children = self.graph.get_children('gene1')
        self.assertEqual(len(children), 1)
        self.assertIn('mRNA1', children)
        
        children = self.graph.get_children('mRNA1')
        self.assertEqual(len(children), 4)  # exon1, exon2, cds1, cds2
        
        # Check non-existent feature
        children = self.graph.get_children('nonexistent')
        self.assertEqual(len(children), 0)
        
    def test_get_parents(self):
        """Test getting parents of a feature."""
        # Parse the GFF3 file
        parser = GFF3Parser()
        parser.parse(self.test_gff)
        
        # Build the graph
        self.graph.build_from_features(parser.features_by_id)
        
        # Check parents
        parents = self.graph.get_parents('mRNA1')
        self.assertEqual(len(parents), 1)
        self.assertIn('gene1', parents)
        
        parents = self.graph.get_parents('exon1')
        self.assertEqual(len(parents), 1)
        self.assertIn('mRNA1', parents)
        
        # Check orphan
        parents = self.graph.get_parents('gene1')
        self.assertEqual(len(parents), 0)
        
        # Check non-existent feature
        parents = self.graph.get_parents('nonexistent')
        self.assertEqual(len(parents), 0)
        
    def test_get_ancestors(self):
        """Test getting ancestors of a feature."""
        # Parse the GFF3 file
        parser = GFF3Parser()
        parser.parse(self.test_gff)
        
        # Build the graph
        self.graph.build_from_features(parser.features_by_id)
        
        # Check ancestors
        ancestors = self.graph.get_ancestors('exon1')
        self.assertEqual(len(ancestors), 2)
        self.assertIn('mRNA1', ancestors)
        self.assertIn('gene1', ancestors)
        
        # Check orphan
        ancestors = self.graph.get_ancestors('gene1')
        self.assertEqual(len(ancestors), 0)
        
        # Check non-existent feature
        ancestors = self.graph.get_ancestors('nonexistent')
        self.assertEqual(len(ancestors), 0)
        
    def test_get_descendants(self):
        """Test getting descendants of a feature."""
        # Parse the GFF3 file
        parser = GFF3Parser()
        parser.parse(self.test_gff)
        
        # Build the graph
        self.graph.build_from_features(parser.features_by_id)
        
        # Check descendants
        descendants = self.graph.get_descendants('gene1')
        self.assertEqual(len(descendants), 5)  # mRNA1, exon1, exon2, cds1, cds2
        
        # Check feature with no children
        descendants = self.graph.get_descendants('exon1')
        self.assertEqual(len(descendants), 0)
        
        # Check non-existent feature
        descendants = self.graph.get_descendants('nonexistent')
        self.assertEqual(len(descendants), 0)
        
    def test_get_feature_subgraph(self):
        """Test getting a subgraph for a feature and its descendants."""
        # Parse the GFF3 file
        parser = GFF3Parser()
        parser.parse(self.test_gff)
        
        # Build the graph
        self.graph.build_from_features(parser.features_by_id)
        
        # Check subgraph
        subgraph = self.graph.get_feature_subgraph('gene1')
        self.assertEqual(subgraph.number_of_nodes(), 6)  # gene1, mRNA1, exon1, exon2, cds1, cds2
        
        # Check feature with no children
        subgraph = self.graph.get_feature_subgraph('exon1')
        self.assertEqual(subgraph.number_of_nodes(), 1)  # just exon1
        
        # Check non-existent feature
        subgraph = self.graph.get_feature_subgraph('nonexistent')
        self.assertEqual(subgraph.number_of_nodes(), 0)
        
    def test_get_feature_by_id(self):
        """Test getting a feature by ID."""
        # Parse the GFF3 file
        parser = GFF3Parser()
        parser.parse(self.test_gff)
        
        # Build the graph
        self.graph.build_from_features(parser.features_by_id)
        
        # Check feature
        feature = self.graph.get_feature_by_id('gene1')
        self.assertIsNotNone(feature)
        self.assertEqual(feature.type, 'gene')
        
        # Check non-existent feature
        feature = self.graph.get_feature_by_id('nonexistent')
        self.assertIsNone(feature)
        
    def test_handle_orphaned_features(self):
        """Test handling of orphaned features (features with no parents)."""
        # Create test features
        features = {
            'gene1': SeqFeature(FeatureLocation(0, 100), type='gene', 
                               qualifiers={'ID': ['gene1']}),
            'mRNA1': SeqFeature(FeatureLocation(0, 100), type='mRNA', 
                               qualifiers={'ID': ['mRNA1'], 'Parent': ['gene1']}),
            'orphan1': SeqFeature(FeatureLocation(200, 300), type='misc_feature', 
                                 qualifiers={'ID': ['orphan1']}),
            'orphan2': SeqFeature(FeatureLocation(400, 500), type='misc_feature', 
                                 qualifiers={'ID': ['orphan2']}),
        }
        
        # Build the graph
        self.graph.build_from_features(features)
        
        # Check orphans
        orphans = self.graph.get_orphans()
        self.assertEqual(len(orphans), 3)  # gene1, orphan1, orphan2
        self.assertIn('gene1', orphans)
        self.assertIn('orphan1', orphans)
        self.assertIn('orphan2', orphans)

if __name__ == '__main__':
    unittest.main()
