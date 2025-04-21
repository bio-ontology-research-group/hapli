"""
Feature relationship builder that creates a graph of parent-child relationships.
Uses a graph structure to represent feature relationships from GFF3 data.
"""
import logging
from typing import Dict, List, Set, Optional, Tuple, Any
import networkx as nx
from Bio.SeqFeature import SeqFeature

class FeatureGraph:
    """
    Builds and manages a graph of feature relationships from GFF3 data.
    Uses NetworkX for the graph structure.
    """
    
    def __init__(self):
        """Initialize the feature relationship graph."""
        self.logger = logging.getLogger(__name__)
        self.graph = nx.DiGraph()  # Directed graph for parent-child relationships
        self.orphans = set()  # Features without parents
        
    def build_from_features(self, features: Dict[str, SeqFeature]) -> nx.DiGraph:
        """
        Build a graph of feature relationships.
        
        Args:
            features: Dictionary of features keyed by feature ID
            
        Returns:
            NetworkX DiGraph representing the feature relationships
        """
        self.logger.info("Building feature relationship graph")
        self.graph.clear()
        self.orphans.clear()
        
        # Add all features as nodes
        for feature_id, feature in features.items():
            self.graph.add_node(feature_id, feature=feature)
            
        # Add edges for parent-child relationships
        for feature_id, feature in features.items():
            if 'Parent' in feature.qualifiers:
                for parent_id in feature.qualifiers['Parent']:
                    if parent_id in features:
                        self.graph.add_edge(parent_id, feature_id)
                    else:
                        self.logger.warning(f"Feature {feature_id} references non-existent parent {parent_id}")
            else:
                self.orphans.add(feature_id)
                
        self.logger.info(f"Built feature graph with {self.graph.number_of_nodes()} nodes, "
                         f"{self.graph.number_of_edges()} edges, and {len(self.orphans)} orphans")
        return self.graph
    
    def get_children(self, feature_id: str) -> List[str]:
        """
        Get all immediate children of a feature.
        
        Args:
            feature_id: ID of the feature
            
        Returns:
            List of child feature IDs
        """
        if feature_id not in self.graph:
            return []
        return list(self.graph.successors(feature_id))
    
    def get_parents(self, feature_id: str) -> List[str]:
        """
        Get all immediate parents of a feature.
        
        Args:
            feature_id: ID of the feature
            
        Returns:
            List of parent feature IDs
        """
        if feature_id not in self.graph:
            return []
        return list(self.graph.predecessors(feature_id))
        
    def get_ancestors(self, feature_id: str) -> Set[str]:
        """
        Get all ancestors (recursive parents) of a feature.
        
        Args:
            feature_id: ID of the feature
            
        Returns:
            Set of ancestor feature IDs
        """
        if feature_id not in self.graph:
            return set()
            
        ancestors = set()
        queue = list(self.graph.predecessors(feature_id))
        while queue:
            current = queue.pop(0)
            if current not in ancestors:
                ancestors.add(current)
                queue.extend(self.graph.predecessors(current))
        return ancestors
    
    def get_descendants(self, feature_id: str) -> Set[str]:
        """
        Get all descendants (recursive children) of a feature.
        
        Args:
            feature_id: ID of the feature
            
        Returns:
            Set of descendant feature IDs
        """
        if feature_id not in self.graph:
            return set()
            
        return set(nx.descendants(self.graph, feature_id))
    
    def get_orphans(self) -> Set[str]:
        """
        Get all orphan features (features without parents).
        
        Returns:
            Set of orphan feature IDs
        """
        return self.orphans
    
    def get_feature_subgraph(self, feature_id: str) -> nx.DiGraph:
        """
        Get a subgraph containing a feature and all its descendants.
        
        Args:
            feature_id: ID of the root feature
            
        Returns:
            NetworkX DiGraph of the feature's hierarchy
        """
        if feature_id not in self.graph:
            return nx.DiGraph()
            
        descendants = self.get_descendants(feature_id)
        descendants.add(feature_id)
        return self.graph.subgraph(descendants)
    
    def get_feature_by_id(self, feature_id: str) -> Optional[SeqFeature]:
        """
        Get a feature by its ID.
        
        Args:
            feature_id: ID of the feature
            
        Returns:
            Feature object, or None if not found
        """
        if feature_id in self.graph:
            return self.graph.nodes[feature_id].get('feature')
        return None
