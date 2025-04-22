"""
RDF formatter for annotation reports.

This module provides the RDFFormatter class for generating
Resource Description Framework (RDF) reports in various serializations.
"""

import os
import logging
from typing import Dict, Any, List, Optional
from datetime import datetime
import uuid

import rdflib
from rdflib import Graph, Literal, URIRef, Namespace, BNode
from rdflib.namespace import RDF, RDFS, XSD, DC, DCTERMS
from pyshex import ShExEvaluator

from src.analysis.summary_generator import AnalysisSummary

logger = logging.getLogger(__name__)

# Define namespaces
HAPLO = Namespace("http://example.org/haplo/")
FALDO = Namespace("http://biohackathon.org/resource/faldo#")
SO = Namespace("http://purl.obolibrary.org/obo/so#")
VG = Namespace("http://example.org/vg/")

class RDFFormatter:
    """
    Formats analysis results as RDF in various serializations.
    
    This class converts AnalysisSummary objects into RDF format
    with support for multiple serializations (Turtle, XML, JSON-LD, N-Triples).
    """
    
    def __init__(self):
        """Initialize the RDF formatter."""
        self.format_map = {
            'turtle': 'turtle',
            'xml': 'xml',
            'json-ld': 'json-ld',
            'ntriples': 'nt',
            'n3': 'n3'
        }
        
    def format(self, summary: AnalysisSummary, output_file: str, rdf_format: str = 'turtle') -> str:
        """
        Format an analysis summary as RDF.
        
        Args:
            summary: Analysis summary to format
            output_file: Path to output file
            rdf_format: RDF serialization format
            
        Returns:
            Path to the generated file
            
        Raises:
            ValueError: If rdf_format is not supported
        """
        if rdf_format not in self.format_map:
            raise ValueError(f"Unsupported RDF format: {rdf_format}")
            
        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        
        # Create RDF graph
        g = Graph()
        
        # Bind namespaces
        g.bind('haplo', HAPLO)
        g.bind('faldo', FALDO)
        g.bind('so', SO)
        g.bind('vg', VG)
        g.bind('dc', DC)
        g.bind('dcterms', DCTERMS)
        
        # Add report metadata
        report_uri = URIRef(f"{HAPLO}report/{uuid.uuid4()}")
        g.add((report_uri, RDF.type, HAPLO.AnnotationReport))
        g.add((report_uri, DCTERMS.created, Literal(datetime.now().isoformat(), datatype=XSD.dateTime)))
        g.add((report_uri, HAPLO.pathId, Literal(summary.path_id)))
        g.add((report_uri, HAPLO.referenceId, Literal(summary.path_id.split('_')[0] if '_' in summary.path_id else "reference")))
        g.add((report_uri, HAPLO.featureCount, Literal(summary.feature_count, datatype=XSD.integer)))
        
        # Add path data
        path_uri = URIRef(f"{VG}path/{summary.path_id}")
        g.add((path_uri, RDF.type, VG.Path))
        g.add((path_uri, RDFS.label, Literal(summary.path_id)))
        g.add((report_uri, HAPLO.analyzedPath, path_uri))
        
        # Add feature data
        for feature_id, feature_summary in summary.feature_summaries.items():
            # Create feature URI
            feature_uri = URIRef(f"{HAPLO}feature/{feature_id}")
            g.add((feature_uri, RDF.type, HAPLO.GenomicFeature))
            g.add((feature_uri, RDFS.label, Literal(feature_id)))
            g.add((feature_uri, HAPLO.featureType, Literal(feature_summary.feature_type)))
            
            # Add to report
            g.add((report_uri, HAPLO.hasFeature, feature_uri))
            
            # Add impact data
            if feature_summary.impact_type:
                impact_node = BNode()
                g.add((impact_node, RDF.type, HAPLO.FeatureImpact))
                g.add((impact_node, HAPLO.impactType, Literal(feature_summary.impact_type.value)))
                
                if hasattr(feature_summary, 'sequence_identity') and feature_summary.sequence_identity is not None:
                    g.add((impact_node, HAPLO.sequenceIdentity, 
                          Literal(feature_summary.sequence_identity, datatype=XSD.decimal)))
                    
                if hasattr(feature_summary, 'coverage') and feature_summary.coverage is not None:
                    g.add((impact_node, HAPLO.coverage, 
                          Literal(feature_summary.coverage, datatype=XSD.decimal)))
                
                g.add((feature_uri, HAPLO.impact, impact_node))
            
            # Add location data using FALDO
            if feature_summary.location:
                location_node = BNode()
                g.add((location_node, RDF.type, FALDO.Region))
                
                # Start position
                start_node = BNode()
                g.add((start_node, RDF.type, FALDO.Position))
                g.add((start_node, FALDO.position, Literal(feature_summary.location[0], datatype=XSD.integer)))
                g.add((location_node, FALDO.begin, start_node))
                
                # End position
                end_node = BNode()
                g.add((end_node, RDF.type, FALDO.Position))
                g.add((end_node, FALDO.position, Literal(feature_summary.location[1], datatype=XSD.integer)))
                g.add((location_node, FALDO.end, end_node))
                
                # Strand
                if feature_summary.location[2] == 1:
                    g.add((start_node, RDF.type, FALDO.ForwardStrandPosition))
                    g.add((end_node, RDF.type, FALDO.ForwardStrandPosition))
                else:
                    g.add((start_node, RDF.type, FALDO.ReverseStrandPosition))
                    g.add((end_node, RDF.type, FALDO.ReverseStrandPosition))
                    
                g.add((feature_uri, FALDO.location, location_node))
            
            # Add variant data
            if hasattr(feature_summary, 'variants') and feature_summary.variants:
                for i, variant in enumerate(feature_summary.variants):
                    variant_uri = URIRef(f"{HAPLO}feature/{feature_id}/variant/{i}")
                    g.add((variant_uri, RDF.type, HAPLO.SequenceVariant))
                    g.add((variant_uri, HAPLO.variantType, Literal(variant.variant_type.value)))
                    g.add((variant_uri, HAPLO.position, Literal(variant.position, datatype=XSD.integer)))
                    g.add((variant_uri, HAPLO.referenceAllele, Literal(variant.reference)))
                    g.add((variant_uri, HAPLO.alternateAllele, Literal(variant.alternate)))
                    g.add((variant_uri, HAPLO.variantLength, Literal(variant.length, datatype=XSD.integer)))
                    if variant.quality is not None:
                        g.add((variant_uri, HAPLO.variantQuality, Literal(variant.quality, datatype=XSD.decimal)))
                    g.add((feature_uri, HAPLO.hasVariant, variant_uri))
            
            # Add parent-child relationships
            if hasattr(feature_summary, 'parent_features') and feature_summary.parent_features:
                for parent_id in feature_summary.parent_features:
                    parent_uri = URIRef(f"{HAPLO}feature/{parent_id}")
                    g.add((feature_uri, HAPLO.hasParent, parent_uri))
        
        # Serialize to the specified format
        rdf_format_str = self.format_map[rdf_format]
        g.serialize(destination=output_file, format=rdf_format_str)
        
        logger.info(f"Generated RDF report in {rdf_format} format at {output_file}")
        return output_file
        
    def format_comparative(self, comparative_data: Dict[str, Any], output_file: str, rdf_format: str = 'turtle') -> str:
        """
        Format comparative analysis as RDF.
        
        Args:
            comparative_data: Comparative data structure
            output_file: Path to output file
            rdf_format: RDF serialization format
            
        Returns:
            Path to the generated file
            
        Raises:
            ValueError: If rdf_format is not supported
        """
        if rdf_format not in self.format_map:
            raise ValueError(f"Unsupported RDF format: {rdf_format}")
            
        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        
        # Create RDF graph
        g = Graph()
        
        # Bind namespaces
        g.bind('haplo', HAPLO)
        g.bind('faldo', FALDO)
        g.bind('so', SO)
        g.bind('vg', VG)
        g.bind('dc', DC)
        g.bind('dcterms', DCTERMS)
        
        # Add report metadata
        report_uri = URIRef(f"{HAPLO}comparative-report/{uuid.uuid4()}")
        g.add((report_uri, RDF.type, HAPLO.ComparativeAnnotationReport))
        g.add((report_uri, DCTERMS.created, Literal(comparative_data['metadata']['timestamp'], datatype=XSD.dateTime)))
        g.add((report_uri, HAPLO.referenceId, Literal(comparative_data['metadata']['reference_id'])))
        g.add((report_uri, HAPLO.featureCount, Literal(comparative_data['metadata']['feature_count'], datatype=XSD.integer)))
        
        # Add path data
        for path_id in comparative_data['metadata']['paths']:
            path_uri = URIRef(f"{VG}path/{path_id}")
            g.add((path_uri, RDF.type, VG.Path))
            g.add((path_uri, RDFS.label, Literal(path_id)))
            g.add((report_uri, HAPLO.comparesPath, path_uri))
            
        # Add comparative feature data
        for feature_id, feature_data in comparative_data['features'].items():
            # Create shared feature node
            feature_uri = URIRef(f"{HAPLO}feature/{feature_id}")
            g.add((feature_uri, RDF.type, HAPLO.GenomicFeature))
            g.add((feature_uri, RDFS.label, Literal(feature_id)))
            g.add((report_uri, HAPLO.hasComparativeFeature, feature_uri))
            
            # Add path-specific data
            for path_id, path_data in feature_data['paths'].items():
                if not path_data or path_data.get('impact_type') == 'ABSENT':
                    continue
                    
                path_uri = URIRef(f"{VG}path/{path_id}")
                
                # Create path-specific feature occurrence
                occurrence_uri = URIRef(f"{HAPLO}feature/{feature_id}/in/{path_id}")
                g.add((occurrence_uri, RDF.type, HAPLO.FeatureOccurrence))
                g.add((occurrence_uri, HAPLO.inPath, path_uri))
                g.add((occurrence_uri, HAPLO.hasFeature, feature_uri))
                g.add((occurrence_uri, HAPLO.impactType, Literal(path_data.get('impact_type'))))
                
                # Add metrics if available
                if 'sequence_identity' in path_data and path_data['sequence_identity'] is not None:
                    g.add((occurrence_uri, HAPLO.sequenceIdentity, 
                          Literal(path_data['sequence_identity'], datatype=XSD.decimal)))
                    
                if 'coverage' in path_data and path_data['coverage'] is not None:
                    g.add((occurrence_uri, HAPLO.coverage, 
                          Literal(path_data['coverage'], datatype=XSD.decimal)))
                
                # Add location if available
                if path_data.get('start') is not None and path_data.get('end') is not None:
                    location_node = BNode()
                    g.add((location_node, RDF.type, FALDO.Region))
                    
                    # Start position
                    start_node = BNode()
                    g.add((start_node, RDF.type, FALDO.Position))
                    g.add((start_node, FALDO.position, Literal(path_data['start'], datatype=XSD.integer)))
                    g.add((location_node, FALDO.begin, start_node))
                    
                    # End position
                    end_node = BNode()
                    g.add((end_node, RDF.type, FALDO.Position))
                    g.add((end_node, FALDO.position, Literal(path_data['end'], datatype=XSD.integer)))
                    g.add((location_node, FALDO.end, end_node))
                    
                    # Strand
                    if path_data.get('strand') == '+' or path_data.get('strand') == 1:
                        g.add((start_node, RDF.type, FALDO.ForwardStrandPosition))
                        g.add((end_node, RDF.type, FALDO.ForwardStrandPosition))
                    elif path_data.get('strand') == '-' or path_data.get('strand') == -1:
                        g.add((start_node, RDF.type, FALDO.ReverseStrandPosition))
                        g.add((end_node, RDF.type, FALDO.ReverseStrandPosition))
                        
                    g.add((occurrence_uri, FALDO.location, location_node))
                
                # Add variant summary if available
                if path_data.get('variants'):
                    variants_node = BNode()
                    g.add((variants_node, RDF.type, HAPLO.VariantSummary))
                    
                    for var_type, count in path_data['variants'].items():
                        if var_type == 'SNP':
                            g.add((variants_node, HAPLO.snpCount, Literal(count, datatype=XSD.integer)))
                        elif var_type == 'INSERTION':
                            g.add((variants_node, HAPLO.insertionCount, Literal(count, datatype=XSD.integer)))
                        elif var_type == 'DELETION':
                            g.add((variants_node, HAPLO.deletionCount, Literal(count, datatype=XSD.integer)))
                        elif var_type == 'COMPLEX':
                            g.add((variants_node, HAPLO.complexCount, Literal(count, datatype=XSD.integer)))
                    
                    g.add((occurrence_uri, HAPLO.variantSummary, variants_node))
        
        # Serialize to the specified format
        rdf_format_str = self.format_map[rdf_format]
        g.serialize(destination=output_file, format=rdf_format_str)
        
        logger.info(f"Generated comparative RDF report in {rdf_format} format at {output_file}")
        return output_file
        
    def validate(self, rdf_file: str, shex_file: str) -> bool:
        """
        Validate an RDF report against a ShEx schema.
        
        Args:
            rdf_file: Path to RDF report file
            shex_file: Path to ShEx schema file
            
        Returns:
            True if validation succeeded, False otherwise
        """
        try:
            # Load RDF data
            g = Graph()
            g.parse(rdf_file)
            
            # Find report node(s)
            report_nodes = []
            for s, p, o in g.triples((None, RDF.type, HAPLO.AnnotationReport)):
                report_nodes.append(s)
            for s, p, o in g.triples((None, RDF.type, HAPLO.ComparativeAnnotationReport)):
                report_nodes.append(s)
                
            if not report_nodes:
                logger.error("No report nodes found in RDF file")
                return False
                
            # Run validation
            validator = ShExEvaluator(g, shex_file)
            for node in report_nodes:
                result = validator.validate(node, HAPLO.AnnotationReport)
                if not result:
                    logger.error(f"Validation failed for node {node}")
                    return False
                    
            logger.info(f"RDF report {rdf_file} successfully validated against {shex_file}")
            return True
            
        except Exception as e:
            logger.error(f"Error validating RDF file: {str(e)}")
            return False
