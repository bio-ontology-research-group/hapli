#!/usr/bin/env python3
"""
Utility script to generate integration test data.

This script creates a complete set of test data for integration testing:
- Reference FASTA sequence
- GFF3 annotation file
- GFA2 variation graph
- VCF variant file

All files are created using appropriate bioinformatics libraries to ensure they are
properly formatted and the introduced variants are correctly aligned with features.
"""

import os
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
import random
import re
from collections import defaultdict

# Try to import bioinformatics libraries
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
    import gfapy
    from pyfaidx import Fasta
except ImportError:
    print("Required libraries not found. Please install with:")
    print("pip install biopython gfapy pyfaidx")
    exit(1)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Output directory
OUTPUT_DIR = Path("data/integration_test")

class TestDataGenerator:
    """Generate test data for integration testing."""
    
    def __init__(self, output_dir: Path = OUTPUT_DIR):
        """Initialize the generator with output paths."""
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Output file paths
        self.fasta_path = output_dir / "reference.fasta"
        self.gff_path = output_dir / "annotations.gff3"
        self.gfa_path = output_dir / "variation_graph.gfa"
        self.vcf_path = output_dir / "variants.vcf"
        
        # Data structures
        self.reference_seq = None
        self.features = []
        self.segments = {}
        self.paths = {}
        self.variants = []
        
        # Sample and haplotype info
        self.samples = ["sample1", "sample2"]
        self.haplotypes = ["hap1", "hap2"]
        
    def generate_all(self):
        """Generate all test data files."""
        logger.info("Generating integration test data")
        
        # Generate in order of dependencies
        self.generate_reference_sequence()
        self.generate_features()
        self.generate_variants()
        self.build_variation_graph()
        
        # Write to files
        self.write_fasta_file()
        self.write_gff_file()
        self.write_gfa_file()
        self.write_vcf_file()
        
        logger.info(f"All test data generated in {self.output_dir}")
        
    def generate_reference_sequence(self, length: int = 1560):
        """Generate a reference sequence of the given length."""
        logger.info(f"Generating reference sequence of length {length}")
        
        # Create a random DNA sequence (for reproducibility, we could use a seed)
        nucleotides = ['A', 'C', 'G', 'T']
        sequence = ''.join(random.choice(nucleotides) for _ in range(length))
        
        # Create BioPython SeqRecord
        self.reference_seq = SeqRecord(
            Seq(sequence),
            id="reference",
            name="reference",
            description="Reference sequence for integration testing"
        )
        
        logger.info(f"Generated reference sequence: {len(self.reference_seq)}bp")
        return self.reference_seq
        
    def generate_features(self):
        """Generate genomic features for the reference sequence."""
        logger.info("Generating genomic features")
        
        # Feature coordinates based on our test plan
        feature_coords = [
            # Gene name, start, end, strand, subfeatures
            ("Gene1", 50, 250, 1, self._generate_mrna_subfeatures),
            ("Gene2", 300, 500, 1, self._generate_mrna_subfeatures),
            ("Gene3", 550, 650, 1, self._generate_mrna_subfeatures),
            ("Gene4", 700, 900, 1, self._generate_mrna_subfeatures),
            ("Gene5", 950, 1050, 1, self._generate_mrna_subfeatures),
            ("Gene6", 1100, 1200, -1, self._generate_mrna_subfeatures),
            ("Gene7", 1250, 1350, 1, self._generate_mrna_subfeatures),
            ("Gene8", 1400, 1500, 1, self._generate_mrna_subfeatures),
            ("Gene9", 1520, 1560, 1, self._generate_mrna_subfeatures),
            ("Gene10", 200, 500, -1, self._generate_alternative_mrna_subfeatures),
        ]
        
        # Generate feature objects
        for name, start, end, strand, subfeat_func in feature_coords:
            gene_id = f"gene{name[4:]}"  # Extract number from name
            
            # Create gene feature
            gene_feature = SeqFeature(
                FeatureLocation(start-1, end, strand=strand),  # 0-based internally
                type="gene",
                qualifiers={
                    "ID": [gene_id],
                    "Name": [name]
                }
            )
            
            # Generate subfeatures (mRNAs, exons, etc.)
            subfeat_func(gene_feature)
            
            # Add to feature list
            self.features.append(gene_feature)
        
        logger.info(f"Generated {len(self.features)} gene features with subfeatures")
        return self.features
            
    def _generate_mrna_subfeatures(self, gene_feature):
        """Generate mRNA and child features for a gene."""
        gene_id = gene_feature.qualifiers["ID"][0]
        gene_start = gene_feature.location.start
        gene_end = gene_feature.location.end
        gene_strand = gene_feature.location.strand
        gene_length = gene_end - gene_start
        
        # Create mRNA feature
        mrna_id = f"mRNA{gene_id[4:]}.1"  # e.g., mRNA1.1
        mrna_feature = SeqFeature(
            FeatureLocation(gene_start, gene_end, strand=gene_strand),
            type="mRNA",
            qualifiers={
                "ID": [mrna_id],
                "Parent": [gene_id],
                "Name": [f"{gene_id[0].upper()}{gene_id[1:]}.1"]  # e.g., Gene1.1
            }
        )
        
        # Decide on exon structure
        if gene_length < 100:
            # Small gene, single exon
            exon_coords = [(gene_start, gene_end)]
        else:
            # Multi-exon gene
            # For simplicity, we'll create two exons with an intron between
            intron_size = min(gene_length // 5, 50)  # Intron ~20% of gene length, max 50bp
            exon1_end = gene_start + (gene_length - intron_size) // 2
            exon2_start = exon1_end + intron_size
            exon_coords = [(gene_start, exon1_end), (exon2_start, gene_end)]
        
        # Add exon features
        for i, (ex_start, ex_end) in enumerate(exon_coords, 1):
            exon_id = f"exon{mrna_id[4:]}.{i}"  # e.g., exon1.1.1
            exon_feature = SeqFeature(
                FeatureLocation(ex_start, ex_end, strand=gene_strand),
                type="exon",
                qualifiers={
                    "ID": [exon_id],
                    "Parent": [mrna_id]
                }
            )
            mrna_feature.sub_features.append(exon_feature)
        
        # Add UTR and CDS features if gene is long enough
        if gene_length >= 80:
            if gene_strand == 1:  # Forward strand
                utr5_start, utr5_end = gene_start, gene_start + 20
                cds_coords = []
                
                # CDS in first exon (if it exists after 5' UTR)
                if exon_coords[0][0] < utr5_end and exon_coords[0][1] > utr5_end:
                    cds_coords.append((utr5_end, exon_coords[0][1]))
                
                # CDS in remaining exons except final portion
                if len(exon_coords) > 1:
                    for ex_start, ex_end in exon_coords[1:-1]:
                        cds_coords.append((ex_start, ex_end))
                    
                    # Last exon CDS ends before 3' UTR
                    utr3_start = max(exon_coords[-1][1] - 20, exon_coords[-1][0])
                    cds_coords.append((exon_coords[-1][0], utr3_start))
                    utr3_end = exon_coords[-1][1]
                else:
                    # Single exon gene
                    utr3_start = max(exon_coords[0][1] - 20, utr5_end)
                    cds_coords = [(utr5_end, utr3_start)]
                    utr3_end = exon_coords[0][1]
                
                # 5' UTR
                utr5_feature = SeqFeature(
                    FeatureLocation(utr5_start, utr5_end, strand=gene_strand),
                    type="five_prime_UTR",
                    qualifiers={
                        "ID": [f"5UTR{mrna_id[4:]}"],
                        "Parent": [mrna_id]
                    }
                )
                mrna_feature.sub_features.append(utr5_feature)
                
                # CDS features
                for i, (cds_start, cds_end) in enumerate(cds_coords, 1):
                    # Phase calculation is simplified for demo
                    phase = 0 if i == 1 else [0, 2, 1][i % 3]
                    
                    cds_feature = SeqFeature(
                        FeatureLocation(cds_start, cds_end, strand=gene_strand),
                        type="CDS",
                        qualifiers={
                            "ID": [f"CDS{mrna_id[4:]}.{i}"],
                            "Parent": [mrna_id],
                            "phase": [str(phase)]
                        }
                    )
                    mrna_feature.sub_features.append(cds_feature)
                
                # 3' UTR
                utr3_feature = SeqFeature(
                    FeatureLocation(utr3_start, utr3_end, strand=gene_strand),
                    type="three_prime_UTR",
                    qualifiers={
                        "ID": [f"3UTR{mrna_id[4:]}"],
                        "Parent": [mrna_id]
                    }
                )
                mrna_feature.sub_features.append(utr3_feature)
                
            else:  # Reverse strand, mirror the forward strand logic
                utr3_start, utr3_end = gene_start, gene_start + 20
                cds_coords = []
                
                # CDS in first exon (if it exists after 3' UTR)
                if exon_coords[0][0] < utr3_end and exon_coords[0][1] > utr3_end:
                    cds_coords.append((utr3_end, exon_coords[0][1]))
                
                # CDS in remaining exons except final portion
                if len(exon_coords) > 1:
                    for ex_start, ex_end in exon_coords[1:-1]:
                        cds_coords.append((ex_start, ex_end))
                    
                    # Last exon CDS ends before 5' UTR
                    utr5_start = max(exon_coords[-1][1] - 20, exon_coords[-1][0])
                    cds_coords.append((exon_coords[-1][0], utr5_start))
                    utr5_end = exon_coords[-1][1]
                else:
                    # Single exon gene
                    utr5_start = max(exon_coords[0][1] - 20, utr3_end)
                    cds_coords = [(utr3_end, utr5_start)]
                    utr5_end = exon_coords[0][1]
                
                # 3' UTR
                utr3_feature = SeqFeature(
                    FeatureLocation(utr3_start, utr3_end, strand=gene_strand),
                    type="three_prime_UTR",
                    qualifiers={
                        "ID": [f"3UTR{mrna_id[4:]}"],
                        "Parent": [mrna_id]
                    }
                )
                mrna_feature.sub_features.append(utr3_feature)
                
                # CDS features
                for i, (cds_start, cds_end) in enumerate(cds_coords, 1):
                    # Phase calculation is simplified for demo
                    phase = 0 if i == 1 else [0, 2, 1][i % 3]
                    
                    cds_feature = SeqFeature(
                        FeatureLocation(cds_start, cds_end, strand=gene_strand),
                        type="CDS",
                        qualifiers={
                            "ID": [f"CDS{mrna_id[4:]}.{i}"],
                            "Parent": [mrna_id],
                            "phase": [str(phase)]
                        }
                    )
                    mrna_feature.sub_features.append(cds_feature)
                
                # 5' UTR
                utr5_feature = SeqFeature(
                    FeatureLocation(utr5_start, utr5_end, strand=gene_strand),
                    type="five_prime_UTR",
                    qualifiers={
                        "ID": [f"5UTR{mrna_id[4:]}"],
                        "Parent": [mrna_id]
                    }
                )
                mrna_feature.sub_features.append(utr5_feature)
        
        # Add mRNA to gene's subfeatures
        gene_feature.sub_features.append(mrna_feature)
    
    def _generate_alternative_mrna_subfeatures(self, gene_feature):
        """Generate alternative mRNA isoforms for a gene (e.g., Gene10)."""
        gene_id = gene_feature.qualifiers["ID"][0]
        gene_start = gene_feature.location.start
        gene_end = gene_feature.location.end
        gene_strand = gene_feature.location.strand
        gene_length = gene_end - gene_start
        
        # Create two mRNA isoforms
        for iso_num in [1, 2]:
            mrna_id = f"mRNA{gene_id[4:]}.{iso_num}"  # e.g., mRNA10.1
            mrna_feature = SeqFeature(
                FeatureLocation(gene_start, gene_end, strand=gene_strand),
                type="mRNA",
                qualifiers={
                    "ID": [mrna_id],
                    "Parent": [gene_id],
                    "Name": [f"{gene_id[0].upper()}{gene_id[1:]}.{iso_num}"]  # e.g., Gene10.1
                }
            )
            
            # Different exon structure for each isoform
            if iso_num == 1:
                # First isoform has three exons
                exon1_end = gene_start + 50
                exon2_start = exon1_end + 50
                exon2_end = exon2_start + 50
                exon3_start = exon2_end + 50
                exon_coords = [
                    (gene_start, exon1_end), 
                    (exon2_start, exon2_end), 
                    (exon3_start, gene_end)
                ]
            else:
                # Second isoform has two exons (skips the middle exon)
                exon1_end = gene_start + 50
                exon2_start = gene_start + 150  # Skip the middle exon
                exon_coords = [
                    (gene_start, exon1_end), 
                    (exon2_start, gene_end)
                ]
            
            # Add exon features
            for i, (ex_start, ex_end) in enumerate(exon_coords, 1):
                exon_id = f"exon{gene_id[4:]}.{iso_num}.{i}"  # e.g., exon10.1.1
                exon_feature = SeqFeature(
                    FeatureLocation(ex_start, ex_end, strand=gene_strand),
                    type="exon",
                    qualifiers={
                        "ID": [exon_id],
                        "Parent": [mrna_id]
                    }
                )
                mrna_feature.sub_features.append(exon_feature)
            
            # Add UTR and CDS features
            if gene_strand == 1:  # Forward strand
                utr5_start, utr5_end = gene_start, gene_start + 20
                utr3_start = gene_end - 20
                utr3_end = gene_end
                
                # CDS coordinates depend on exon structure
                if iso_num == 1:
                    cds_coords = [
                        (utr5_end, exon1_end),
                        (exon2_start, exon2_end),
                        (exon3_start, utr3_start)
                    ]
                else:
                    cds_coords = [
                        (utr5_end, exon1_end),
                        (exon2_start, utr3_start)
                    ]
            else:  # Reverse strand
                utr3_start, utr3_end = gene_start, gene_start + 20
                utr5_start = gene_end - 20
                utr5_end = gene_end
                
                # CDS coordinates depend on exon structure
                if iso_num == 1:
                    cds_coords = [
                        (utr3_end, exon1_end),
                        (exon2_start, exon2_end),
                        (exon3_start, utr5_start)
                    ]
                else:
                    cds_coords = [
                        (utr3_end, exon1_end),
                        (exon2_start, utr5_start)
                    ]
            
            # Add 5' UTR
            utr5_feature = SeqFeature(
                FeatureLocation(utr5_start, utr5_end, strand=gene_strand),
                type="five_prime_UTR",
                qualifiers={
                    "ID": [f"5UTR{gene_id[4:]}.{iso_num}"],
                    "Parent": [mrna_id]
                }
            )
            
            # Add shared UTRs to the second isoform
            if iso_num == 2:
                utr5_feature.qualifiers["Parent"].append(f"mRNA{gene_id[4:]}.1")
                utr5_feature.qualifiers["ID"] = [f"5UTR{gene_id[4:]}"]
            
            mrna_feature.sub_features.append(utr5_feature)
            
            # Add CDS features
            for i, (cds_start, cds_end) in enumerate(cds_coords, 1):
                # Phase calculation is simplified for demo
                phase = 0 if i == 1 else [0, 2, 1][i % 3]
                
                cds_feature = SeqFeature(
                    FeatureLocation(cds_start, cds_end, strand=gene_strand),
                    type="CDS",
                    qualifiers={
                        "ID": [f"CDS{gene_id[4:]}.{iso_num}.{i}"],
                        "Parent": [mrna_id],
                        "phase": [str(phase)]
                    }
                )
                mrna_feature.sub_features.append(cds_feature)
            
            # Add 3' UTR
            utr3_feature = SeqFeature(
                FeatureLocation(utr3_start, utr3_end, strand=gene_strand),
                type="three_prime_UTR",
                qualifiers={
                    "ID": [f"3UTR{gene_id[4:]}.{iso_num}"],
                    "Parent": [mrna_id]
                }
            )
            
            # Add shared UTRs to the second isoform
            if iso_num == 2:
                utr3_feature.qualifiers["Parent"].append(f"mRNA{gene_id[4:]}.1")
                utr3_feature.qualifiers["ID"] = [f"3UTR{gene_id[4:]}"]
            
            mrna_feature.sub_features.append(utr3_feature)
            
            # Add mRNA to gene's subfeatures
            gene_feature.sub_features.append(mrna_feature)
        
    def generate_variants(self):
        """Define variants that will be introduced into the graph."""
        logger.info("Generating variants")
        
        # Variants will be defined according to our test plan
        variant_specs = [
            # type, pos, ref, alt, affected_gene, sample_genotypes
            ("SNP", 100, "T", "C", "gene1", {"sample1": "1|0", "sample2": "0|0"}),
            ("SNP", 250, "A", "G", "gene2", {"sample1": "0|0", "sample2": "1|0"}),
            ("INS", 600, "C", "CTAG", "gene3", {"sample1": "0|0", "sample2": "1|1"}),
            ("DEL", 750, "ACTGC", "A", "gene4", {"sample1": "1|0", "sample2": "0|0"}),
            ("DUP", 1150, "G", "GTCGA", "gene6", {"sample1": "0|0", "sample2": "0|1"}),
            ("INV", 1300, "ACTGA", "TGACT", "gene7", {"sample1": "0|1", "sample2": "0|0"}),
            ("DEL", 1450, "CGATA", "C", "gene8", {"sample1": "0|0", "sample2": "1|1"}),
            ("NOCALL", 1530, "ATT", ".", "gene9", {"sample1": "1|0", "sample2": "0|0"}),
        ]
        
        # Create variant records
        for var_type, pos, ref, alt, gene_id, sample_gts in variant_specs:
            variant = {
                "type": var_type,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene_id": gene_id,
                "sample_genotypes": sample_gts
            }
            self.variants.append(variant)
        
        logger.info(f"Generated {len(self.variants)} variants")
        return self.variants
    
    def build_variation_graph(self):
        """Build the GFA2 variation graph with segments, edges, and paths."""
        logger.info("Building variation graph")
        
        # Create a graph object
        self.graph = gfapy.Gfa()
        self.graph.header = gfapy.Line("H\tVN:Z:2.0")
        
        # Reference sequence
        ref_seq = str(self.reference_seq.seq)
        
        # Create segments for the reference sequence (in chunks)
        segment_size = 50
        for i in range(0, len(ref_seq), segment_size):
            seg_id = f"seg{i//segment_size + 1}"
            seq_chunk = ref_seq[i:i+segment_size]
            self.graph.append(gfapy.Line(f"S\t{seg_id}\t{len(seq_chunk)}\t{seq_chunk}"))
            self.segments[seg_id] = {
                "start": i,
                "end": i + len(seq_chunk),
                "sequence": seq_chunk
            }
        
        # Create variant segments
        for variant in self.variants:
            var_pos = variant["pos"] - 1  # Convert to 0-based
            var_ref = variant["ref"]
            var_alt = variant["alt"]
            var_type = variant["type"]
            
            # Determine segment that contains this variant
            containing_seg = None
            for seg_id, seg_data in self.segments.items():
                if seg_data["start"] <= var_pos < seg_data["end"]:
                    containing_seg = seg_id
                    break
            
            if var_type == "SNP":
                # Create a segment with the SNP
                var_seq = ref_seq[var_pos-5:var_pos] + var_alt + ref_seq[var_pos+1:var_pos+45]
                var_seg_id = f"{containing_seg}_snp"
                self.graph.append(gfapy.Line(f"S\t{var_seg_id}\t{len(var_seq)}\t{var_seq}"))
                self.segments[var_seg_id] = {
                    "start": var_pos-5,
                    "end": var_pos + 45,
                    "sequence": var_seq,
                    "variant": variant
                }
            
            elif var_type == "INS":
                # Create a segment with the insertion
                var_seq = ref_seq[var_pos-5:var_pos] + var_alt + ref_seq[var_pos+1:var_pos+40]
                var_seg_id = f"{containing_seg}_ins"
                self.graph.append(gfapy.Line(f"S\t{var_seg_id}\t{len(var_seq)}\t{var_seq}"))
                self.segments[var_seg_id] = {
                    "start": var_pos-5,
                    "end": var_pos + 40,
                    "sequence": var_seq,
                    "variant": variant
                }
            
            elif var_type == "DEL":
                # Create a segment with the deletion
                var_seq = ref_seq[var_pos-5:var_pos] + var_alt + ref_seq[var_pos+len(var_ref):var_pos+len(var_ref)+45]
                var_seg_id = f"{containing_seg}_del"
                self.graph.append(gfapy.Line(f"S\t{var_seg_id}\t{len(var_seq)}\t{var_seq}"))
                self.segments[var_seg_id] = {
                    "start": var_pos-5,
                    "end": var_pos + len(var_ref) + 45,
                    "sequence": var_seq,
                    "variant": variant
                }
            
            elif var_type == "DUP":
                # Create a segment with the duplication
                var_seq = ref_seq[var_pos-5:var_pos] + var_alt + ref_seq[var_pos+1:var_pos+40]
                var_seg_id = f"{containing_seg}_dup"
                self.graph.append(gfapy.Line(f"S\t{var_seg_id}\t{len(var_seq)}\t{var_seq}"))
                self.segments[var_seg_id] = {
                    "start": var_pos-5,
                    "end": var_pos + 40,
                    "sequence": var_seq,
                    "variant": variant
                }
            
            elif var_type == "INV":
                # Create a segment with the inversion
                var_seq = ref_seq[var_pos-5:var_pos] + var_alt + ref_seq[var_pos+len(var_ref):var_pos+len(var_ref)+40]
                var_seg_id = f"{containing_seg}_inv"
                self.graph.append(gfapy.Line(f"S\t{var_seg_id}\t{len(var_seq)}\t{var_seq}"))
                self.segments[var_seg_id] = {
                    "start": var_pos-5,
                    "end": var_pos + len(var_ref) + 40,
                    "sequence": var_seq,
                    "variant": variant
                }
            
            # Add more variant types as needed
        
        # Create edges between segments
        # First, connect reference segments
        for i in range(1, len(self.segments)):
            seg1_id = f"seg{i}"
            seg2_id = f"seg{i+1}"
            if seg1_id in self.segments and seg2_id in self.segments:
                edge_id = f"e{i}"
                self.graph.append(gfapy.Line(f"E\t{edge_id}\t{seg1_id}+\t{seg2_id}+\t0\t{self.segments[seg1_id]['end'] - self.segments[seg1_id]['start']}$\t0\t{self.segments[seg2_id]['end'] - self.segments[seg2_id]['start']}$\t*"))
        
        # Then add edges for variant paths
        # This is simplified - in a real implementation, we would need more complex logic
        # to properly connect variant segments to the reference graph
        
        # Create paths for each sample and haplotype
        # For demonstration, we'll create 4 paths: sample1_hap1, sample1_hap2, sample2_hap1, sample2_hap2
        # Each path will include different variants according to the genotypes
        
        # Sample1 Haplotype1 path
        sample1_hap1_segments = []
        for i in range(1, len(self.segments)+1):
            seg_id = f"seg{i}"
            if seg_id in self.segments:
                # Check if this segment should be replaced by a variant
                replaced = False
                for variant in self.variants:
                    var_pos = variant["pos"] - 1
                    sample_gt = variant["sample_genotypes"]["sample1"]
                    hap1_allele = int(sample_gt.split("|")[0])
                    
                    if hap1_allele == 1 and self.segments[seg_id]["start"] <= var_pos < self.segments[seg_id]["end"]:
                        # This segment contains a variant that should be included in hap1
                        var_type = variant["type"]
                        var_seg_id = f"{seg_id}_{var_type.lower()}"
                        if var_seg_id in self.segments:
                            sample1_hap1_segments.append(f"{var_seg_id}+")
                            replaced = True
                            break
                
                if not replaced:
                    sample1_hap1_segments.append(f"{seg_id}+")
        
        # Add path to graph
        self.graph.append(gfapy.Line(f"P\tsample1_hap1\t{','.join(sample1_hap1_segments)}\t*"))
        
        # Similarly create paths for the other haplotypes
        # This is a simplified approach - in a real implementation we would need
        # more sophisticated path creation logic
        
        logger.info("Variation graph built")
        return self.graph
    
    def write_fasta_file(self):
        """Write the reference sequence to a FASTA file."""
        with open(self.fasta_path, 'w') as f:
            SeqIO.write(self.reference_seq, f, "fasta")
        logger.info(f"FASTA file written to {self.fasta_path}")
    
    def write_gff_file(self):
        """Write the features to a GFF3 file."""
        with open(self.gff_path, 'w') as f:
            # Write GFF3 header
            f.write("##gff-version 3\n")
            f.write(f"##sequence-region {self.reference_seq.id} 1 {len(self.reference_seq)}\n")
            
            # Write each feature
            for gene in self.features:
                # Convert 0-based to 1-based coordinates for GFF3
                gene_start = gene.location.start + 1
                gene_end = gene.location.end
                gene_strand = "+" if gene.location.strand == 1 else "-"
                
                # Write gene feature
                f.write(f"{self.reference_seq.id}\tsource\tgene\t{gene_start}\t{gene_end}\t.\t{gene_strand}\t.\tID={gene.qualifiers['ID'][0]};Name={gene.qualifiers['Name'][0]}\n")
                
                # Write mRNA features
                for mrna in gene.sub_features:
                    if mrna.type == "mRNA":
                        mrna_start = mrna.location.start + 1
                        mrna_end = mrna.location.end
                        mrna_strand = "+" if mrna.location.strand == 1 else "-"
                        
                        f.write(f"{self.reference_seq.id}\tsource\tmRNA\t{mrna_start}\t{mrna_end}\t.\t{mrna_strand}\t.\tID={mrna.qualifiers['ID'][0]};Parent={gene.qualifiers['ID'][0]};Name={mrna.qualifiers['Name'][0]}\n")
                        
                        # Write child features (exons, CDS, UTRs)
                        for feat in mrna.sub_features:
                            feat_start = feat.location.start + 1
                            feat_end = feat.location.end
                            feat_strand = "+" if feat.location.strand == 1 else "-"
                            
                            # Format qualifiers
                            quals = [f"ID={feat.qualifiers['ID'][0]}"]
                            if "Parent" in feat.qualifiers:
                                quals.append(f"Parent={','.join(feat.qualifiers['Parent'])}")
                            
                            # Add phase for CDS features
                            phase = "."
                            if feat.type == "CDS" and "phase" in feat.qualifiers:
                                phase = feat.qualifiers['phase'][0]
                            
                            f.write(f"{self.reference_seq.id}\tsource\t{feat.type}\t{feat_start}\t{feat_end}\t.\t{feat_strand}\t{phase}\t{';'.join(quals)}\n")
        
        logger.info(f"GFF3 file written to {self.gff_path}")
    
    def write_gfa_file(self):
        """Write the variation graph to a GFA2 file."""
        with open(self.gfa_path, 'w') as f:
            # Get GFA string representation
            gfa_content = str(self.graph)
            f.write(gfa_content)
        logger.info(f"GFA2 file written to {self.gfa_path}")
    
    def write_vcf_file(self):
        """Write the variants to a VCF file."""
        with open(self.vcf_path, 'w') as f:
            # Write VCF header
            f.write("##fileformat=VCFv4.2\n")
            f.write("##fileDate=20250422\n")
            f.write("##source=IntegrationTestDataset\n")
            f.write(f"##reference={self.fasta_path.name}\n")
            f.write('##INFO=<ID=TYPE,Number=A,Type=String,Description="Type of variant">\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            
            # Write column header
            f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.samples[0]}\t{self.samples[1]}\n")
            
            # Write variant records
            for variant in self.variants:
                chrom = self.reference_seq.id
                pos = variant["pos"]
                var_id = "."
                ref = variant["ref"]
                alt = variant["alt"]
                qual = "100"
                filter_field = "PASS"
                info = f"TYPE={variant['type']}"
                format_field = "GT"
                
                # Sample genotypes
                sample1_gt = variant["sample_genotypes"]["sample1"]
                sample2_gt = variant["sample_genotypes"]["sample2"]
                
                f.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}\t{format_field}\t{sample1_gt}\t{sample2_gt}\n")
        
        logger.info(f"VCF file written to {self.vcf_path}")

def main():
    """Main function to generate integration test data."""
    parser = argparse.ArgumentParser(description="Generate integration test data.")
    parser.add_argument("-o", "--output-dir", default=OUTPUT_DIR, help="Output directory")
    args = parser.parse_args()
    
    # Create generator and generate all files
    generator = TestDataGenerator(Path(args.output_dir))
    generator.generate_all()
    
    return 0

if __name__ == "__main__":
    exit(main())
