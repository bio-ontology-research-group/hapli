#!/usr/bin/env python3
"""
Generate genome features in GFF3 format from a reference FASTA file.

This script creates hierarchical genome features including genes, exons, introns,
and splice sites, ensuring no overlaps between features.
"""

import argparse
import logging
import random
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


@dataclass
class Feature:
    """Represents a genome feature."""
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: Dict[str, str]
    
    def to_gff3_line(self) -> str:
        """Convert feature to GFF3 format line."""
        attr_str = ';'.join([f"{k}={v}" for k, v in self.attributes.items()])
        return f"{self.seqid}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{attr_str}"


class GenomeFeatureGenerator:
    """Generate genome features for a reference genome."""
    
    def __init__(self, seed: int = None):
        """Initialize the generator with optional random seed."""
        if seed is not None:
            random.seed(seed)
        self.features = []
        self.used_regions = {}  # chromosome -> list of (start, end) tuples
        
    def load_reference_genome(self, fasta_path: Path) -> Dict[str, SeqRecord]:
        """Load reference genome from FASTA file."""
        logging.info(f"Loading reference genome from {fasta_path}")
        genome = {}
        with open(fasta_path, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                genome[record.id] = record
                self.used_regions[record.id] = []
        logging.info(f"Loaded {len(genome)} chromosomes")
        return genome
    
    def is_region_available(self, chromosome: str, start: int, end: int, buffer: int = 100) -> bool:
        """Check if a region is available (not overlapping with existing features)."""
        for used_start, used_end in self.used_regions[chromosome]:
            # Check for overlap with buffer
            if not (end + buffer < used_start or start - buffer > used_end):
                return False
        return True
    
    def reserve_region(self, chromosome: str, start: int, end: int) -> None:
        """Reserve a region to prevent overlaps."""
        self.used_regions[chromosome].append((start, end))
        # Keep regions sorted for efficient checking
        self.used_regions[chromosome].sort()
    
    def find_available_region(self, chromosome: str, length: int, max_attempts: int = 1000) -> Tuple[int, int]:
        """Find an available region of specified length on a chromosome."""
        chr_length = len(self.genome[chromosome].seq)
        
        for _ in range(max_attempts):
            start = random.randint(1, max(1, chr_length - length))
            end = start + length - 1
            
            if self.is_region_available(chromosome, start, end):
                return start, end
        
        raise ValueError(f"Could not find available region of length {length} on {chromosome}")
    
    def generate_gene_structure(self, chromosome: str, gene_id: str, gene_start: int, gene_end: int, strand: str) -> List[Feature]:
        """Generate a complete gene structure with exons and introns."""
        features = []
        
        # Gene feature
        gene_attrs = {
            'ID': gene_id,
            'Name': f"gene_{gene_id}",
            'biotype': 'protein_coding'
        }
        gene_feature = Feature(
            seqid=chromosome,
            source='synthetic',
            type='gene',
            start=gene_start,
            end=gene_end,
            score='.',
            strand=strand,
            phase='.',
            attributes=gene_attrs
        )
        features.append(gene_feature)
        
        # mRNA feature
        mrna_id = f"{gene_id}.1"
        mrna_attrs = {
            'ID': mrna_id,
            'Parent': gene_id,
            'Name': f"transcript_{mrna_id}"
        }
        mrna_feature = Feature(
            seqid=chromosome,
            source='synthetic',
            type='mRNA',
            start=gene_start,
            end=gene_end,
            score='.',
            strand=strand,
            phase='.',
            attributes=mrna_attrs
        )
        features.append(mrna_feature)
        
        # Generate exons and introns
        gene_length = gene_end - gene_start + 1
        num_exons = random.randint(2, min(8, max(2, gene_length // 500)))  # 2-8 exons
        
        # Calculate exon and intron sizes
        total_exon_length = int(gene_length * random.uniform(0.6, 0.8))  # 60-80% exons
        avg_exon_length = total_exon_length // num_exons
        
        exon_positions = []
        current_pos = gene_start
        
        for i in range(num_exons):
            if i == num_exons - 1:
                # Last exon goes to gene end
                exon_start = current_pos
                exon_end = gene_end
            else:
                # Random exon length around average
                exon_length = max(50, int(avg_exon_length * random.uniform(0.5, 1.5)))
                exon_start = current_pos
                exon_end = min(gene_end, exon_start + exon_length - 1)
                
                # Add intron
                intron_length = random.randint(100, 1000)
                current_pos = exon_end + intron_length + 1
                
                # Ensure we don't exceed gene boundaries
                if current_pos >= gene_end:
                    exon_end = gene_end
                    current_pos = gene_end + 1
            
            exon_positions.append((exon_start, exon_end))
            if current_pos > gene_end:
                break
        
        # Create exon features
        for i, (exon_start, exon_end) in enumerate(exon_positions):
            exon_id = f"{mrna_id}.exon{i+1}"
            exon_attrs = {
                'ID': exon_id,
                'Parent': mrna_id,
                'exon_number': str(i + 1)
            }
            exon_feature = Feature(
                seqid=chromosome,
                source='synthetic',
                type='exon',
                start=exon_start,
                end=exon_end,
                score='.',
                strand=strand,
                phase='.' if i > 0 else '0',  # First exon starts at phase 0
                attributes=exon_attrs
            )
            features.append(exon_feature)
            
            # Add CDS feature (same as exon for simplicity)
            cds_id = f"{mrna_id}.cds{i+1}"
            cds_attrs = {
                'ID': cds_id,
                'Parent': mrna_id
            }
            cds_feature = Feature(
                seqid=chromosome,
                source='synthetic',
                type='CDS',
                start=exon_start,
                end=exon_end,
                score='.',
                strand=strand,
                phase='.' if i > 0 else '0',
                attributes=cds_attrs
            )
            features.append(cds_feature)
        
        # Create intron features
        for i in range(len(exon_positions) - 1):
            intron_start = exon_positions[i][1] + 1
            intron_end = exon_positions[i + 1][0] - 1
            
            if intron_start <= intron_end:
                intron_id = f"{mrna_id}.intron{i+1}"
                intron_attrs = {
                    'ID': intron_id,
                    'Parent': mrna_id,
                    'intron_number': str(i + 1)
                }
                intron_feature = Feature(
                    seqid=chromosome,
                    source='synthetic',
                    type='intron',
                    start=intron_start,
                    end=intron_end,
                    score='.',
                    strand=strand,
                    phase='.',
                    attributes=intron_attrs
                )
                features.append(intron_feature)
                
                # Add splice sites
                if strand == '+':
                    # Donor site (GT) at start of intron
                    donor_attrs = {
                        'ID': f"{intron_id}.donor",
                        'Parent': intron_id,
                        'splice_type': 'donor'
                    }
                    donor_feature = Feature(
                        seqid=chromosome,
                        source='synthetic',
                        type='splice_donor_site',
                        start=intron_start,
                        end=intron_start + 1,
                        score='.',
                        strand=strand,
                        phase='.',
                        attributes=donor_attrs
                    )
                    features.append(donor_feature)
                    
                    # Acceptor site (AG) at end of intron
                    acceptor_attrs = {
                        'ID': f"{intron_id}.acceptor",
                        'Parent': intron_id,
                        'splice_type': 'acceptor'
                    }
                    acceptor_feature = Feature(
                        seqid=chromosome,
                        source='synthetic',
                        type='splice_acceptor_site',
                        start=intron_end - 1,
                        end=intron_end,
                        score='.',
                        strand=strand,
                        phase='.',
                        attributes=acceptor_attrs
                    )
                    features.append(acceptor_feature)
                else:  # negative strand
                    # Acceptor site (CT) at start of intron (genomic coordinates)
                    acceptor_attrs = {
                        'ID': f"{intron_id}.acceptor",
                        'Parent': intron_id,
                        'splice_type': 'acceptor'
                    }
                    acceptor_feature = Feature(
                        seqid=chromosome,
                        source='synthetic',
                        type='splice_acceptor_site',
                        start=intron_start,
                        end=intron_start + 1,
                        score='.',
                        strand=strand,
                        phase='.',
                        attributes=acceptor_attrs
                    )
                    features.append(acceptor_feature)
                    
                    # Donor site (AC) at end of intron (genomic coordinates)
                    donor_attrs = {
                        'ID': f"{intron_id}.donor",
                        'Parent': intron_id,
                        'splice_type': 'donor'
                    }
                    donor_feature = Feature(
                        seqid=chromosome,
                        source='synthetic',
                        type='splice_donor_site',
                        start=intron_end - 1,
                        end=intron_end,
                        score='.',
                        strand=strand,
                        phase='.',
                        attributes=donor_attrs
                    )
                    features.append(donor_feature)
        
        return features
    
    def generate_features_for_chromosome(self, chromosome: str, num_genes: int) -> List[Feature]:
        """Generate features for a single chromosome."""
        logging.info(f"Generating {num_genes} genes for chromosome {chromosome}")
        features = []
        
        chr_length = len(self.genome[chromosome].seq)
        
        for gene_num in range(num_genes):
            try:
                # Random gene length between 1kb and 50kb
                gene_length = random.randint(1000, 50000)
                gene_start, gene_end = self.find_available_region(chromosome, gene_length)
                
                # Random strand
                strand = random.choice(['+', '-'])
                
                # Generate gene ID
                gene_id = f"gene_{chromosome}_{gene_num + 1:04d}"
                
                # Generate gene structure
                gene_features = self.generate_gene_structure(chromosome, gene_id, gene_start, gene_end, strand)
                features.extend(gene_features)
                
                # Reserve the region
                self.reserve_region(chromosome, gene_start, gene_end)
                
            except ValueError as e:
                logging.warning(f"Could not place gene {gene_num + 1} on {chromosome}: {e}")
                continue
        
        logging.info(f"Generated {len([f for f in features if f.type == 'gene'])} genes on {chromosome}")
        return features
    
    def generate_all_features(self, genome: Dict[str, SeqRecord], genes_per_mb: float = 20.0) -> List[Feature]:
        """Generate features for all chromosomes."""
        logging.info("Generating genome features")
        self.genome = genome
        all_features = []
        
        for chromosome, record in genome.items():
            chr_length = len(record.seq)
            num_genes = max(1, int((chr_length / 1_000_000) * genes_per_mb))
            
            chr_features = self.generate_features_for_chromosome(chromosome, num_genes)
            all_features.extend(chr_features)
        
        # Sort features by chromosome and position
        all_features.sort(key=lambda f: (f.seqid, f.start, f.end))
        
        logging.info(f"Generated {len(all_features)} total features")
        return all_features
    
    def write_gff3(self, features: List[Feature], output_path: Path) -> None:
        """Write features to GFF3 file."""
        logging.info(f"Writing GFF3 file to {output_path}")
        
        with open(output_path, 'w') as f:
            # Write GFF3 header
            f.write("##gff-version 3\n")
            
            # Write sequence regions
            for chromosome, record in self.genome.items():
                f.write(f"##sequence-region {chromosome} 1 {len(record.seq)}\n")
            
            # Write features
            for feature in features:
                f.write(feature.to_gff3_line() + "\n")
        
        logging.info(f"Wrote {len(features)} features to {output_path}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Generate genome features in GFF3 format from reference FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        'reference',
        type=Path,
        help='Input reference FASTA file'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=Path,
        default='genome_features.gff3',
        help='Output GFF3 file'
    )
    
    parser.add_argument(
        '--genes-per-mb',
        type=float,
        default=20.0,
        help='Number of genes per megabase'
    )
    
    parser.add_argument(
        '--seed',
        type=int,
        help='Random seed for reproducible results'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Validate inputs
    if not args.reference.exists():
        logging.error(f"Reference file not found: {args.reference}")
        sys.exit(1)
    
    try:
        # Initialize generator
        generator = GenomeFeatureGenerator(seed=args.seed)
        
        # Load reference genome
        genome = generator.load_reference_genome(args.reference)
        
        # Generate features
        features = generator.generate_all_features(genome, args.genes_per_mb)
        
        # Write GFF3 file
        generator.write_gff3(features, args.output)
        
        logging.info("Feature generation completed successfully")
        
    except Exception as e:
        logging.error(f"Error generating features: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
