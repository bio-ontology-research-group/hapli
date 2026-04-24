import random
from dataclasses import dataclass
from typing import List, Tuple
from pathlib import Path

@dataclass
class Exon:
    start: int
    end: int

@dataclass
class Transcript:
    id: str
    parent_id: str
    start: int
    end: int
    exons: List[Exon]
    strand: str = "+"

@dataclass
class Gene:
    id: str
    name: str
    chrom: str
    start: int
    end: int
    transcripts: List[Transcript]
    strand: str = "+"
    biotype: str = "protein_coding"

class GFFGenerator:
    def __init__(self, output_path: Path):
        self.output_path = output_path
        self.genes = []

    def generate_random_gene(self, chrom: str, start_range: Tuple[int, int], gene_id_prefix: str = "GENE", gene_idx: int = 1) -> Gene:
        """
        Generates a random gene structure within the start_range.
        Note: Checks should be done by caller to ensure no overlap with other genes if desired.
        """
        start = random.randint(start_range[0], start_range[1])
        # Gene length 1kb to 10kb
        length = random.randint(1000, 10000)
        end = start + length
        strand = random.choice(["+", "-"])
        
        gene_id = f"{gene_id_prefix}{gene_idx:05d}"
        gene_name = f"Gene_{gene_idx}"
        
        # Generate 1-3 transcripts
        num_transcripts = random.randint(1, 3)
        transcripts = []
        
        for t_i in range(num_transcripts):
            t_id = f"{gene_id}.t{t_i+1}"
            
            # Generate exons (3 to 10)
            num_exons = random.randint(3, 10)
            exons = []
            
            # Create random exon structure within gene limits
            # Simple strategy: divide gene into (num_exons * 2 + 1) parts (exon, intron, exon...)
            # This is naive but sufficient for topology testing.
            
            current_pos = start + random.randint(0, 50) # Slight offset from gene start
            remaining_len = end - current_pos - 50
            
            # Average step
            step = remaining_len // (num_exons * 2)
            if step < 50: step = 50 # min size
            
            for e_i in range(num_exons):
                if current_pos >= end: break
                
                # Exon len
                exon_len = random.randint(50, max(51, step))
                exon_end = current_pos + exon_len
                if exon_end > end: exon_end = end
                
                exons.append(Exon(current_pos, exon_end))
                
                # Intron gap
                intron_len = random.randint(50, max(51, step))
                current_pos = exon_end + intron_len
            
            if not exons:
                continue

            t_start = exons[0].start
            t_end = exons[-1].end
            
            transcripts.append(Transcript(t_id, gene_id, t_start, t_end, exons, strand))
            
        return Gene(gene_id, gene_name, chrom, start, end, transcripts, strand)

    def write(self, genes: List[Gene]):
        with open(self.output_path, 'w') as f:
            f.write("##gff-version 3\n")
            for gene in genes:
                # Write Gene
                f.write(f"{gene.chrom}\tSim\tgene\t{gene.start}\t{gene.end}\t.\t{gene.strand}\t.\tID={gene.id};Name={gene.name};biotype={gene.biotype}\n")
                
                for t in gene.transcripts:
                    # Write mRNA
                    f.write(f"{gene.chrom}\tSim\tmRNA\t{t.start}\t{t.end}\t.\t{t.strand}\t.\tID={t.id};Parent={gene.id}\n")
                    
                    # Write Exons
                    for i, exon in enumerate(t.exons):
                        # Calculate phase for CDS if we wanted to be fancy, but simple exon/CDS suffices
                        f.write(f"{gene.chrom}\tSim\texon\t{exon.start}\t{exon.end}\t.\t{t.strand}\t.\tID={t.id}.exon{i+1};Parent={t.id}\n")
                        # Assume whole exon is CDS for simplicity in this synthetic data
                        f.write(f"{gene.chrom}\tSim\tCDS\t{exon.start}\t{exon.end}\t.\t{t.strand}\t0\tID={t.id}.cds{i+1};Parent={t.id}\n")
