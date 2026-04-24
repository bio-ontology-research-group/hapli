import pytest
import pysam
import gffutils
from pathlib import Path

@pytest.fixture
def data_dir(tmp_path):
    return tmp_path / "data"

@pytest.fixture
def dummy_fasta(tmp_path):
    """Creates a dummy reference FASTA file."""
    p = tmp_path / "ref.fa"
    # Create a complex sequence to avoid poly-A/T issues
    # chr1: 100bp
    # 0-49:  ACGT...
    # 50-99: TGCA...
    seq = "ACGT" * 12 + "AC" + "TGCA" * 12 + "TG"
    with open(p, "w") as f:
        f.write(f">chr1\n{seq}\n")
    
    # Index it (requires pysam/samtools)
    pysam.faidx(str(p))
    return p

@pytest.fixture
def dummy_gff(tmp_path):
    """Creates a dummy GFF3 file."""
    p = tmp_path / "features.gff3"
    content = """##gff-version 3
chr1\ttest\tgene\t10\t90\t.\t+\t.\tID=gene1;Name=TestGene;gene_name=TestGene
chr1\ttest\tmRNA\t10\t90\t.\t+\t.\tID=mrna1;Parent=gene1;Name=TestTranscript
chr1\ttest\texon\t10\t40\t.\t+\t.\tID=exon1;Parent=mrna1
chr1\ttest\texon\t60\t90\t.\t+\t.\tID=exon2;Parent=mrna1
"""
    with open(p, "w") as f:
        f.write(content)
    return p

@pytest.fixture
def dummy_vcf(tmp_path):
    """Creates a dummy VCF file."""
    p = tmp_path / "variants.vcf"
    content = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t20\t.\tT\tC\t.\tPASS\t.\tGT\t0|1
chr1\t70\t.\tA\tG\t.\tPASS\t.\tGT\t1|1
"""
    with open(p, "w") as f:
        f.write(content)
    
    # Compress and index (requires tabix/bgzip usually, but pysam can write directly)
    # Using pysam to write a proper VCF usually safer but text is easier for simple tests
    # We need to bgzip and tabix it for pysam.VariantFile to work well with fetch if we use it that way?
    # Actually pysam.VariantFile can read uncompressed VCFs too, but fetch might require index.
    
    # Let's use pysam to create it properly if possible, or just gzip it.
    # pysam.tabix_index can index vcf.gz
    
    # Write text VCF
    with open(p, "w") as f:
        f.write(content)
    
    # Compress
    p_gz = tmp_path / "variants.vcf.gz"
    pysam.tabix_compress(str(p), str(p_gz), force=True)
    pysam.tabix_index(str(p_gz), preset="vcf", force=True)
    
    return p_gz
