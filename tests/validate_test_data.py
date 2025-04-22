#!/usr/bin/env python3
"""
Validation script for integration test data.
Checks that all files exist, are formatted correctly, and contain required elements.
"""

import os
import sys
import re
import logging
from typing import Dict, List, Set
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Path to the integration test data directory
TEST_DATA_DIR = Path("data/integration_test")
REQUIRED_FILES = [
    "reference.fasta",
    "annotations.gff3",
    "variation_graph.gfa",
    "variants.vcf",
    "README.md"
]

def check_files_exist() -> bool:
    """Check if all required files exist in the test data directory."""
    logger.info("Checking for required files...")
    
    all_exist = True
    for filename in REQUIRED_FILES:
        file_path = TEST_DATA_DIR / filename
        if file_path.exists():
            logger.info(f"✓ Found {filename}")
        else:
            logger.error(f"✗ Missing {filename}")
            all_exist = False
    
    return all_exist

def validate_fasta_file() -> bool:
    """Validate the FASTA file format and contents."""
    logger.info("Validating FASTA file...")
    
    fasta_path = TEST_DATA_DIR / "reference.fasta"
    if not fasta_path.exists():
        logger.error("FASTA file not found")
        return False
    
    valid = True
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
        
        # Check for FASTA header
        if not lines or not lines[0].startswith('>'):
            logger.error("FASTA file missing header line (should start with '>')")
            valid = False
        
        # Check sequence lines
        seq_lines = [line.strip() for line in lines[1:] if line.strip()]
        for i, line in enumerate(seq_lines):
            if not re.match(r'^[ACGTNacgtn]+$', line):
                logger.error(f"Invalid FASTA sequence at line {i+2}: contains non-ACGTN characters")
                valid = False
    
    if valid:
        logger.info("✓ FASTA file is valid")
    return valid

def validate_gff3_file() -> bool:
    """Validate the GFF3 file format and contents."""
    logger.info("Validating GFF3 file...")
    
    gff_path = TEST_DATA_DIR / "annotations.gff3"
    if not gff_path.exists():
        logger.error("GFF3 file not found")
        return False
    
    valid = True
    with open(gff_path, 'r') as f:
        lines = f.readlines()
        
        # Check for GFF3 header
        if not lines or not any(line.strip() == "##gff-version 3" for line in lines):
            logger.error("GFF3 file missing header line (##gff-version 3)")
            valid = False
        
        # Check for required feature types
        feature_types = set()
        feature_ids = set()
        parent_ids = set()
        
        for line in lines:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) != 9:
                logger.error(f"Invalid GFF3 line: {line.strip()} (should have 9 columns)")
                valid = False
                continue
            
            feature_types.add(parts[2])
            
            # Extract IDs and Parent references
            attrs = parts[8].split(';')
            id_match = next((a for a in attrs if a.startswith('ID=')), None)
            parent_match = next((a for a in attrs if a.startswith('Parent=')), None)
            
            if id_match:
                feature_id = id_match.split('=')[1]
                feature_ids.add(feature_id)
            
            if parent_match:
                parent_id = parent_match.split('=')[1]
                # Handle multiple parents
                for pid in parent_id.split(','):
                    parent_ids.add(pid)
        
        # Check for required feature types
        required_types = {'gene', 'mRNA', 'exon', 'CDS'}
        for req_type in required_types:
            if req_type not in feature_types:
                logger.error(f"Missing required feature type: {req_type}")
                valid = False
        
        # Check for UTR features
        if not any(ft.endswith('UTR') for ft in feature_types):
            logger.error("Missing UTR feature types (5' or 3' UTR)")
            valid = False
        
        # Check parent-child relationships
        for parent in parent_ids:
            if parent not in feature_ids:
                logger.error(f"Referenced parent ID '{parent}' does not exist in GFF3")
                valid = False
    
    if valid:
        logger.info("✓ GFF3 file is valid")
    return valid

def validate_gfa_file() -> bool:
    """Validate the GFA file format and contents."""
    logger.info("Validating GFA file...")
    
    gfa_path = TEST_DATA_DIR / "variation_graph.gfa"
    if not gfa_path.exists():
        logger.error("GFA file not found")
        return False
    
    valid = True
    segment_ids = set()
    path_ids = set()
    
    with open(gfa_path, 'r') as f:
        lines = f.readlines()
        
        # Check for GFA2 header
        if not lines or not any(line.strip().startswith('H\tVN:Z:2.0') for line in lines):
            logger.error("GFA file missing GFA2 header (H\tVN:Z:2.0)")
            valid = False
        
        # Check segments, edges, and paths
        has_segments = False
        has_edges = False
        has_paths = False
        
        for line in lines:
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            if not parts:
                continue
                
            record_type = parts[0]
            
            if record_type == 'S':  # Segment
                has_segments = True
                if len(parts) < 3:
                    logger.error(f"Invalid segment line: {line.strip()}")
                    valid = False
                else:
                    segment_ids.add(parts[1])
            
            elif record_type == 'E':  # Edge
                has_edges = True
                if len(parts) < 5:
                    logger.error(f"Invalid edge line: {line.strip()}")
                    valid = False
            
            elif record_type == 'P':  # Path
                has_paths = True
                if len(parts) < 3:
                    logger.error(f"Invalid path line: {line.strip()}")
                    valid = False
                else:
                    path_ids.add(parts[1])
                    # Validate path segments
                    path_segments = parts[2].split(',')
                    for path_segment in path_segments:
                        # Remove orientation (+/-) from segment ID
                        segment_id = path_segment[:-1]
                        if segment_id not in segment_ids:
                            logger.error(f"Path references undefined segment: {segment_id}")
                            valid = False
        
        # Check for required elements
        if not has_segments:
            logger.error("GFA file has no segments (S lines)")
            valid = False
        if not has_edges:
            logger.error("GFA file has no edges (E lines)")
            valid = False
        if not has_paths:
            logger.error("GFA file has no paths (P lines)")
            valid = False
        
        # Check for required paths (sample haplotypes)
        required_paths = {'sample1_hap1', 'sample1_hap2', 'sample2_hap1', 'sample2_hap2'}
        for req_path in required_paths:
            if req_path not in path_ids:
                logger.error(f"Missing required path: {req_path}")
                valid = False
    
    if valid:
        logger.info("✓ GFA file is valid")
    return valid

def validate_vcf_file() -> bool:
    """Validate the VCF file format and contents."""
    logger.info("Validating VCF file...")
    
    vcf_path = TEST_DATA_DIR / "variants.vcf"
    if not vcf_path.exists():
        logger.error("VCF file not found")
        return False
    
    valid = True
    with open(vcf_path, 'r') as f:
        lines = f.readlines()
        
        # Check for VCF header
        if not lines or not any(line.strip() == "##fileformat=VCFv4.2" for line in lines):
            logger.error("VCF file missing header line (##fileformat=VCFv4.2)")
            valid = False
        
        # Check for column header
        column_header = None
        for line in lines:
            if line.startswith('#CHROM'):
                column_header = line
                break
        
        if not column_header:
            logger.error("VCF file missing column header line (#CHROM POS ID REF...)")
            valid = False
        else:
            columns = column_header.strip().split('\t')
            required_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
            for col in required_columns:
                if col not in columns:
                    logger.error(f"VCF column header missing required column: {col}")
                    valid = False
            
            # Check for sample columns
            if len(columns) <= 9:
                logger.error("VCF file missing sample columns")
                valid = False
            elif 'sample1' not in columns or 'sample2' not in columns:
                logger.error("VCF file missing required sample columns (sample1, sample2)")
                valid = False
        
        # Check variant lines
        variant_lines = [line for line in lines if not line.startswith('#')]
        if not variant_lines:
            logger.error("VCF file has no variant lines")
            valid = False
        
        # Check for different variant types
        variant_types = set()
        for line in variant_lines:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                logger.error(f"Invalid VCF line: {line.strip()}")
                valid = False
                continue
            
            ref = parts[3]
            alt = parts[4]
            
            if len(ref) < len(alt) and alt.startswith(ref):
                variant_types.add('INS')
            elif len(ref) > len(alt) and ref.startswith(alt):
                variant_types.add('DEL')
            elif len(ref) == 1 and len(alt) == 1:
                variant_types.add('SNP')
            elif alt == '.':
                variant_types.add('NOCALL')
            else:
                variant_types.add('COMPLEX')
        
        # Check for diverse variant types
        required_types = {'SNP', 'INS', 'DEL'}
        for req_type in required_types:
            if req_type not in variant_types:
                logger.error(f"Missing required variant type: {req_type}")
                valid = False
    
    if valid:
        logger.info("✓ VCF file is valid")
    return valid

def cross_validate_files() -> bool:
    """Cross-validate relationships between the different files."""
    logger.info("Cross-validating files...")
    
    valid = True
    
    # Check for consistency in sample IDs between VCF and GFA
    vcf_path = TEST_DATA_DIR / "variants.vcf"
    gfa_path = TEST_DATA_DIR / "variation_graph.gfa"
    
    # Extract sample IDs from VCF
    vcf_samples = set()
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                vcf_samples = set(columns[9:])  # Sample columns start at index 9
                break
    
    # Extract path prefixes from GFA
    gfa_samples = set()
    with open(gfa_path, 'r') as f:
        for line in f:
            if line.startswith('P'):
                parts = line.strip().split('\t')
                path_id = parts[1]
                if '_' in path_id:
                    sample_id = path_id.split('_')[0]
                    gfa_samples.add(sample_id)
    
    # Compare sample IDs
    if vcf_samples != gfa_samples and len(vcf_samples) >= 2 and len(gfa_samples) >= 2:
        logger.warning(f"Sample IDs in VCF ({vcf_samples}) don't match samples in GFA ({gfa_samples})")
    
    # Verify variants in GFA and VCF
    gfa_variants = set()
    with open(gfa_path, 'r') as f:
        for line in f:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                seg_id = parts[1]
                if '_' in seg_id and not seg_id.startswith('sample'):
                    var_type = seg_id.split('_')[1].lower()
                    gfa_variants.add(var_type)
    
    vcf_variants = set()
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                info = parts[7]
                if "TYPE=" in info:
                    var_type = info.split('TYPE=')[1].split(';')[0].lower()
                    vcf_variants.add(var_type)
    
    # Check if variants in GFA match those in VCF
    missing_variants = vcf_variants - gfa_variants
    if missing_variants:
        logger.warning(f"Variants in VCF not found in GFA: {missing_variants}")
        valid = False
    
    return valid

def main():
    """Run validation checks on the integration test dataset."""
    logger.info("Starting validation of integration test dataset...")
    
    if not TEST_DATA_DIR.exists():
        logger.error(f"Test data directory not found: {TEST_DATA_DIR}")
        return False
    
    # Run all validation checks
    files_exist = check_files_exist()
    if not files_exist:
        logger.error("Some required files are missing - aborting validation")
        return False
    
    fasta_valid = validate_fasta_file()
    gff_valid = validate_gff3_file()
    gfa_valid = validate_gfa_file()
    vcf_valid = validate_vcf_file()
    cross_valid = cross_validate_files()
    
    all_valid = fasta_valid and gff_valid and gfa_valid and vcf_valid and cross_valid
    
    if all_valid:
        logger.info("✓ All validation checks passed! The integration test dataset is valid.")
        return True
    else:
        logger.error("✗ Some validation checks failed. See the logs above for details.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
