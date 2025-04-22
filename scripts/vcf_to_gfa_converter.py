#!/usr/bin/env python3
"""
Script to convert VCF files to GFA format.
Allows converting multiple VCF files either jointly or separately.
"""

import argparse
import logging
import os
import sys
import subprocess
from pathlib import Path
import tempfile

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def check_prerequisites():
    """
    Check if required external tools are installed.
    """
    required_tools = ['vg']
    missing_tools = []
    
    for tool in required_tools:
        try:
            subprocess.run([tool, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logger.info(f"Found {tool}")
        except FileNotFoundError:
            missing_tools.append(tool)
    
    if missing_tools:
        logger.error(f"Missing required tools: {', '.join(missing_tools)}")
        logger.error("Please install these tools before proceeding.")
        return False
    
    return True

def convert_vcf_to_gfa_vg(vcf_paths, reference_path, output_path, region=None, memory_gb=16):
    """
    Convert VCF to GFA using the vg toolkit.
    
    Args:
        vcf_paths: List of VCF file paths
        reference_path: Path to reference FASTA
        output_path: Output GFA file path
        region: Optional region to restrict conversion (e.g., "chr1:1000-2000")
        memory_gb: Memory limit for vg in gigabytes
    """
    try:
        # Create temporary directory for intermediate files
        with tempfile.TemporaryDirectory() as temp_dir:
            # Construct vg command for multiple VCFs
            vcf_args = []
            for vcf_path in vcf_paths:
                vcf_args.extend(['-v', str(vcf_path)])
            
            region_args = []
            if region:
                region_args = ['-R', region]
            
            # Step 1: Construct the graph using vg construct
            vg_path = os.path.join(temp_dir, "graph.vg")
            construct_cmd = [
                'vg', 'construct',
                '-r', str(reference_path),
                *vcf_args,
                '-m', str(memory_gb),
                *region_args,
                '-t', '8',  # Threads
                '-o', vg_path
            ]
            
            logger.info(f"Running: {' '.join(construct_cmd)}")
            subprocess.run(construct_cmd, check=True)
            
            # Step 2: Convert to GFA
            gfa_cmd = [
                'vg', 'view',
                vg_path,
                '-F', '-g',  # Output in GFA format
                '-o', str(output_path)
            ]
            
            logger.info(f"Running: {' '.join(gfa_cmd)}")
            subprocess.run(gfa_cmd, check=True)
            
            logger.info(f"Successfully created GFA file: {output_path}")
            return True
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running vg: {e}")
        return False
    except Exception as e:
        logger.error(f"Error converting VCF to GFA: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Convert VCF files to GFA format")
    parser.add_argument('--vcf', '-v', required=True, action='append',
                        help='VCF file paths (can specify multiple times)')
    parser.add_argument('--reference', '-r', required=True,
                        help='Reference genome FASTA file')
    parser.add_argument('--output', '-o', required=True,
                        help='Output GFA file or directory')
    parser.add_argument('--mode', choices=['joint', 'separate'], default='joint',
                        help='Convert VCFs jointly or separately')
    parser.add_argument('--region', type=str, default=None,
                        help='Restrict to region (e.g., "chr1:1000-2000")')
    parser.add_argument('--memory', type=int, default=16,
                        help='Memory limit in GB for vg')
    
    args = parser.parse_args()
    
    # Check if required tools are installed
    if not check_prerequisites():
        sys.exit(1)
    
    # Validate input files
    for vcf_path in args.vcf:
        if not os.path.exists(vcf_path):
            logger.error(f"VCF file does not exist: {vcf_path}")
            sys.exit(1)
    
    if not os.path.exists(args.reference):
        logger.error(f"Reference file does not exist: {args.reference}")
        sys.exit(1)
    
    # Process based on mode
    if args.mode == 'joint':
        # Convert all VCFs together
        output_path = Path(args.output)
        if output_path.suffix != '.gfa':
            output_path = output_path / "joint_output.gfa"
            
        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        success = convert_vcf_to_gfa_vg(
            args.vcf,
            args.reference,
            output_path,
            region=args.region,
            memory_gb=args.memory
        )
        
        if not success:
            logger.error("Joint conversion failed")
            sys.exit(1)
            
    else:  # separate mode
        output_dir = Path(args.output)
        if output_dir.suffix == '.gfa':
            # If output is specified as a file, use its parent directory
            output_dir = output_dir.parent
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Process each VCF separately
        for vcf_path in args.vcf:
            vcf_name = Path(vcf_path).stem
            if vcf_name.endswith('.vcf'):
                vcf_name = vcf_name[:-4]
                
            output_path = output_dir / f"{vcf_name}.gfa"
            
            success = convert_vcf_to_gfa_vg(
                [vcf_path],
                args.reference,
                output_path,
                region=args.region,
                memory_gb=args.memory
            )
            
            if not success:
                logger.error(f"Conversion failed for {vcf_path}")
                continue
    
    logger.info("Conversion completed successfully")

if __name__ == "__main__":
    main()
