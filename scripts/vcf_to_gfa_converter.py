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
            # Use a simple command like 'version' which should exist for vg
            subprocess.run([tool, 'version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            logger.info(f"Found {tool}")
        except (FileNotFoundError, subprocess.CalledProcessError):
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
                # vg construct expects region like 'chr:start-end' or just 'chr'
                # Ensure the region format is compatible if needed, or pass directly
                region_args = ['-R', region]
            
            # Step 1: Construct the graph using vg construct
            vg_path = os.path.join(temp_dir, "graph.vg")
            construct_cmd = [
                'vg', 'construct',
                '-r', str(reference_path),
                *vcf_args,
                '-m', str(memory_gb),
                *region_args,
                '-t', str(os.cpu_count() or 8),  # Use available CPUs or default to 8
                # '-a' # Add -a if you want all alleles included, otherwise it's ref+alt
            ]
            
            logger.info(f"Running: {' '.join(construct_cmd)} > {vg_path}")
            # Redirect stdout to the vg_path file
            with open(vg_path, 'w') as vg_file:
                result = subprocess.run(construct_cmd, stdout=vg_file, stderr=subprocess.PIPE, text=True, check=True)
                if result.stderr:
                    logger.info(f"vg construct stderr: {result.stderr}") # Log stderr for info

            # Check if the vg file was created and is not empty
            if not os.path.exists(vg_path) or os.path.getsize(vg_path) == 0:
                logger.error(f"vg construct failed to create a valid graph file: {vg_path}")
                # Attempt to log stderr if available from a previous failure capture
                if 'result' in locals() and hasattr(result, 'stderr'):
                     logger.error(f"vg construct stderr: {result.stderr}")
                return False

            # Step 2: Convert to GFA using vg view
            # Ensure output directory exists before writing the GFA file
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            
            gfa_cmd = [
                'vg', 'view',
                vg_path,
                '-F', '-g',  # Output in GFA format
            ]
            
            logger.info(f"Running: {' '.join(gfa_cmd)} > {output_path}")
            # Redirect stdout to the final output_path GFA file
            with open(output_path, 'w') as gfa_file:
                 result_view = subprocess.run(gfa_cmd, stdout=gfa_file, stderr=subprocess.PIPE, text=True, check=True)
                 if result_view.stderr:
                     logger.info(f"vg view stderr: {result_view.stderr}") # Log stderr for info

            # Check if the GFA file was created
            if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
                 logger.error(f"vg view failed to create a valid GFA file: {output_path}")
                 if 'result_view' in locals() and hasattr(result_view, 'stderr'):
                     logger.error(f"vg view stderr: {result_view.stderr}")
                 return False

            logger.info(f"Successfully created GFA file: {output_path}")
            return True
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running vg command: {' '.join(e.cmd)}")
        logger.error(f"Return code: {e.returncode}")
        if e.stdout:
            logger.error(f"stdout: {e.stdout}")
        if e.stderr:
            logger.error(f"stderr: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"An unexpected error occurred during VCF to GFA conversion: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def main():
    parser = argparse.ArgumentParser(description="Convert VCF files to GFA format using vg")
    parser.add_argument('--vcf', '-v', required=True, action='append',
                        help='VCF file paths (can specify multiple times). Files should be bgzipped and indexed (.tbi).')
    parser.add_argument('--reference', '-r', required=True,
                        help='Reference genome FASTA file (should be indexed with faidx).')
    parser.add_argument('--output', '-o', required=True,
                        help='Output GFA file path (for joint mode) or output directory (for separate mode).')
    parser.add_argument('--mode', choices=['joint', 'separate'], default='joint',
                        help='Mode of conversion: "joint" combines all VCFs into one GFA, "separate" creates one GFA per VCF.')
    parser.add_argument('--region', type=str, default=None,
                        help='Genomic region to process (e.g., "chr1", "chr1:1000000-2000000"). Required if VCFs cover multiple contigs.')
    parser.add_argument('--memory', type=int, default=16,
                        help='Memory limit in GB for vg construct.')
    parser.add_argument('--threads', type=int, default=os.cpu_count() or 8,
                        help='Number of threads to use for vg construct.') # Added threads argument

    args = parser.parse_args()

    # Check prerequisites
    if not check_prerequisites():
        sys.exit(1)

    # Validate inputs
    for vcf_path in args.vcf:
        if not Path(vcf_path).exists():
            logger.error(f"Input VCF file not found: {vcf_path}")
            sys.exit(1)
        # Optional: Check for index file if needed by vg construct (often requires .tbi)
        # index_path = Path(f"{vcf_path}.tbi")
        # if not index_path.exists():
        #     logger.warning(f"Index file (.tbi) not found for {vcf_path}. vg construct might require it.")

    ref_path = Path(args.reference)
    if not ref_path.exists():
        logger.error(f"Reference FASTA file not found: {args.reference}")
        sys.exit(1)
    # Optional: Check for faidx index
    # faidx_path = Path(f"{args.reference}.fai")
    # if not faidx_path.exists():
    #     logger.warning(f"FASTA index file (.fai) not found for {args.reference}. vg construct might require it.")


    # --- Execution Logic ---
    all_success = True
    if args.mode == 'joint':
        output_path = Path(args.output)
        # If output is a directory, create a default filename inside it
        if output_path.is_dir() or not output_path.suffix:
             output_path.mkdir(parents=True, exist_ok=True) # Ensure dir exists
             output_path = output_path / "joint_graph.gfa"
        else:
             # Ensure parent directory exists if output is a file path
             output_path.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Starting joint conversion for {len(args.vcf)} VCF files...")
        success = convert_vcf_to_gfa_vg(
            vcf_paths=args.vcf,
            reference_path=args.reference,
            output_path=output_path,
            region=args.region,
            memory_gb=args.memory
            # threads=args.threads # Pass threads if function signature is updated
        )
        if not success:
            logger.error("Joint conversion failed.")
            all_success = False

    else: # mode == 'separate'
        output_dir = Path(args.output)
        # If output exists and is a file, use its parent directory
        if output_dir.exists() and output_dir.is_file():
             logger.warning(f"Output path '{args.output}' is a file, using its parent directory '{output_dir.parent}' for separate outputs.")
             output_dir = output_dir.parent

        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Starting separate conversion for {len(args.vcf)} VCF files into directory: {output_dir}")

        for vcf_path_str in args.vcf:
            vcf_path = Path(vcf_path_str)
            # Create a descriptive name for the output GFA
            base_name = vcf_path.name.replace('.vcf.gz', '').replace('.vcf', '')
            output_gfa_path = output_dir / f"{base_name}.gfa"

            logger.info(f"Converting {vcf_path.name} -> {output_gfa_path}...")
            success = convert_vcf_to_gfa_vg(
                vcf_paths=[vcf_path_str], # Pass only one VCF at a time
                reference_path=args.reference,
                output_path=output_gfa_path,
                region=args.region,
                memory_gb=args.memory
                # threads=args.threads # Pass threads if function signature is updated
            )
            if not success:
                logger.error(f"Failed to convert {vcf_path.name}")
                all_success = False
            else:
                logger.info(f"Successfully converted {vcf_path.name}")

    # --- Final Status ---
    if all_success:
        logger.info("Workflow finished successfully.")
        sys.exit(0)
    else:
        logger.error("Workflow finished with errors.")
        sys.exit(1)

if __name__ == "__main__":
    main()
