#!/usr/bin/env python3
"""
Script to convert a FASTA file to GFA format.
This creates a linear GFA representation of the reference genome.
"""

import argparse
import logging
import sys
import subprocess
from pathlib import Path
from Bio import SeqIO

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def fasta_to_gfa(fasta_path, gfa_path):
    """
    Convert a FASTA file to GFA format.
    
    Each sequence in the FASTA file becomes a segment in the GFA.
    """
    try:
        with open(gfa_path, 'w') as gfa_file:
            # Write GFA header
            gfa_file.write("H\tVN:Z:1.0\n")
            
            # Process each sequence in the FASTA file
            for record in SeqIO.parse(fasta_path, "fasta"):
                seq_id = record.id
                sequence = str(record.seq)
                
                # Write segment line
                gfa_file.write(f"S\t{seq_id}\t{sequence}\n")
                
                # Optionally add path for each chromosome
                gfa_file.write(f"P\t{seq_id}_path\t{seq_id}+\t*\n")
            
            logger.info(f"Converted {fasta_path} to GFA format at {gfa_path}")
            return True
    except Exception as e:
        logger.error(f"Error converting FASTA to GFA: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Convert FASTA to GFA format")
    parser.add_argument('--input', '-i', required=True, type=str,
                        help='Input FASTA file path')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output GFA file path (default: input_path.gfa)')
    parser.add_argument('--use-external', action='store_true', 
                        help='Use external tool for conversion if available')
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file does not exist: {input_path}")
        sys.exit(1)
    
    output_path = args.output
    if output_path is None:
        output_path = input_path.with_suffix('.gfa')
    else:
        output_path = Path(output_path)
    
    if args.use_external:
        # Try to use vg or other external tools if they're installed
        try:
            # Check if vg is installed
            try:
                subprocess.run(['vg', '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                logger.info("Using vg for FASTA to GFA conversion")
                
                # Use vg to convert FASTA to GFA
                cmd = [
                    'vg', 'construct',
                    '-r', str(input_path),
                    '-m', '32',  # Memory limit in GB
                    '>', str(output_path)
                ]
                subprocess.run(' '.join(cmd), shell=True, check=True)
                logger.info(f"Successfully converted {input_path} to {output_path} using vg")
                return
            except (subprocess.SubprocessError, FileNotFoundError):
                logger.warning("vg not found or error running vg. Falling back to internal conversion.")
                
            # Check for any other tools here if desired
            
        except Exception as e:
            logger.warning(f"Error using external tool: {e}. Falling back to internal conversion.")
    
    # Use internal conversion if external tools not available or not requested
    success = fasta_to_gfa(input_path, output_path)
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
