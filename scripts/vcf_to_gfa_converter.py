#!/usr/bin/env python3
"""
Script to convert VCF files to GFA format using vg.
Allows converting multiple VCF files either jointly or separately.
Includes fail-fast error handling and resume capability.
"""

import argparse
import logging
import os
import sys
import subprocess
from pathlib import Path
import tempfile
import shutil # Import shutil for which

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def check_prerequisites():
    """
    Check if required external tools are installed.
    """
    required_tools = ['vg', 'samtools'] # Added samtools
    missing_tools = []

    for tool in required_tools:
        tool_path = shutil.which(tool)
        if tool_path:
            logger.info(f"Found {tool} at {tool_path}")
        else:
            logger.warning(f"{tool} not found in PATH.")
            missing_tools.append(tool)

    if missing_tools:
        logger.error(f"Missing required tools: {', '.join(missing_tools)}")
        logger.error("Please install these tools and ensure they are in your PATH.")
        if 'vg' in missing_tools:
             logger.error("  'vg' can be installed from https://github.com/vgteam/vg")
        if 'samtools' in missing_tools:
             logger.error("  'samtools' can be installed from https://github.com/samtools/samtools")
        return False

    return True

def convert_vcf_to_gfa_vg(vcf_paths, reference_path, output_path, region=None, memory_gb=16, threads=None):
    """
    Convert VCF to GFA using the vg toolkit.

    Args:
        vcf_paths: List of VCF file paths (bgzipped and indexed).
        reference_path: Path to reference FASTA (indexed).
        output_path: Output GFA file path.
        region: Optional region to restrict conversion (e.g., "chr1:1000-2000").
        memory_gb: Memory limit for vg construct in gigabytes.
        threads: Number of threads for vg construct.

    Returns:
        True if conversion was successful, False otherwise.
    """
    temp_dir_base = None # Initialize to prevent UnboundLocalError in finally
    try:
        # Ensure output directory exists before potentially creating temp files or output
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

        # Create temporary directory for intermediate files
        # Use a directory relative to the output path for better cleanup/permissions
        temp_dir_base = Path(output_path).parent / "tmp_vg_convert"
        temp_dir_base.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=f"{Path(output_path).stem}_", dir=temp_dir_base) as temp_dir:
            logger.info(f"Using temporary directory: {temp_dir}")
            # Construct vg command for multiple VCFs
            vcf_args = []
            for vcf_path in vcf_paths:
                # Check VCF index exists
                tbi_path = Path(f"{vcf_path}.tbi")
                if not tbi_path.exists():
                     logger.error(f"Required index file (.tbi) not found for {vcf_path}")
                     return False
                vcf_args.extend(['-v', str(vcf_path)])

            # Check reference index exists
            fai_path = Path(f"{reference_path}.fai")
            if not fai_path.exists():
                 logger.warning(f"Reference index file (.fai) not found for {reference_path}. Attempting to generate it...")
                 try:
                     samtools_exe = shutil.which('samtools')
                     if not samtools_exe:
                          logger.error("samtools command not found in PATH. Cannot generate FASTA index.")
                          return False
                     subprocess.run([samtools_exe, 'faidx', str(reference_path)], check=True, capture_output=True)
                     logger.info(f"Successfully generated index for {reference_path}")
                 except (subprocess.CalledProcessError) as idx_err:
                     logger.error(f"Failed to generate FASTA index for {reference_path}. 'samtools faidx' failed.")
                     logger.error(f"Error details: {idx_err.stderr.decode(errors='replace') if idx_err.stderr else 'No stderr'}")
                     return False


            region_args = []
            if region:
                # vg construct expects region like 'chr:start-end' or just 'chr'
                region_args = ['-R', region]

            # Determine number of threads
            num_threads = threads if threads is not None else (os.cpu_count() or 8)

            # Step 1: Construct the graph using vg construct
            vg_path = os.path.join(temp_dir, "graph.vg")
            vg_exe = shutil.which('vg') # Find vg executable
            if not vg_exe:
                 logger.error("vg command not found in PATH.")
                 return False

            construct_cmd = [
                vg_exe, 'construct',
                '-r', str(reference_path),
                *vcf_args,
                '-m', str(memory_gb), # Memory in GB - vg construct expects an integer value
                *region_args,
                '-t', str(num_threads),
                # '-a' # Add -a if you want all alleles included, otherwise it's ref+alt
                # Consider adding --validate to check graph structure
            ]

            logger.info(f"Running vg construct: {' '.join(construct_cmd)} > {vg_path}")
            # Redirect stdout to the vg_path file
            try:
                with open(vg_path, 'wb') as vg_file: # Write vg output as bytes
                    result = subprocess.run(construct_cmd, stdout=vg_file, stderr=subprocess.PIPE, check=False) # Capture stderr as bytes
                    # Check return code explicitly
                    if result.returncode != 0:
                         # Decode stderr for logging, replacing errors
                         stderr_decoded = result.stderr.decode(errors='replace') if result.stderr else ""
                         raise subprocess.CalledProcessError(result.returncode, construct_cmd, output=None, stderr=stderr_decoded)
                    # Log stderr even on success (might contain warnings)
                    if result.stderr:
                        logger.info(f"vg construct stderr: {result.stderr.decode(errors='replace').strip()}")

            except subprocess.CalledProcessError as e:
                 logger.error(f"vg construct command failed: {' '.join(e.cmd)}")
                 logger.error(f"Return code: {e.returncode}")
                 # stderr is already decoded in the raised exception if possible
                 if e.stderr:
                     logger.error(f"stderr: {e.stderr.strip()}")
                 # Attempt to provide common troubleshooting tips based on error
                 if "variant/reference sequence mismatch" in str(e.stderr):
                     logger.error("Potential Issue: VCF coordinates might not match the reference FASTA assembly version.")
                 elif "index file" in str(e.stderr) and "not found" in str(e.stderr):
                     logger.error("Potential Issue: Ensure VCF files are bgzipped and indexed (.tbi) and reference FASTA is indexed (.fai).")
                 return False # Exit function on command failure

            # Check if the vg file was created and is not empty
            if not os.path.exists(vg_path) or os.path.getsize(vg_path) == 0:
                logger.error(f"vg construct failed to create a valid graph file: {vg_path}")
                # Attempt to log stderr if available from a previous failure capture
                if 'result' in locals() and hasattr(result, 'stderr') and result.stderr:
                     logger.error(f"vg construct stderr: {result.stderr.decode(errors='replace').strip()}")
                return False

            # Step 1.5: Validate the intermediate graph
            validate_cmd = [vg_exe, 'validate', vg_path]
            logger.info(f"Running vg validate: {' '.join(validate_cmd)}")
            try:
                result_validate = subprocess.run(validate_cmd, capture_output=True, check=False) # Capture output
                # Decode stderr for logging, replacing errors
                stderr_decoded = result_validate.stderr.decode(errors='replace') if result_validate.stderr else ""
                if result_validate.returncode != 0:
                    logger.error(f"vg validate command failed for {vg_path}")
                    logger.error(f"Return code: {result_validate.returncode}")
                    if stderr_decoded:
                        logger.error(f"stderr: {stderr_decoded.strip()}")
                    # Optionally log stdout as well if it contains useful info on failure
                    stdout_decoded = result_validate.stdout.decode(errors='replace') if result_validate.stdout else ""
                    if stdout_decoded:
                         logger.error(f"stdout: {stdout_decoded.strip()}")
                    return False # Exit if validation fails
                else:
                    logger.info(f"Intermediate graph {vg_path} validated successfully.")
                    if stderr_decoded: # Log warnings from validate even on success
                         logger.info(f"vg validate stderr: {stderr_decoded.strip()}")

            except Exception as e:
                 logger.error(f"An unexpected error occurred during vg validate: {e}")
                 return False


            # Step 2: Convert to GFA using vg view
            gfa_cmd = [
                vg_exe, 'view',
                vg_path,
                '-F', '-g',  # Output in GFA format
                # Consider adding --validate to check GFA structure
            ]

            logger.info(f"Running vg view: {' '.join(gfa_cmd)} > {output_path}")
            # Redirect stdout to the final output_path GFA file
            try:
                # Write stdout directly as bytes, capture stderr as bytes
                with open(output_path, 'wb') as gfa_file:
                     result_view = subprocess.run(gfa_cmd, stdout=gfa_file, stderr=subprocess.PIPE, check=False) # text=False is default
                     # Decode stderr for logging, replacing errors
                     stderr_decoded = result_view.stderr.decode(errors='replace') if result_view.stderr else ""
                     if result_view.returncode != 0:
                          raise subprocess.CalledProcessError(result_view.returncode, gfa_cmd, output=None, stderr=stderr_decoded)
                     if stderr_decoded:
                         logger.info(f"vg view stderr: {stderr_decoded.strip()}") # Log non-error stderr for info

            except subprocess.CalledProcessError as e:
                 logger.error(f"vg view command failed: {' '.join(e.cmd)}")
                 logger.error(f"Return code: {e.returncode}")
                 # stderr is already decoded in the raised exception if possible
                 if e.stderr:
                     logger.error(f"stderr: {e.stderr.strip()}")
                 # Clean up potentially incomplete output file
                 if os.path.exists(output_path):
                     os.remove(output_path)
                     logger.info(f"Removed potentially incomplete output file: {output_path}")
                 return False # Exit function on command failure


            # Check if the GFA file was created and is not empty
            if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
                 logger.error(f"vg view failed to create a valid GFA file: {output_path}")
                 if 'result_view' in locals() and hasattr(result_view, 'stderr') and result_view.stderr:
                     logger.error(f"vg view stderr: {result_view.stderr.decode(errors='replace').strip()}")
                 # Clean up potentially incomplete output file
                 if os.path.exists(output_path):
                     os.remove(output_path)
                     logger.info(f"Removed empty output file: {output_path}")
                 return False

            logger.info(f"Successfully created GFA file: {output_path}")
            return True

    except Exception as e:
        logger.error(f"An unexpected error occurred during VCF to GFA conversion for output {output_path}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        # Clean up potentially incomplete output file
        if 'output_path' in locals() and os.path.exists(output_path):
             # Check if it was created by this function run or existed before
             # This is tricky, maybe just attempt removal if size is 0 or it was touched recently?
             # For simplicity, let's remove if it exists after an error here.
             try:
                 os.remove(output_path)
                 logger.info(f"Removed potentially incomplete output file due to error: {output_path}")
             except OSError as rm_err:
                 logger.warning(f"Could not remove potentially incomplete file {output_path}: {rm_err}")
        return False
    finally:
         # Clean up the base temporary directory if it exists and is empty
         if temp_dir_base and temp_dir_base.exists():
             try:
                 # Check if empty before removing
                 if not any(temp_dir_base.iterdir()):
                     temp_dir_base.rmdir()
                     logger.debug(f"Removed empty base temporary directory: {temp_dir_base}")
             except OSError as e:
                 logger.warning(f"Could not remove base temporary directory {temp_dir_base}: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert VCF files to GFA format using vg. Supports joint or separate conversion, fail-fast, and resume.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # Show default values
        )
    parser.add_argument('--vcf', '-v', required=True, action='append',
                        help='Input VCF file path(s). Can be specified multiple times. Files should be bgzipped and indexed (.tbi).')
    parser.add_argument('--reference', '-r', required=True,
                        help='Reference genome FASTA file path. Should be indexed (.fai) or indexable by samtools.')
    parser.add_argument('--output', '-o', required=True,
                        help='Output GFA file path (for joint mode) or output directory (for separate mode).')
    parser.add_argument('--mode', choices=['joint', 'separate'], default='joint',
                        help='Conversion mode: "joint" combines all VCFs into one GFA, "separate" creates one GFA per VCF.')
    parser.add_argument('--region', type=str, default=None,
                        help='Optional: Genomic region to process (e.g., "chr1", "chr1:1000000-2000000"). Required by vg construct if VCFs cover multiple contigs without explicit contig lines.')
    parser.add_argument('--memory', type=int, default=16,
                        help='Memory limit in GB for vg construct.')
    parser.add_argument('--threads', type=int, default=os.cpu_count() or 8,
                        help='Number of threads to use for vg construct.')
    # Add skip-existing flag if needed by the caller (workflow manager handles this now)
    # parser.add_argument('--skip-existing', action='store_true', help='Skip conversion if output file already exists.')


    args = parser.parse_args()

    # --- Prerequisites Check ---
    if not check_prerequisites():
        logger.error("Prerequisite check failed. Exiting.")
        sys.exit(1)

    # --- Input Validation ---
    valid_inputs = True
    for vcf_path_str in args.vcf:
        vcf_path = Path(vcf_path_str)
        if not vcf_path.exists():
            logger.error(f"Input VCF file not found: {vcf_path}")
            valid_inputs = False
        # Check for index file - vg construct usually requires .tbi for region processing
        index_path = Path(f"{vcf_path_str}.tbi")
        if not index_path.exists():
            logger.warning(f"Index file (.tbi) not found for {vcf_path}. vg construct might require it, especially with --region.")
            # Consider making this an error if --region is specified:
            # if args.region:
            #     logger.error(f"Index file (.tbi) is required for {vcf_path} when using --region.")
            #     valid_inputs = False

    ref_path = Path(args.reference)
    if not ref_path.exists():
        logger.error(f"Reference FASTA file not found: {args.reference}")
        valid_inputs = False
    # Check for faidx index - vg construct requires it
    faidx_path = Path(f"{args.reference}.fai")
    if not faidx_path.exists():
        logger.warning(f"FASTA index file (.fai) not found for {args.reference}. vg construct requires it. Will attempt to generate using samtools.")
        # The convert function will handle generation, but warn here.

    if not valid_inputs:
        logger.error("Input validation failed. Exiting.")
        sys.exit(1)


    # --- Execution Logic ---
    if args.mode == 'joint':
        output_path = Path(args.output)
        # If output is specified as a directory, create a default filename inside it
        if output_path.is_dir() or (not output_path.suffix and not output_path.exists()):
             output_path.mkdir(parents=True, exist_ok=True) # Ensure dir exists
             output_path = output_path / "joint_graph.gfa"
             logger.info(f"Output path is a directory, will write joint GFA to: {output_path}")
        else:
             # Ensure parent directory exists if output is a file path
             output_path.parent.mkdir(parents=True, exist_ok=True)

        # Resume Check (Handled by Workflow Manager now)
        # if args.skip_existing and output_path.exists() and output_path.is_file() and os.path.getsize(output_path) > 0:
        #     logger.info(f"Output file {output_path} already exists and is not empty. Skipping joint conversion.")
        #     sys.exit(0) # Exit successfully as the target file exists

        logger.info(f"Starting joint conversion for {len(args.vcf)} VCF files...")
        success = convert_vcf_to_gfa_vg(
            vcf_paths=args.vcf,
            reference_path=args.reference,
            output_path=output_path,
            region=args.region,
            memory_gb=args.memory,
            threads=args.threads
        )
        if not success:
            logger.error(f"Joint conversion failed. See logs above for details. Exiting.")
            sys.exit(1) # Exit immediately on failure

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
            # Handle potential double extensions like .vcf.gz
            base_name = vcf_path.name
            if base_name.endswith(".vcf.gz"):
                base_name = base_name[:-len(".vcf.gz")]
            elif base_name.endswith(".vcf"):
                 base_name = base_name[:-len(".vcf")]
            output_gfa_path = output_dir / f"{base_name}.gfa"

            # Resume Check (Handled by Workflow Manager now)
            # if args.skip_existing and output_gfa_path.exists() and output_gfa_path.is_file() and os.path.getsize(output_gfa_path) > 0:
            #     logger.info(f"Output file {output_gfa_path} already exists and is not empty. Skipping conversion for {vcf_path.name}.")
            #     continue # Skip to the next file

            logger.info(f"Converting {vcf_path.name} -> {output_gfa_path}...")
            success = convert_vcf_to_gfa_vg(
                vcf_paths=[vcf_path_str], # Pass only one VCF at a time
                reference_path=args.reference,
                output_path=output_gfa_path,
                region=args.region,
                memory_gb=args.memory,
                threads=args.threads
            )
            if not success:
                logger.error(f"Failed to convert {vcf_path.name}. See logs above for details. Exiting.")
                sys.exit(1) # Exit immediately on failure
            else:
                logger.info(f"Successfully converted {vcf_path.name}")

    # --- Final Status ---
    logger.info("Workflow step finished successfully.") # Changed message slightly
    sys.exit(0) # Exit successfully if all steps completed or were skipped

if __name__ == "__main__":
    main()
