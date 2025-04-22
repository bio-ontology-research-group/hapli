# src/converters/vcf_to_gfa_cli.py
import argparse
import logging
import sys
import os
import time
from .vcf_to_gfa import VCFtoGFAConverter, VCFtoGFAConversionError
from .reference_handler import ReferenceHandlerError # Import specific error

# Setup basic logging config (can be overridden by args)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger("vcf_to_gfa_cli") # Specific logger name

def main():
    """Command-line interface for VCF to GFA conversion."""
    parser = argparse.ArgumentParser(
        description="Convert a VCF file to GFA format, constructing haplotype paths.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # Show defaults in help
    )

    # Required arguments
    parser.add_argument("-v", "--vcf", required=True, help="Path to the input VCF file (.vcf, .vcf.gz, .bcf)")
    parser.add_argument("-f", "--fasta", required=True, help="Path to the reference FASTA file (.fa, .fa.gz)")
    parser.add_argument("-o", "--output", required=True, help="Path for the output GFA file")

    # Optional arguments
    parser.add_argument(
        "--path-template",
        default="{sample}_hap{hap}",
        help="Template string for naming haplotype paths. Use {sample} and {hap} placeholders."
    )
    parser.add_argument(
        "--unphased-strategy",
        choices=['ref', 'alt', 'skip'], # Removed 'bubble' as it's not implemented
        default='ref',
        help="Strategy for handling unphased or missing genotypes: "
             "'ref' (use reference allele), "
             "'alt' (use first alternate allele if available, else ref), "
             "'skip' (effectively use reference, skip variant allele)."
    )
    parser.add_argument(
        "--region",
        default=None,
        help="Process only a specific region (e.g., 'chr1', 'chr1:1000-2000'). Requires VCF index (.tbi/.csi). Coordinates are 1-based, inclusive."
    )
    # TODO: Add filtering options (quality, type) if needed in the future

    parser.add_argument(
        "--log-level",
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help="Set the logging verbosity."
    )

    args = parser.parse_args()

    # --- Configure Logging ---
    try:
        log_level_numeric = getattr(logging, args.log_level.upper())
        # Set level for all loggers
        root_logger = logging.getLogger()
        root_logger.setLevel(log_level_numeric)
        # If handlers exist, update their level too (e.g., basicConfig handler)
        for handler in root_logger.handlers:
             handler.setLevel(log_level_numeric)
        logger.info(f"Logging level set to {args.log_level}")
    except AttributeError:
        logger.error(f"Invalid log level: {args.log_level}. Using INFO.")
        logging.getLogger().setLevel(logging.INFO)


    logger.info("Starting VCF to GFA conversion CLI...")
    start_time = time.time()
    logger.info(f"  VCF Input: {args.vcf}")
    logger.info(f"  FASTA Reference: {args.fasta}")
    logger.info(f"  GFA Output: {args.output}")
    logger.info(f"  Path Template: {args.path_template}")
    logger.info(f"  Unphased Strategy: {args.unphased_strategy}")
    logger.info(f"  Region: {args.region or 'All Contigs'}")


    # --- Input File Validation ---
    if not os.path.exists(args.vcf):
        logger.error(f"Input VCF file not found: {args.vcf}")
        sys.exit(1)
    if not os.path.exists(args.fasta):
        logger.error(f"Input FASTA file not found: {args.fasta}")
        sys.exit(1)
    # Check for FASTA index (pyfaidx will create it, but good to inform user)
    if not os.path.exists(args.fasta + ".fai"):
         # Check if directory is writable for index creation
         fasta_dir = os.path.dirname(os.path.abspath(args.fasta))
         if not os.access(fasta_dir, os.W_OK):
              logger.warning(f"FASTA index (.fai) not found for {args.fasta}, and directory {fasta_dir} may not be writable. Index creation might fail.")
         else:
              logger.info(f"FASTA index (.fai) not found for {args.fasta}. pyfaidx will attempt to create it.")
    # Check for VCF index if needed (pysam requires it for region fetching)
    if args.region:
         has_index = False
         # Check for both .tbi and .csi index files
         for ext in ['.tbi', '.csi']:
             if os.path.exists(args.vcf + ext):
                 has_index = True
                 logger.info(f"Found VCF index: {args.vcf + ext}")
                 break
         if not has_index:
              logger.error(f"Region ('{args.region}') specified, but required VCF index (.tbi or .csi) not found for {args.vcf}.")
              sys.exit(1)


    # --- Run Conversion ---
    exit_code = 0
    try:
        logger.info("Initializing converter...")
        # Use context manager for automatic resource cleanup
        with VCFtoGFAConverter(
            vcf_filepath=args.vcf,
            fasta_filepath=args.fasta,
            output_gfa_filepath=args.output,
            path_template=args.path_template,
            unphased_strategy=args.unphased_strategy
        ) as converter:
            logger.info("Starting conversion process...")
            converter.convert(region=args.region)
        elapsed_time = time.time() - start_time
        logger.info(f"Conversion finished successfully in {elapsed_time:.2f} seconds.")

    except (VCFtoGFAConversionError, ReferenceHandlerError, PhasingError, ValueError, FileNotFoundError) as e:
        logger.error(f"Conversion failed: {e}")
        exit_code = 1
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}", exc_info=True) # Show traceback for unexpected errors
        exit_code = 1
    finally:
        elapsed_time = time.time() - start_time
        if exit_code == 0:
             logger.debug(f"Total execution time: {elapsed_time:.2f} seconds.")
        else:
             logger.error(f"Execution failed after {elapsed_time:.2f} seconds.")
        sys.exit(exit_code)


if __name__ == "__main__":
    # Ensure the script can find modules in the 'src' directory when run directly
    script_dir = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.dirname(os.path.dirname(script_dir)) # Assumes src is two levels up
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)
        logger.debug(f"Added {src_dir} to sys.path")

    main()
