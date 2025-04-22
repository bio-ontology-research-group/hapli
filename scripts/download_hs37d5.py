#!/usr/bin/env python3
"""
Downloads the hs37d5 reference genome FASTA file from the 1000 Genomes FTP site
and generates the index locally using samtools.
"""

import argparse
import ftplib
import gzip
import logging
import os
import shutil
import sys
import subprocess # Added for samtools
from pathlib import Path
from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

FTP_SERVER = "ftp.1000genomes.ebi.ac.uk"
FTP_DIR = "/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/"
FASTA_FILENAME_GZ = "hs37d5.fa.gz"
FASTA_FILENAME = "hs37d5.fa"
# INDEX_FILENAME = "hs37d5.fa.fai" # No longer downloading index

def check_samtools_available():
    """Check if samtools is available in the path."""
    samtools_path = shutil.which("samtools")
    if samtools_path:
        logger.info(f"Found samtools at {samtools_path}")
        return True
    else:
        logger.error("samtools not found in PATH. It is required to generate the FASTA index.")
        return False

def download_file_with_progress(server: str, directory: str, filename: str, output_path: Path):
    """Downloads a file from an FTP server with a progress bar."""
    try:
        ftp = ftplib.FTP(server)
        ftp.login()  # Anonymous login
        logger.info(f"Connected to FTP server: {server}")
        ftp.cwd(directory)
        logger.info(f"Changed directory to: {directory}")

        # Get file size for progress bar
        try:
            file_size = ftp.size(filename)
            if file_size is None:
                 logger.warning(f"Could not determine size of {filename}. Progress bar may not be accurate.")
                 file_size = 0 # Set to 0 if size cannot be determined
        except ftplib.error_perm as e:
             logger.warning(f"FTP error getting size for {filename}: {e}. Progress bar may not be accurate.")
             file_size = 0


        logger.info(f"Downloading {filename} ({file_size / (1024*1024):.2f} MB) to {output_path}...")

        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'wb') as fp, \
             tqdm(total=file_size, unit='B', unit_scale=True, unit_divisor=1024,
                  desc=f"Downloading {filename}", leave=False, disable=None) as pbar:
            try:
                def callback(data):
                    fp.write(data)
                    pbar.update(len(data))

                ftp.retrbinary(f'RETR {filename}', callback)
            except ftplib.error_perm as e:
                 logger.error(f"FTP error during download of {filename}: {e}")
                 ftp.quit()
                 # Clean up potentially incomplete file
                 if output_path.exists():
                     output_path.unlink()
                 return False
            except Exception as e:
                 logger.error(f"Unexpected error during download of {filename}: {e}")
                 ftp.quit()
                 # Clean up potentially incomplete file
                 if output_path.exists():
                     output_path.unlink()
                 return False


        ftp.quit()
        logger.info(f"Successfully downloaded {filename}")
        return True

    except ftplib.all_errors as e:
        logger.error(f"FTP connection or command failed: {e}")
        # Clean up potentially incomplete file
        if output_path.exists():
            output_path.unlink()
        return False
    except Exception as e:
        logger.error(f"An unexpected error occurred during FTP download: {e}")
        # Clean up potentially incomplete file
        if output_path.exists():
            output_path.unlink()
        return False

def extract_gzip(gzip_path: Path, output_path: Path, remove_gz: bool = True):
    """Extracts a gzipped file."""
    logger.info(f"Extracting {gzip_path} to {output_path}...")
    try:
        with gzip.open(gzip_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        logger.info(f"Successfully extracted {output_path}")
        if remove_gz:
            gzip_path.unlink()
            logger.info(f"Removed original gzip file: {gzip_path}")
        return True
    except gzip.BadGzipFile:
        logger.error(f"Error: {gzip_path} is not a valid gzip file or is corrupted.")
        # Clean up potentially incomplete output file
        if output_path.exists():
            output_path.unlink()
        return False
    except Exception as e:
        logger.error(f"Failed to extract {gzip_path}: {e}")
        # Clean up potentially incomplete output file
        if output_path.exists():
            output_path.unlink()
        return False

def generate_fasta_index(fasta_path: Path) -> bool:
    """Generates a FASTA index (.fai) using samtools faidx."""
    index_path = fasta_path.with_suffix(fasta_path.suffix + '.fai')
    logger.info(f"Generating FASTA index for {fasta_path.name} -> {index_path.name}...")
    try:
        cmd = ['samtools', 'faidx', str(fasta_path)]
        logger.debug(f"Running command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        if not index_path.exists() or index_path.stat().st_size == 0:
             logger.error(f"samtools faidx command finished but index file {index_path} is missing or empty.")
             return False

        logger.info(f"Successfully generated FASTA index: {index_path.name}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"samtools faidx failed for {fasta_path.name}: {e.stderr}")
        return False
    except FileNotFoundError:
        # This should be caught by check_samtools_available, but check again
        logger.error("samtools command not found. Cannot generate index.")
        return False
    except Exception as e:
        logger.error(f"An unexpected error occurred during index generation: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description=f"Download the {FASTA_FILENAME} reference genome and generate its index.")
    parser.add_argument('--output-dir', '-o', required=True, help='Directory to save the downloaded files.')
    parser.add_argument('--extract', action='store_true', help='Extract the downloaded FASTA file (required for indexing).')
    parser.add_argument('--skip-existing', action='store_true', help='Skip download/extraction/indexing if final files exist.')

    args = parser.parse_args()

    # --- Prerequisite Check ---
    if not check_samtools_available():
        sys.exit(1)

    if not args.extract:
        logger.warning("--extract flag is not set. FASTA index cannot be generated without extracting the FASTA file. Please add --extract.")
        # Decide whether to exit or proceed without indexing
        # For now, let's exit as the goal is usually to have the index.
        sys.exit(1)


    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_output_path_gz = output_dir / FASTA_FILENAME_GZ
    fasta_output_path = output_dir / FASTA_FILENAME
    index_output_path = fasta_output_path.with_suffix(fasta_output_path.suffix + '.fai') # Define expected index path

    # --- Skip Logic ---
    # Check if both the final FASTA (.fa) and its index (.fa.fai) exist
    final_fasta_exists = fasta_output_path.exists() and fasta_output_path.is_file() and os.path.getsize(fasta_output_path) > 0
    final_index_exists = index_output_path.exists() and index_output_path.is_file() and os.path.getsize(index_output_path) > 0

    if args.skip_existing and final_fasta_exists and final_index_exists:
        logger.info(f"Final files {fasta_output_path.name} and {index_output_path.name} already exist in {output_dir}. Skipping.")
        sys.exit(0)

    # --- Download FASTA ---
    fasta_download_target = fasta_output_path_gz
    fasta_available_for_indexing = final_fasta_exists # Is the .fa file already there?

    if final_fasta_exists:
        logger.info(f"Skipping FASTA download as extracted file {fasta_output_path.name} exists.")
        fasta_download_success = True
    elif fasta_output_path_gz.exists():
         logger.info(f"Compressed FASTA file {fasta_output_path_gz.name} exists. Will proceed to extraction.")
         fasta_download_success = True # Download itself is skipped, but file is present
    else:
        logger.info("Downloading FASTA file...")
        fasta_download_success = download_file_with_progress(FTP_SERVER, FTP_DIR, FASTA_FILENAME_GZ, fasta_download_target)
        if not fasta_download_success:
            logger.error("FASTA download failed. Exiting.")
            sys.exit(1)

    # --- Extract FASTA ---
    # Extraction is required for indexing, args.extract is checked at the start
    if not final_fasta_exists: # Only extract if .fa doesn't exist
        if fasta_output_path_gz.exists():
             extract_success = extract_gzip(fasta_output_path_gz, fasta_output_path, remove_gz=True)
             if not extract_success:
                 logger.error("FASTA extraction failed. Exiting.")
                 sys.exit(1)
             fasta_available_for_indexing = True
        else:
             # This case should ideally not be reached if download succeeded
             logger.error(f"Cannot extract FASTA. Gzip file {fasta_output_path_gz.name} not found.")
             sys.exit(1)
    else:
         logger.info(f"Skipping extraction as {fasta_output_path.name} already exists.")
         fasta_available_for_indexing = True


    # --- Generate Index ---
    if fasta_available_for_indexing:
        if final_index_exists:
            logger.info(f"Skipping index generation as {index_output_path.name} already exists.")
        else:
            index_success = generate_fasta_index(fasta_output_path)
            if not index_success:
                logger.error("FASTA index generation failed. Exiting.")
                sys.exit(1)
    else:
        # Should not happen if logic above is correct
        logger.error("FASTA file is not available for indexing. Exiting.")
        sys.exit(1)


    logger.info("Reference genome download, extraction, and indexing completed successfully.")
    sys.exit(0)


if __name__ == "__main__":
    main()
