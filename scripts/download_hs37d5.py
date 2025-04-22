#!/usr/bin/env python3
"""
Downloads the hs37d5 reference genome FASTA file and its index
from the 1000 Genomes FTP site.
"""

import argparse
import ftplib
import gzip
import logging
import os
import shutil
import sys
from pathlib import Path
from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

FTP_SERVER = "ftp.1000genomes.ebi.ac.uk"
FTP_DIR = "/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/"
FASTA_FILENAME_GZ = "hs37d5.fa.gz"
FASTA_FILENAME = "hs37d5.fa"
INDEX_FILENAME = "hs37d5.fa.fai"

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

def main():
    parser = argparse.ArgumentParser(description=f"Download the {FASTA_FILENAME} reference genome and index.")
    parser.add_argument('--output-dir', '-o', required=True, help='Directory to save the downloaded files.')
    parser.add_argument('--extract', action='store_true', help='Extract the downloaded FASTA file.')
    parser.add_argument('--skip-existing', action='store_true', help='Skip download/extraction if final files exist.')


    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_output_path_gz = output_dir / FASTA_FILENAME_GZ
    fasta_output_path = output_dir / FASTA_FILENAME
    index_output_path = output_dir / INDEX_FILENAME

    # --- Skip Logic ---
    final_fasta_exists = fasta_output_path.exists() and fasta_output_path.is_file() and os.path.getsize(fasta_output_path) > 0
    final_index_exists = index_output_path.exists() and index_output_path.is_file() and os.path.getsize(index_output_path) > 0

    if args.skip_existing and final_fasta_exists and final_index_exists:
        logger.info(f"Final files {fasta_output_path.name} and {index_output_path.name} already exist in {output_dir}. Skipping download and extraction.")
        sys.exit(0)
    elif args.skip_existing and args.extract and final_fasta_exists:
         logger.info(f"Final FASTA file {fasta_output_path.name} already exists.")
         # Still need to check/download index below
    elif args.skip_existing and not args.extract and fasta_output_path_gz.exists():
         logger.info(f"Compressed FASTA file {fasta_output_path_gz.name} already exists.")
         # Still need to check/download index below


    # --- Download FASTA ---
    fasta_download_target = fasta_output_path_gz
    if args.extract and final_fasta_exists:
        logger.info(f"Skipping FASTA download as extracted file {fasta_output_path.name} exists.")
        fasta_download_success = True # Treat as success for logic flow
    elif not args.extract and fasta_output_path_gz.exists():
         logger.info(f"Skipping FASTA download as compressed file {fasta_output_path_gz.name} exists.")
         fasta_download_success = True # Treat as success for logic flow
    else:
        logger.info("Downloading FASTA file...")
        fasta_download_success = download_file_with_progress(FTP_SERVER, FTP_DIR, FASTA_FILENAME_GZ, fasta_download_target)
        if not fasta_download_success:
            logger.error("FASTA download failed. Exiting.")
            sys.exit(1)

    # --- Download Index ---
    if final_index_exists:
        logger.info(f"Skipping index download as file {index_output_path.name} exists.")
        index_download_success = True
    else:
        logger.info("Downloading FASTA index file...")
        index_download_success = download_file_with_progress(FTP_SERVER, FTP_DIR, INDEX_FILENAME, index_output_path)
        if not index_download_success:
            logger.error("FASTA index download failed. Exiting.")
            sys.exit(1)


    # --- Extract FASTA ---
    if args.extract:
        if final_fasta_exists:
             logger.info(f"Skipping extraction as {fasta_output_path.name} already exists.")
        elif fasta_output_path_gz.exists():
             extract_success = extract_gzip(fasta_output_path_gz, fasta_output_path, remove_gz=True)
             if not extract_success:
                 logger.error("FASTA extraction failed. Exiting.")
                 sys.exit(1)
        elif not fasta_download_success: # Should not happen if download failed earlier, but check anyway
             logger.error("Cannot extract FASTA as download failed or was skipped and file doesn't exist.")
             sys.exit(1)
        else:
             logger.warning(f"Cannot extract FASTA. Gzip file {fasta_output_path_gz.name} not found (and extracted version doesn't exist).")
             # Exit with error if extraction was requested but couldn't happen
             sys.exit(1)


    logger.info("Reference genome download (and extraction, if requested) completed successfully.")
    sys.exit(0)


if __name__ == "__main__":
    main()
