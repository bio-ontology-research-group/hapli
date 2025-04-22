#!/usr/bin/env python3
"""
Downloads GRCh37 GFF3 annotation files from Ensembl (Release 113, Build 87 annotations).
Currently downloads only chromosome 22.
"""

import argparse
import gzip
import logging
import os
import shutil
import sys
import requests
from pathlib import Path
from tqdm import tqdm
from urllib.parse import urlparse

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Ensembl Release 113 (GRCh37, Build 87 annotations) - Chromosome 22 specific
ENSEMBL_BASE_URL = "http://ftp.ensembl.org/pub/grch37/release-113/gff3/homo_sapiens/"
# TODO: Add argument to specify chromosome if needed later
GFF3_FILENAME_GZ = "Homo_sapiens.GRCh37.87.chromosome.22.gff3.gz"
GFF3_FILENAME = "Homo_sapiens.GRCh37.87.chromosome.22.gff3"

def download_file_with_progress(url: str, output_path: Path):
    """Downloads a file from HTTP/FTP with a progress bar."""
    logger.info(f"Attempting to download from: {url}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        response = requests.get(url, stream=True, timeout=60) # Added timeout
        response.raise_for_status()  # Raise an exception for bad status codes

        total_size = int(response.headers.get('content-length', 0))
        block_size = 8192 # 8 Kibibytes

        with open(output_path, 'wb') as f, \
             tqdm(total=total_size, unit='B', unit_scale=True, unit_divisor=1024,
                  desc=f"Downloading {output_path.name}", leave=False, disable=None) as pbar:
            for chunk in response.iter_content(chunk_size=block_size):
                f.write(chunk)
                pbar.update(len(chunk))

        # Verify download size
        if total_size != 0 and pbar.n != total_size:
             logger.error(f"Download incomplete: Expected {total_size} bytes, got {pbar.n} bytes.")
             # Clean up incomplete file
             if output_path.exists():
                 output_path.unlink()
             return False

        logger.info(f"Successfully downloaded {output_path.name}")
        return True

    except requests.exceptions.RequestException as e:
        logger.error(f"Error downloading {url}: {e}")
        # Clean up potentially incomplete file
        if output_path.exists():
            output_path.unlink()
        return False
    except Exception as e:
        logger.error(f"An unexpected error occurred during download: {e}")
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
    parser = argparse.ArgumentParser(description=f"Download GRCh37 GFF3 annotation for Chromosome 22 from Ensembl Release 113.")
    parser.add_argument('--output-dir', '-o', required=True, help='Directory to save the downloaded GFF3 file.')
    parser.add_argument('--extract', action='store_true', help='Extract the downloaded GFF3 file.')
    parser.add_argument('--skip-existing', action='store_true', help='Skip download/extraction if final file exists.')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gff3_url = f"{ENSEMBL_BASE_URL}{GFF3_FILENAME_GZ}"
    gff3_output_path_gz = output_dir / GFF3_FILENAME_GZ
    gff3_output_path = output_dir / GFF3_FILENAME

    # --- Skip Logic ---
    final_gff3_exists = gff3_output_path.exists() and gff3_output_path.is_file() and os.path.getsize(gff3_output_path) > 0

    if args.skip_existing:
        if args.extract and final_gff3_exists:
            logger.info(f"Final extracted file {gff3_output_path.name} already exists in {output_dir}. Skipping.")
            sys.exit(0)
        elif not args.extract and gff3_output_path_gz.exists() and os.path.getsize(gff3_output_path_gz) > 0:
             logger.info(f"Final compressed file {gff3_output_path_gz.name} already exists in {output_dir}. Skipping.")
             sys.exit(0)
        # If only one exists but the other is needed, proceed below

    # --- Download GFF3 ---
    download_needed = True
    if args.extract and final_gff3_exists:
        logger.info(f"Skipping GFF3 download as extracted file {gff3_output_path.name} exists.")
        download_needed = False
    elif not args.extract and gff3_output_path_gz.exists():
         logger.info(f"Skipping GFF3 download as compressed file {gff3_output_path_gz.name} exists.")
         download_needed = False

    if download_needed:
        logger.info("Downloading GFF3 file...")
        download_success = download_file_with_progress(gff3_url, gff3_output_path_gz)
        if not download_success:
            logger.error("GFF3 download failed. Exiting.")
            sys.exit(1)
    else:
        download_success = True # Indicate success for logic flow if skipped


    # --- Extract GFF3 ---
    if args.extract:
        if final_gff3_exists:
             logger.info(f"Skipping extraction as {gff3_output_path.name} already exists.")
        elif gff3_output_path_gz.exists():
             extract_success = extract_gzip(gff3_output_path_gz, gff3_output_path, remove_gz=True)
             if not extract_success:
                 logger.error("GFF3 extraction failed. Exiting.")
                 sys.exit(1)
        elif not download_success: # Should not happen if download failed earlier
             logger.error("Cannot extract GFF3 as download failed or was skipped and file doesn't exist.")
             sys.exit(1)
        else:
             logger.warning(f"Cannot extract GFF3. Gzip file {gff3_output_path_gz.name} not found (and extracted version doesn't exist).")
             # Exit with error if extraction was requested but couldn't happen
             sys.exit(1)

    logger.info("GFF3 download (and extraction, if requested) completed successfully.")
    sys.exit(0)

if __name__ == "__main__":
    main()
