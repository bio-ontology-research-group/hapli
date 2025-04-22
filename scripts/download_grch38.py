#!/usr/bin/env python3
"""
Script to download the GRCh38 reference genome FASTA file.
"""

import argparse
import os
import logging
import ftplib
import gzip
import shutil
from pathlib import Path
from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Define FTP server details for Ensembl
ENSEMBL_FTP_SERVER = "ftp.ensembl.org"
ENSEMBL_FTP_DIR = "/pub/release-109/fasta/homo_sapiens/dna/"
ENSEMBL_FILENAME = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

# Define FTP server details for NCBI
NCBI_FTP_SERVER = "ftp.ncbi.nlm.nih.gov"
NCBI_FTP_DIR = "/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
NCBI_FILENAME = "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

def download_file_with_progress(server, directory, filename, output_path):
    """
    Download a file from an FTP server with a progress bar.
    """
    try:
        # Get file size first
        ftp = ftplib.FTP(server)
        ftp.login()
        ftp.cwd(directory)
        file_size = ftp.size(filename)
        ftp.quit()
        
        # Now download with progress tracking
        with open(output_path, 'wb') as f:
            logger.info(f"Connecting to {server}")
            ftp = ftplib.FTP(server)
            ftp.login()
            ftp.cwd(directory)
            
            # Create progress bar
            with tqdm(total=file_size, unit='B', unit_scale=True, desc=filename) as pbar:
                def callback(data):
                    f.write(data)
                    pbar.update(len(data))
                
                ftp.retrbinary(f"RETR {filename}", callback)
            
            ftp.quit()
            
        logger.info(f"Downloaded {filename} to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Error downloading file: {e}")
        return False

def extract_gzip(gzip_path, output_path):
    """
    Extract a gzip file to the specified output path.
    """
    try:
        with gzip.open(gzip_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                logger.info(f"Extracting {gzip_path} to {output_path}")
                shutil.copyfileobj(f_in, f_out)
        logger.info(f"Extraction complete")
        return True
    except Exception as e:
        logger.error(f"Error extracting file: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Download GRCh38 reference genome FASTA file")
    parser.add_argument('--output-dir', type=str, default='data/reference',
                        help='Directory to save the downloaded FASTA file')
    parser.add_argument('--extract', action='store_true', 
                        help='Extract the downloaded gzip file')
    parser.add_argument('--source', type=str, choices=['ensembl', 'ncbi'], default='ensembl',
                        help='Source to download from (ensembl or ncbi)')
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if args.source == 'ensembl':
        server = ENSEMBL_FTP_SERVER
        directory = ENSEMBL_FTP_DIR
        filename = ENSEMBL_FILENAME
    else:  # ncbi
        server = NCBI_FTP_SERVER
        directory = NCBI_FTP_DIR
        filename = NCBI_FILENAME
    
    # Download file
    gz_path = output_dir / filename
    success = download_file_with_progress(server, directory, filename, str(gz_path))
    
    if success and args.extract:
        # Extract the gzip file
        output_path = output_dir / filename.replace('.gz', '')
        extract_gzip(gz_path, output_path)
        
        # Optionally remove the gz file after extraction
        if os.path.exists(output_path):
            os.remove(gz_path)
            logger.info(f"Removed compressed file {gz_path}")

if __name__ == "__main__":
    main()
