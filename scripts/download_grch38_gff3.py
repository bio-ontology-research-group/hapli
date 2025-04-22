#!/usr/bin/env python3
"""
Script to download GFF3 annotation files matching the GRCh38 reference genome.
"""

import argparse
import os
import logging
import ftplib
import gzip
import shutil
from pathlib import Path
from tqdm import tqdm
import requests

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Define FTP server details for Ensembl
ENSEMBL_FTP_SERVER = "ftp.ensembl.org"
ENSEMBL_FTP_DIR = "/pub/release-109/gff3/homo_sapiens/"
ENSEMBL_FILENAME = "Homo_sapiens.GRCh38.109.gff3.gz"
ENSEMBL_CHR_DIR = "/pub/release-109/gff3/homo_sapiens/"
ENSEMBL_CHR_PATTERN = "Homo_sapiens.GRCh38.109.chromosome.{}.gff3.gz"  # Format with chromosome number

# Define FTP server details for NCBI
NCBI_FTP_SERVER = "ftp.ncbi.nlm.nih.gov"
NCBI_FTP_DIR = "/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.gff.gz"
NCBI_FILENAME = "GCA_000001405.15_GRCh38_genomic.gff.gz"

# Define HTTP URL for Gencode
GENCODE_HTTP_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gff3.gz"
GENCODE_FILENAME = "gencode.v44.annotation.gff3.gz"

def download_file_with_progress(url, output_path, is_ftp=True):
    """
    Download a file with a progress bar.
    Works with both FTP and HTTP/HTTPS URLs.
    """
    try:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        
        if is_ftp:
            # Parse FTP URL
            server = url.split('/')[2]
            path = '/' + '/'.join(url.split('/')[3:])
            directory = os.path.dirname(path)
            filename = os.path.basename(path)
            
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
        else:
            # HTTP/HTTPS download
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            # Get file size for progress bar
            file_size = int(response.headers.get('content-length', 0))
            
            # Download with progress bar
            with open(output_path, 'wb') as f:
                with tqdm(total=file_size, unit='B', unit_scale=True, desc=os.path.basename(url)) as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
        
        logger.info(f"Downloaded to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Error downloading file: {e}")
        return False

def extract_gzip(gzip_path, output_path=None, remove_gz=False):
    """
    Extract a gzip file to the specified output path.
    
    Args:
        gzip_path: Path to the gzipped file
        output_path: Path to extract to (if None, removes only the .gz extension)
        remove_gz: Whether to remove the original gzipped file after extraction
    """
    try:
        # If output_path is not specified, remove only the .gz extension
        if output_path is None:
            output_path = str(gzip_path).replace('.gz', '')
        
        with gzip.open(gzip_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                logger.info(f"Extracting {gzip_path} to {output_path}")
                shutil.copyfileobj(f_in, f_out)
        
        if remove_gz and os.path.exists(output_path):
            os.remove(gzip_path)
            logger.info(f"Removed compressed file {gzip_path}")
            
        return True
    except Exception as e:
        logger.error(f"Error extracting file: {e}")
        return False

def download_chromosome_specific_gff(output_dir, source='ensembl', chromosomes=None, extract=False):
    """
    Download chromosome-specific GFF files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if chromosomes is None:
        # Default to a few common chromosomes
        chromosomes = ['1', '2', '3', '21', '22', 'X', 'Y']
    
    downloaded_files = []
    
    if source == 'ensembl':
        server = ENSEMBL_FTP_SERVER
        directory = ENSEMBL_CHR_DIR
        
        for chrom in chromosomes:
            filename = ENSEMBL_CHR_PATTERN.format(chrom)
            output_path = output_dir / f"chr{chrom}_{filename}"
            
            # Check if file exists (for resumption)
            if os.path.exists(output_path) or os.path.exists(output_path.with_suffix('').with_suffix('')):
                logger.info(f"File already exists: {output_path}")
                downloaded_files.append(output_path)
                continue
                
            url = f"ftp://{server}{directory}/{filename}"
            success = download_file_with_progress(url, output_path, is_ftp=True)
            
            if success:
                downloaded_files.append(output_path)
                if extract:
                    # Only remove .gz suffix, keeping the .gff3 extension
                    extract_path = output_path.with_suffix('')
                    extract_gzip(output_path, extract_path, remove_gz=False)
    
    return downloaded_files

def download_gff3(output_dir, source='ensembl', extract=False):
    """
    Download GFF3 annotation file from the specified source.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if source == 'ensembl':
        server = ENSEMBL_FTP_SERVER
        directory = ENSEMBL_FTP_DIR
        filename = ENSEMBL_FILENAME
        url = f"ftp://{server}{directory}/{filename}"
        is_ftp = True
    elif source == 'ncbi':
        url = f"ftp://{NCBI_FTP_SERVER}{NCBI_FTP_DIR}"
        filename = NCBI_FILENAME
        is_ftp = True
    elif source == 'gencode':
        url = GENCODE_HTTP_URL
        filename = GENCODE_FILENAME
        is_ftp = False
    else:
        logger.error(f"Unknown source: {source}")
        return None
    
    output_path = output_dir / filename
    
    # Check if file exists (for resumption)
    if os.path.exists(output_path):
        logger.info(f"File already exists: {output_path}")
        return output_path
    
    success = download_file_with_progress(url, output_path, is_ftp=is_ftp)
    
    if success and extract:
        # Only remove .gz suffix, keeping the .gff3 extension
        extract_path = output_path.with_suffix('')  # Remove just .gz
        extract_gzip(output_path, extract_path, remove_gz=False)
        return extract_path
    
    return output_path if success else None

def main():
    parser = argparse.ArgumentParser(description="Download GFF3 annotation files for GRCh38")
    parser.add_argument('--output-dir', type=str, default='data/annotation',
                        help='Directory to save the downloaded files')
    parser.add_argument('--extract', action='store_true', 
                        help='Extract the downloaded gzip files')
    parser.add_argument('--source', type=str, choices=['ensembl', 'ncbi', 'gencode'], default='ensembl',
                        help='Source to download from (ensembl, ncbi, or gencode)')
    parser.add_argument('--by-chromosome', action='store_true',
                        help='Download chromosome-specific files (Ensembl only)')
    parser.add_argument('--chromosomes', type=str, nargs='+',
                        help='Specific chromosomes to download (e.g., 1 2 X Y)')
    
    args = parser.parse_args()
    
    if args.by_chromosome and args.source == 'ensembl':
        downloaded_files = download_chromosome_specific_gff(
            args.output_dir,
            source=args.source,
            chromosomes=args.chromosomes,
            extract=args.extract
        )
        
        logger.info(f"Downloaded {len(downloaded_files)} chromosome-specific GFF files:")
        for f in downloaded_files:
            logger.info(f"  - {f}")
    else:
        output_path = download_gff3(
            args.output_dir,
            source=args.source,
            extract=args.extract
        )
        
        if output_path:
            logger.info(f"GFF3 file downloaded successfully to {output_path}")
        else:
            logger.error("Failed to download GFF3 file")
            return 1
    
    return 0

if __name__ == "__main__":
    main()
