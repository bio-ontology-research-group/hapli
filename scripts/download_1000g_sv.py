#!/usr/bin/env python3
"""
Script to download structural variant (SV) VCF files from the 1000 Genomes Project.
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

# 1000 Genomes FTP server details
FTP_SERVER = "ftp.1000genomes.ebi.ac.uk"
FTP_SV_DIR = "/vol1/ftp/phase3/integrated_sv_map"
SV_FILENAME_PATTERN = "ALL.wgs.mergedSV.v8.20130502.svs.genotypes.{}.vcf.gz"

# Smaller chromosomes for faster download
DEFAULT_CHROMOSOMES = ['21', '22']  # Smaller chromosomes

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
            import requests
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

def download_sv_vcfs(output_dir, num_samples=3, chromosomes=None, extract=False):
    """
    Download structural variant VCF files from the 1000 Genomes Project.
    
    Args:
        output_dir: Directory to save downloaded files
        num_samples: Number of chromosomes to download (using different chromosomes as "samples")
        chromosomes: Specific chromosomes to download
        extract: Whether to extract the downloaded files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if chromosomes is None:
        # Use DEFAULT_CHROMOSOMES first, then add more if needed
        chromosomes = DEFAULT_CHROMOSOMES.copy()
        if num_samples > len(chromosomes):
            # Add more chromosomes if needed
            additional_chroms = ['20', '19', '18', '17', '16']
            chromosomes.extend(additional_chroms[:num_samples - len(chromosomes)])
    
    # Make sure we only process the requested number of samples
    chromosomes = chromosomes[:num_samples]
    
    downloaded_files = []
    
    for i, chrom in enumerate(chromosomes):
        sv_filename = SV_FILENAME_PATTERN.format(chrom)
        output_path = output_dir / f"sv_sample_{i+1}_chr{chrom}_{sv_filename}"
        
        # Build FTP URL
        url = f"ftp://{FTP_SERVER}{FTP_SV_DIR}/{sv_filename}"
        
        # Download file
        success = download_file_with_progress(url, output_path, is_ftp=True)
        
        if success:
            logger.info(f"Successfully downloaded SV file for chromosome {chrom}")
            if extract:
                extract_path = str(output_path).replace('.gz', '')
                if extract_gzip(output_path, extract_path, remove_gz=False):
                    downloaded_files.append(extract_path)
                else:
                    downloaded_files.append(output_path)
            else:
                downloaded_files.append(output_path)
        else:
            logger.error(f"Failed to download SV file for chromosome {chrom}")
    
    return downloaded_files

def main():
    parser = argparse.ArgumentParser(description="Download structural variant VCF files from the 1000 Genomes Project")
    parser.add_argument('--output-dir', type=str, default='data/structural_variants',
                        help='Directory to save the downloaded SV VCF files')
    parser.add_argument('--samples', type=int, default=3,
                        help='Number of SV files to download (using different chromosomes)')
    parser.add_argument('--chromosomes', type=str, nargs='+',
                        help='Specific chromosomes to download (overrides --samples)')
    parser.add_argument('--extract', action='store_true', 
                        help='Extract the downloaded gzip files')
    args = parser.parse_args()
    
    downloaded_files = download_sv_vcfs(
        args.output_dir,
        num_samples=args.samples,
        chromosomes=args.chromosomes,
        extract=args.extract
    )
    
    logger.info("Download summary:")
    logger.info(f"Structural Variant VCFs ({len(downloaded_files)}):")
    for vcf in downloaded_files:
        logger.info(f"  - {vcf}")

if __name__ == "__main__":
    main()
