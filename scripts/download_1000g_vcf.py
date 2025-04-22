#!/usr/bin/env python3
"""
Script to download VCF files from the 1000 Genomes Project.
Downloads up to 3 phased and 3 unphased VCF files.
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

# 1000 Genomes FTP server details
FTP_SERVER = "ftp.1000genomes.ebi.ac.uk"

# Define some specific files that we know are phased/unphased
KNOWN_PHASED_VCFS = [
    "/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
]

KNOWN_UNPHASED_VCFS = [
    # Using smaller chromosomes/regions for testing
    "/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom22.ILLUMINA.bwa.GBR.low_coverage.20130415.bam.genotypes.vcf.gz",
    "/vol1/ftp/phase3/data/HG00097/alignment/HG00097.chrom22.ILLUMINA.bwa.GBR.low_coverage.20130415.bam.genotypes.vcf.gz",
    "/vol1/ftp/phase3/data/HG00099/alignment/HG00099.chrom22.ILLUMINA.bwa.GBR.low_coverage.20130415.bam.genotypes.vcf.gz"
]

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

def extract_gzip(gzip_path, output_path, remove_gz=False):
    """
    Extract a gzip file to the specified output path.
    """
    try:
        with gzip.open(gzip_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                logger.info(f"Extracting {gzip_path}")
                shutil.copyfileobj(f_in, f_out)
        
        if remove_gz and os.path.exists(output_path):
            os.remove(gzip_path)
            logger.info(f"Removed compressed file {gzip_path}")
            
        return True
    except Exception as e:
        logger.error(f"Error extracting file: {e}")
        return False

def download_vcfs(output_dir, phased_count=3, unphased_count=3, extract=False):
    """
    Download specified number of phased and unphased VCF files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Download phased VCFs
    phased_vcfs = []
    for i in range(min(phased_count, len(KNOWN_PHASED_VCFS))):
        vcf_path = KNOWN_PHASED_VCFS[i]
        filename = os.path.basename(vcf_path)
        output_path = output_dir / f"phased_{i+1}_{filename}"
        
        url = f"ftp://{FTP_SERVER}{vcf_path}"
        success = download_file_with_progress(url, output_path, is_ftp=True)
        
        if success and extract:
            extract_path = output_path.with_suffix('').with_suffix('')  # Remove .vcf.gz
            extract_gzip(output_path, extract_path, remove_gz=False)
            phased_vcfs.append(extract_path)
        else:
            phased_vcfs.append(output_path)
    
    # Download unphased VCFs
    unphased_vcfs = []
    for i in range(min(unphased_count, len(KNOWN_UNPHASED_VCFS))):
        vcf_path = KNOWN_UNPHASED_VCFS[i]
        filename = os.path.basename(vcf_path)
        output_path = output_dir / f"unphased_{i+1}_{filename}"
        
        url = f"ftp://{FTP_SERVER}{vcf_path}"
        success = download_file_with_progress(url, output_path, is_ftp=True)
        
        if success and extract:
            extract_path = output_path.with_suffix('').with_suffix('')  # Remove .vcf.gz
            extract_gzip(output_path, extract_path, remove_gz=False)
            unphased_vcfs.append(extract_path)
        else:
            unphased_vcfs.append(output_path)
    
    return phased_vcfs, unphased_vcfs

def main():
    parser = argparse.ArgumentParser(description="Download VCF files from the 1000 Genomes Project")
    parser.add_argument('--output-dir', type=str, default='data/vcf',
                        help='Directory to save the downloaded VCF files')
    parser.add_argument('--phased', type=int, default=3,
                        help='Number of phased VCF files to download (max 3)')
    parser.add_argument('--unphased', type=int, default=3,
                        help='Number of unphased VCF files to download (max 3)')
    parser.add_argument('--extract', action='store_true', 
                        help='Extract the downloaded gzip files')
    args = parser.parse_args()
    
    phased_vcfs, unphased_vcfs = download_vcfs(
        args.output_dir, 
        phased_count=args.phased, 
        unphased_count=args.unphased,
        extract=args.extract
    )
    
    logger.info("Download summary:")
    logger.info(f"Phased VCFs ({len(phased_vcfs)}):")
    for vcf in phased_vcfs:
        logger.info(f"  - {vcf}")
    
    logger.info(f"Unphased VCFs ({len(unphased_vcfs)}):")
    for vcf in unphased_vcfs:
        logger.info(f"  - {vcf}")

if __name__ == "__main__":
    main()
