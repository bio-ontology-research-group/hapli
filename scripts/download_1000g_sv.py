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
SV_FILENAME = "ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"

# No longer using chromosomes

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
            except ftplib.error_perm as e:
                logger.error(f"FTP error: {e}")
                # Try alternative pattern as fallback
                if "mergedSV" not in filename:
                    alt_filename = filename.replace("integrated_sv_map_v2", "mergedSV.v8")
                    logger.info(f"Trying alternative filename: {alt_filename}")
                    try:
                        ftp = ftplib.FTP(server)
                        ftp.login()
                        ftp.cwd(directory)
                        file_size = ftp.size(alt_filename)
                        ftp.quit()
                        
                        with open(output_path, 'wb') as f:
                            logger.info(f"Connecting to {server}")
                            ftp = ftplib.FTP(server)
                            ftp.login()
                            ftp.cwd(directory)
                            
                            # Create progress bar
                            with tqdm(total=file_size, unit='B', unit_scale=True, desc=alt_filename) as pbar:
                                def callback(data):
                                    f.write(data)
                                    pbar.update(len(data))
                                
                                ftp.retrbinary(f"RETR {alt_filename}", callback)
                            
                            ftp.quit()
                    except Exception as e2:
                        logger.error(f"Error with alternative filename: {e2}")
                        raise
                else:
                    raise
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
        num_samples: Number of samples to download (multiple copies with different names)
        chromosomes: Not used - kept for backwards compatibility
        extract: Whether to extract the downloaded files
    """
    # Alternative locations to try if main location fails
    ALTERNATIVE_SV_DIRS = [
        "/vol1/ftp/phase3/integrated_sv_map",
        "/vol1/ftp/phase3/integrated_sv_map/ALL",
        "/vol1/ftp/release/20130502/supporting/sv_deletions"
    ]
    
    # Alternative filenames to try
    ALTERNATIVE_SV_FILES = [
        "ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz",
        "ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz"
    ]
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    downloaded_files = []
    found = False
    
    # Try the primary location and filename first
    for i in range(num_samples):
        # Use sample_index for naming multiple copies
        sample_index = i + 1
        output_path = output_dir / f"sv_sample_{sample_index}_{SV_FILENAME}"
        
        # Only download once but create multiple symlinks or copies
        if i == 0 or not found:
            # Build FTP URL
            url = f"ftp://{FTP_SERVER}{FTP_SV_DIR}/{SV_FILENAME}"
            logger.info(f"Trying to download from {url}")
            
            # Download file
            success = download_file_with_progress(url, output_path, is_ftp=True)
            
            if success:
                found = True
                logger.info(f"Successfully downloaded SV file")
                if extract:
                    extract_path = str(output_path).replace('.gz', '')
                    if extract_gzip(output_path, extract_path, remove_gz=False):
                        downloaded_files.append(extract_path)
                    else:
                        downloaded_files.append(output_path)
                else:
                    downloaded_files.append(output_path)
        else:
            # For additional samples, create a copy (or symlink) of the first file
            if i > 0 and found:
                try:
                    # Use shutil.copy2 to preserve metadata
                    first_file = downloaded_files[0]
                    if extract and first_file.endswith('.vcf'):
                        # Copy the extracted file
                        extract_path = str(output_path).replace('.gz', '')
                        shutil.copy2(first_file, extract_path)
                        downloaded_files.append(extract_path)
                        # Also copy the compressed file
                        shutil.copy2(str(first_file) + '.gz', output_path)
                    else:
                        # Copy the compressed file
                        shutil.copy2(first_file, output_path)
                        downloaded_files.append(output_path)
                        # If extract is True, also extract this copy
                        if extract:
                            extract_path = str(output_path).replace('.gz', '')
                            if extract_gzip(output_path, extract_path, remove_gz=False):
                                downloaded_files.append(extract_path)
                
                except Exception as e:
                    logger.error(f"Error creating additional sample copy: {e}")
    
    # If primary location fails, try alternatives
    if not found:
        logger.info("Trying alternative locations")
        
        # Try alternative directories and filenames
        for alt_dir in ALTERNATIVE_SV_DIRS:
            for alt_filename in ALTERNATIVE_SV_FILES:
                if found:
                    break
                    
                alt_output_path = output_dir / f"sv_sample_1_{alt_filename}"
                
                alt_url = f"ftp://{FTP_SERVER}{alt_dir}/{alt_filename}"
                logger.info(f"Trying alternative URL: {alt_url}")
                
                alt_success = download_file_with_progress(alt_url, alt_output_path, is_ftp=True)
                
                if alt_success:
                    found = True
                    logger.info(f"Successfully downloaded SV file from alternative location")
                    if extract:
                        alt_extract_path = str(alt_output_path).replace('.gz', '')
                        if extract_gzip(alt_output_path, alt_extract_path, remove_gz=False):
                            downloaded_files.append(alt_extract_path)
                        else:
                            downloaded_files.append(alt_output_path)
                    else:
                        downloaded_files.append(alt_output_path)
        
        if not found:
            logger.error(f"Failed to download SV file after trying all alternatives")
    
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
