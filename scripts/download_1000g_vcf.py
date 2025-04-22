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
import subprocess
import random
from pathlib import Path
from tqdm import tqdm
import requests

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 1000 Genomes FTP server details
FTP_SERVER = "ftp.1000genomes.ebi.ac.uk"

def check_bcftools_available():
    """Check if bcftools is available in the path."""
    try:
        result = subprocess.run(["bcftools", "--version"], 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE, 
                                text=True)
        if result.returncode == 0:
            logger.info("bcftools is available: " + result.stdout.splitlines()[0])
            return True
        else:
            logger.warning("bcftools command failed")
            return False
    except FileNotFoundError:
        logger.warning("bcftools not found in PATH")
        return False

def get_samples_from_vcf(vcf_path):
    """Extract list of sample IDs from a VCF file using bcftools."""
    try:
        result = subprocess.run(
            ["bcftools", "query", "-l", vcf_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        if result.returncode == 0:
            samples = result.stdout.strip().split('\n')
            logger.info(f"Found {len(samples)} samples in VCF file")
            return samples
        else:
            logger.error(f"Error getting samples: {result.stderr}")
            return []
    except Exception as e:
        logger.error(f"Error executing bcftools: {e}")
        return []

def extract_sample_from_vcf(input_vcf, output_vcf, sample_id):
    """Extract a single sample from a VCF file using bcftools."""
    try:
        logger.info(f"Extracting sample {sample_id} from {input_vcf} to {output_vcf}")
        result = subprocess.run(
            ["bcftools", "view", "-s", sample_id, "-Ov", "-o", output_vcf, input_vcf],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        if result.returncode == 0:
            logger.info(f"Successfully extracted sample {sample_id}")
            return True
        else:
            logger.error(f"Error extracting sample: {result.stderr}")
            return False
    except Exception as e:
        logger.error(f"Error executing bcftools: {e}")
        return False

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

def download_vcfs(output_dir, phased_count=3, unphased_count=3, extract=False, force_unphased=False, specific_samples=None):
    """
    Download specified number of phased and unphased VCF files.
    
    Args:
        output_dir: Directory to save downloads
        phased_count: Number of phased VCFs to download
        unphased_count: Number of unphased VCFs to download
        extract: Whether to extract the downloaded files
        force_unphased: Force download of unphased VCFs even if some fail
        specific_samples: List of specific sample IDs to extract
    """
    # Check if bcftools is available
    bcftools_available = check_bcftools_available()
    if not bcftools_available and specific_samples:
        logger.warning("bcftools not available - will download multi-sample VCF files but cannot extract individual samples")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Download phased VCFs
    phased_vcfs = []
    
    # We'll only download the first phased VCF and then extract specific samples from it
    if phased_count > 0:
        vcf_path = KNOWN_PHASED_VCFS[0]
        filename = os.path.basename(vcf_path)
        
        # Download the main multi-sample phased VCF
        output_path = output_dir / f"phased_full_{filename}"
        main_extracted_path = None
        
        url = f"ftp://{FTP_SERVER}{vcf_path}"
        success = download_file_with_progress(url, output_path, is_ftp=True)
        
        if success:
            logger.info(f"Successfully downloaded phased VCF file")
            
            # Extract the gzipped file if requested
            if extract:
                main_extracted_path = str(output_path).replace('.gz', '')
                if extract_gzip(output_path, main_extracted_path, remove_gz=False):
                    logger.info(f"Extracted to {main_extracted_path}")
                else:
                    logger.error("Failed to extract phased VCF file")
                    main_extracted_path = None
            
            # Path to use for bcftools
            vcf_for_samples = main_extracted_path if main_extracted_path else str(output_path)
            
            # If bcftools is available and sample IDs are specified, extract individual samples
            if bcftools_available and specific_samples:
                # Get list of samples from the VCF
                all_samples = get_samples_from_vcf(vcf_for_samples)
                
                if all_samples:
                    logger.info(f"Phased VCF contains {len(all_samples)} samples")
                    
                    # Determine which samples to extract
                    samples_to_extract = []
                    if specific_samples:
                        # Use specified sample IDs if provided
                        samples_to_extract = [s for s in specific_samples if s in all_samples]
                        if len(samples_to_extract) < len(specific_samples):
                            logger.warning(f"Some specified samples not found in phased VCF file")
                            # Add random samples to make up the number
                            missing = min(phased_count, len(specific_samples)) - len(samples_to_extract)
                            remaining_samples = [s for s in all_samples if s not in samples_to_extract]
                            if remaining_samples and missing > 0:
                                extra_samples = random.sample(remaining_samples, min(missing, len(remaining_samples)))
                                logger.info(f"Adding random samples to meet count: {', '.join(extra_samples)}")
                                samples_to_extract.extend(extra_samples)
                    else:
                        # Select random samples
                        samples_to_extract = random.sample(all_samples, min(phased_count, len(all_samples)))
                    
                    logger.info(f"Extracting {len(samples_to_extract)} samples from phased VCF: {', '.join(samples_to_extract)}")
                    
                    # Extract each sample
                    for i, sample_id in enumerate(samples_to_extract):
                        sample_output = output_dir / f"phased_{i+1}_{sample_id}.vcf"
                        if extract_sample_from_vcf(vcf_for_samples, sample_output, sample_id):
                            phased_vcfs.append(str(sample_output))
                else:
                    logger.error("Could not extract sample information from phased VCF")
                    # Add the full VCF as fallback
                    if main_extracted_path:
                        phased_vcfs.append(main_extracted_path)
                    else:
                        phased_vcfs.append(str(output_path))
            else:
                # Without bcftools or specific samples, download multiple phased VCFs directly
                if specific_samples:
                    logger.warning("Cannot extract specific samples without bcftools. Using full VCF files.")
                
                # Add the first VCF we already downloaded
                if main_extracted_path:
                    phased_vcfs.append(main_extracted_path)
                else:
                    phased_vcfs.append(str(output_path))
                
                # Download additional phased VCFs if needed
                for i in range(1, min(phased_count, len(KNOWN_PHASED_VCFS))):
                    vcf_path = KNOWN_PHASED_VCFS[i]
                    filename = os.path.basename(vcf_path)
                    add_output_path = output_dir / f"phased_{i+1}_{filename}"
                    
                    url = f"ftp://{FTP_SERVER}{vcf_path}"
                    success = download_file_with_progress(url, add_output_path, is_ftp=True)
                    
                    if success and extract:
                        extract_path = add_output_path.with_suffix('').with_suffix('')  # Remove .vcf.gz
                        extract_gzip(add_output_path, extract_path, remove_gz=False)
                        phased_vcfs.append(extract_path)
                    elif success:
                        phased_vcfs.append(add_output_path)
    
    # Alternative unphased VCFs - more likely to download successfully
    ALTERNATIVE_UNPHASED_VCFS = [
        "/vol1/ftp/phase3/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam.genotypes.vcf.gz",
        "/vol1/ftp/phase3/data/NA12878/alignment/NA12878.chrom21.ILLUMINA.bwa.CEU.low_coverage.20121211.bam.genotypes.vcf.gz",
        "/vol1/ftp/phase3/data/NA12878/alignment/NA12878.chrom22.ILLUMINA.bwa.CEU.low_coverage.20121211.bam.genotypes.vcf.gz"
    ]
    
    # Download unphased VCFs
    unphased_vcfs = []
    
    # If we have specific samples and they're in our alternative list, use those files directly
    sample_to_vcf_map = {
        "NA12878": [ALTERNATIVE_UNPHASED_VCFS[0], ALTERNATIVE_UNPHASED_VCFS[1], ALTERNATIVE_UNPHASED_VCFS[2]]
    }
    
    sample_specific_downloads = False
    if specific_samples and any(s in sample_to_vcf_map for s in specific_samples):
        logger.info("Using sample-specific VCF sources for unphased VCFs")
        sample_specific_downloads = True
        
        sample_idx = 0
        for sample in specific_samples:
            if sample in sample_to_vcf_map and sample_idx < unphased_count:
                for vcf_path in sample_to_vcf_map[sample]:
                    filename = os.path.basename(vcf_path)
                    output_path = output_dir / f"unphased_{sample_idx+1}_{sample}_{filename}"
                    
                    url = f"ftp://{FTP_SERVER}{vcf_path}"
                    success = download_file_with_progress(url, output_path, is_ftp=True)
                    
                    if success:
                        if extract:
                            extract_path = output_path.with_suffix('').with_suffix('')  # Remove .vcf.gz
                            extract_gzip(output_path, extract_path, remove_gz=False)
                            unphased_vcfs.append(extract_path)
                        else:
                            unphased_vcfs.append(output_path)
                        
                        sample_idx += 1
                        if sample_idx >= unphased_count:
                            break
    
    # If we don't have specific samples or haven't downloaded enough samples yet, use the standard approach
    unphased_count_remaining = unphased_count - len(unphased_vcfs)
    
    if not sample_specific_downloads and unphased_count_remaining > 0:
        # Try the first source - for a multi-sample VCF if possible
        vcf_path = KNOWN_UNPHASED_VCFS[0]
        filename = os.path.basename(vcf_path)
        main_output_path = output_dir / f"unphased_full_{filename}"
        
        url = f"ftp://{FTP_SERVER}{vcf_path}"
        success = download_file_with_progress(url, main_output_path, is_ftp=True)
        
        if success:
            logger.info(f"Successfully downloaded unphased VCF file")
            main_extracted_path = None
            
            # Extract the gzipped file if requested
            if extract:
                main_extracted_path = str(main_output_path).replace('.gz', '')
                if extract_gzip(main_output_path, main_extracted_path, remove_gz=False):
                    logger.info(f"Extracted to {main_extracted_path}")
                else:
                    logger.error("Failed to extract unphased VCF file")
                    main_extracted_path = None
            
            # Path to use for bcftools
            vcf_for_samples = main_extracted_path if main_extracted_path else str(main_output_path)
            
            # If bcftools is available and sample IDs are specified, extract individual samples
            if bcftools_available and specific_samples:
                # Get list of samples from the VCF
                all_samples = get_samples_from_vcf(vcf_for_samples)
                
                if all_samples:
                    logger.info(f"Unphased VCF contains {len(all_samples)} samples")
                    
                    # Determine which samples to extract
                    samples_to_extract = []
                    if specific_samples:
                        # Use specified sample IDs if provided
                        samples_to_extract = [s for s in specific_samples if s in all_samples]
                        if len(samples_to_extract) < len(specific_samples):
                            logger.warning(f"Some specified samples not found in unphased VCF file")
                            # Add random samples to make up the number
                            missing = min(unphased_count, len(specific_samples)) - len(samples_to_extract)
                            remaining_samples = [s for s in all_samples if s not in samples_to_extract]
                            if remaining_samples and missing > 0:
                                extra_samples = random.sample(remaining_samples, min(missing, len(remaining_samples)))
                                logger.info(f"Adding random samples to meet count: {', '.join(extra_samples)}")
                                samples_to_extract.extend(extra_samples)
                    else:
                        # Select random samples
                        samples_to_extract = random.sample(all_samples, min(unphased_count, len(all_samples)))
                    
                    logger.info(f"Extracting {len(samples_to_extract)} samples from unphased VCF: {', '.join(samples_to_extract)}")
                    
                    # Extract each sample
                    for i, sample_id in enumerate(samples_to_extract):
                        sample_output = output_dir / f"unphased_{i+1}_{sample_id}.vcf"
                        if extract_sample_from_vcf(vcf_for_samples, sample_output, sample_id):
                            unphased_vcfs.append(str(sample_output))
                else:
                    logger.error("Could not extract sample information from unphased VCF")
                    # Add the full VCF as fallback
                    if main_extracted_path:
                        unphased_vcfs.append(main_extracted_path)
                    else:
                        unphased_vcfs.append(str(main_output_path))
            else:
                # Without bcftools or specific samples, download individual VCFs directly
                if specific_samples:
                    logger.warning("Cannot extract specific samples without bcftools. Using full VCF files.")
                
                # Add the first VCF we already downloaded
                if main_extracted_path:
                    unphased_vcfs.append(main_extracted_path)
                else:
                    unphased_vcfs.append(str(main_output_path))
                
                # Download additional unphased VCFs if needed
                unphased_count_remaining = unphased_count - 1
        
        # If we need more unphased VCFs, try additional sources
        if unphased_count_remaining > 0:
            # First try the original sources (after the first one)
            for i in range(1, min(unphased_count_remaining + 1, len(KNOWN_UNPHASED_VCFS))):
                vcf_path = KNOWN_UNPHASED_VCFS[i]
                filename = os.path.basename(vcf_path)
                output_path = output_dir / f"unphased_{len(unphased_vcfs)+1}_{filename}"
                
                url = f"ftp://{FTP_SERVER}{vcf_path}"
                success = download_file_with_progress(url, output_path, is_ftp=True)
                
                if success:
                    if extract:
                        extract_path = output_path.with_suffix('').with_suffix('')  # Remove .vcf.gz
                        extract_gzip(output_path, extract_path, remove_gz=False)
                        unphased_vcfs.append(extract_path)
                    else:
                        unphased_vcfs.append(output_path)
                    unphased_count_remaining -= 1
                    
                    if unphased_count_remaining <= 0:
                        break
            
            # If we still need more and force_unphased is enabled, try alternative sources
            if unphased_count_remaining > 0 and force_unphased:
                logger.info("Using alternative unphased VCF sources")
                for i in range(min(unphased_count_remaining, len(ALTERNATIVE_UNPHASED_VCFS))):
                    vcf_path = ALTERNATIVE_UNPHASED_VCFS[i]
                    filename = os.path.basename(vcf_path)
                    output_path = output_dir / f"unphased_alt_{len(unphased_vcfs)+1}_{filename}"
                    
                    url = f"ftp://{FTP_SERVER}{vcf_path}"
                    success = download_file_with_progress(url, output_path, is_ftp=True)
                    
                    if success:
                        if extract:
                            extract_path = output_path.with_suffix('').with_suffix('')  # Remove .vcf.gz
                            extract_gzip(output_path, extract_path, remove_gz=False)
                            unphased_vcfs.append(extract_path)
                        else:
                            unphased_vcfs.append(output_path)
                        unphased_count_remaining -= 1
                        
                        if unphased_count_remaining <= 0:
                            break
    
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
    parser.add_argument('--force-unphased', action='store_true',
                        help='Force download of unphased VCFs using alternative sources if needed')
    parser.add_argument('--sample-ids', type=str, nargs='+',
                        help='Specific sample IDs to extract (e.g., NA12878 HG00096)')
    args = parser.parse_args()
    
    phased_vcfs, unphased_vcfs = download_vcfs(
        args.output_dir, 
        phased_count=args.phased, 
        unphased_count=args.unphased,
        extract=args.extract,
        force_unphased=args.force_unphased,
        specific_samples=args.sample_ids
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
