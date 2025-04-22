#!/usr/bin/env python3
"""
Script to download VCF files from the 1000 Genomes Project.
Downloads specified VCFs, optionally extracts specific samples into separate
compressed VCF files (.vcf.gz).
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
from typing import List, Tuple, Optional

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
                                check=True, text=True)
        if result.returncode == 0:
            logger.info("bcftools is available: " + result.stdout.splitlines()[0])
            return True
        else:
            # This path might not be reached due to check=True
            logger.warning("bcftools command failed")
            return False
    except FileNotFoundError:
        logger.warning("bcftools not found in PATH")
        return False
    except subprocess.CalledProcessError as e:
        logger.warning(f"bcftools check failed: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"Error checking for bcftools: {e}")
        return False


def get_samples_from_vcf(vcf_path: str) -> List[str]:
    """Extract list of sample IDs from a VCF file (can be .vcf or .vcf.gz) using bcftools."""
    try:
        # bcftools query -l works on both compressed and uncompressed VCFs
        cmd = ["bcftools", "query", "-l", vcf_path]
        logger.debug(f"Running command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True, text=True
        )
        samples = result.stdout.strip().split('\n')
        # Filter out empty strings if any
        samples = [s for s in samples if s]
        logger.info(f"Found {len(samples)} samples in VCF file: {vcf_path}")
        return samples
    except subprocess.CalledProcessError as e:
        logger.error(f"Error getting samples from {vcf_path}: {e.stderr}")
        return []
    except FileNotFoundError:
        logger.error(f"bcftools not found when trying to get samples from {vcf_path}")
        return []
    except Exception as e:
        logger.error(f"Error executing bcftools query: {e}")
        return []

def extract_sample_from_vcf(input_vcf: str, output_vcf_gz: str, sample_id: str) -> bool:
    """
    Extract a single sample from a VCF file using bcftools.
    Outputs a compressed VCF file (.vcf.gz).

    Args:
        input_vcf: Path to the input VCF file (can be .vcf or .vcf.gz).
        output_vcf_gz: Path to the output compressed VCF file (.vcf.gz).
        sample_id: The sample ID to extract.

    Returns:
        True if extraction was successful, False otherwise.
    """
    try:
        # Use -Oz to output compressed VCF
        cmd = [
            "bcftools", "view",
            "-s", sample_id,      # Select sample
            "-Oz",                # Output compressed VCF (gzipped)
            "-o", output_vcf_gz,  # Output file path
            input_vcf             # Input file path
        ]
        logger.info(f"Extracting sample {sample_id} from {Path(input_vcf).name} to {Path(output_vcf_gz).name}")
        logger.debug(f"Running command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE, # Capture stdout/stderr even on success
            stderr=subprocess.PIPE,
            check=True, text=True
        )
        # Log stderr for potential warnings from bcftools
        if result.stderr:
            logger.debug(f"bcftools view stderr for {sample_id}: {result.stderr.strip()}")

        # Verify output file was created
        if not Path(output_vcf_gz).exists() or Path(output_vcf_gz).stat().st_size == 0:
             logger.error(f"bcftools view command succeeded but output file {output_vcf_gz} is missing or empty.")
             return False

        logger.info(f"Successfully extracted sample {sample_id} to {output_vcf_gz}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error extracting sample {sample_id}: {e.stderr}")
        # Clean up potentially empty/corrupt output file
        if Path(output_vcf_gz).exists():
            try:
                os.remove(output_vcf_gz)
            except OSError:
                pass
        return False
    except FileNotFoundError:
        logger.error(f"bcftools not found when trying to extract sample {sample_id}")
        return False
    except Exception as e:
        logger.error(f"Error executing bcftools view for sample {sample_id}: {e}")
        return False


# Define some specific files that we know are phased/unphased
KNOWN_PHASED_VCFS = [
    "/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
]

# Common trio samples with reliable data availability
KNOWN_SAMPLE_IDS = ["NA12878", "NA12891", "NA12892"]

# Sample-specific unphased VCFs mapped by sample ID - more reliable than general unphased VCFs
SAMPLE_TO_UNPHASED_VCF = {
    "NA12878": [
        "/vol1/ftp/phase3/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam.genotypes.vcf.gz",
        "/vol1/ftp/phase3/data/NA12878/alignment/NA12878.chrom21.ILLUMINA.bwa.CEU.low_coverage.20121211.bam.genotypes.vcf.gz",
        "/vol1/ftp/phase3/data/NA12878/alignment/NA12878.chrom22.ILLUMINA.bwa.CEU.low_coverage.20121211.bam.genotypes.vcf.gz"
    ],
    "NA12891": [
        "/vol1/ftp/phase3/data/NA12891/alignment/NA12891.chrom20.ILLUMINA.bwa.CEU.low_coverage.20130415.bam.genotypes.vcf.gz",
        "/vol1/ftp/phase3/data/NA12891/alignment/NA12891.chrom21.ILLUMINA.bwa.CEU.low_coverage.20130415.bam.genotypes.vcf.gz",
        "/vol1/ftp/phase3/data/NA12891/alignment/NA12891.chrom22.ILLUMINA.bwa.CEU.low_coverage.20130415.bam.genotypes.vcf.gz"
    ],
    "NA12892": [
        "/vol1/ftp/phase3/data/NA12892/alignment/NA12892.chrom20.ILLUMINA.bwa.CEU.low_coverage.20130415.bam.genotypes.vcf.gz",
        "/vol1/ftp/phase3/data/NA12892/alignment/NA12892.chrom21.ILLUMINA.bwa.CEU.low_coverage.20130415.bam.genotypes.vcf.gz",
        "/vol1/ftp/phase3/data/NA12892/alignment/NA12892.chrom22.ILLUMINA.bwa.CEU.low_coverage.20130415.bam.genotypes.vcf.gz"
    ],
    # Add more samples if needed, e.g., HG01383 (from workflow example)
    "HG01383": [ # Assuming similar paths exist for this sample
        "/vol1/ftp/phase3/data/HG01383/alignment/HG01383.chrom20.ILLUMINA.bwa.PUR.low_coverage.20130415.bam.genotypes.vcf.gz",
        "/vol1/ftp/phase3/data/HG01383/alignment/HG01383.chrom21.ILLUMINA.bwa.PUR.low_coverage.20130415.bam.genotypes.vcf.gz",
        "/vol1/ftp/phase3/data/HG01383/alignment/HG01383.chrom22.ILLUMINA.bwa.PUR.low_coverage.20130415.bam.genotypes.vcf.gz",
    ]
}

# Fallback general unphased VCFs (multi-sample) - Use release files which are more standard
KNOWN_UNPHASED_VCFS_FALLBACK = [
    "/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", # v5a is often considered unphased relative to v5b
    "/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
]


def download_file_with_progress(url: str, output_path: str, is_ftp: bool = True) -> bool:
    """
    Download a file with a progress bar.
    Works with both FTP and HTTP/HTTPS URLs.
    """
    try:
        output_path_obj = Path(output_path)
        output_path_obj.parent.mkdir(parents=True, exist_ok=True)

        if is_ftp:
            # Parse FTP URL
            parts = url.split('/')
            server = parts[2]
            path = '/' + '/'.join(parts[3:])
            directory = os.path.dirname(path)
            filename = os.path.basename(path)

            logger.info(f"Connecting to FTP server: {server}")
            with ftplib.FTP(server) as ftp:
                ftp.login()
                ftp.cwd(directory)
                file_size = ftp.size(filename)

                if file_size is None:
                    logger.warning(f"Could not determine size of FTP file: {filename}")
                    file_size = 0 # Set to 0 to allow download attempt

                logger.info(f"Downloading FTP file: {filename} ({file_size} bytes)")
                with open(output_path_obj, 'wb') as f:
                    with tqdm(total=file_size, unit='B', unit_scale=True, desc=filename, disable=None) as pbar:
                        total_written = 0
                        def callback(data):
                            nonlocal total_written
                            f.write(data)
                            pbar.update(len(data))
                            total_written += len(data)

                        ftp.retrbinary(f"RETR {filename}", callback)

                # Verify downloaded size if possible
                if file_size > 0 and output_path_obj.stat().st_size != file_size:
                     logger.warning(f"Downloaded size ({output_path_obj.stat().st_size}) does not match expected size ({file_size}) for {filename}")

        else: # HTTP/HTTPS
            logger.info(f"Downloading HTTP(S) file: {os.path.basename(url)}")
            response = requests.get(url, stream=True, timeout=60) # Added timeout
            response.raise_for_status()

            file_size = int(response.headers.get('content-length', 0))

            with open(output_path_obj, 'wb') as f:
                with tqdm(total=file_size, unit='B', unit_scale=True, desc=os.path.basename(url), disable=None) as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))

            # Verify downloaded size if possible
            if file_size > 0 and output_path_obj.stat().st_size != file_size:
                 logger.warning(f"Downloaded size ({output_path_obj.stat().st_size}) does not match expected size ({file_size}) for {os.path.basename(url)}")


        logger.info(f"Successfully downloaded to {output_path}")
        return True
    except ftplib.all_errors as e:
        logger.error(f"FTP Error downloading {url}: {e}")
        return False
    except requests.exceptions.RequestException as e:
        logger.error(f"HTTP Error downloading {url}: {e}")
        return False
    except Exception as e:
        logger.error(f"Error downloading file {url}: {e}", exc_info=True)
        # Clean up partial download
        if Path(output_path).exists():
            try:
                os.remove(output_path)
            except OSError:
                pass
        return False

def extract_gzip(gzip_path: str, output_path: str, remove_gz: bool = False) -> bool:
    """
    Extract a gzip file to the specified output path.
    """
    try:
        logger.info(f"Extracting {gzip_path} to {output_path}")
        with gzip.open(gzip_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Verify extraction by checking output file existence and size
        if not Path(output_path).exists() or Path(output_path).stat().st_size == 0:
             # Check if input was empty
             if Path(gzip_path).stat().st_size > 0:
                 logger.error(f"Extraction failed: Output file {output_path} is missing or empty.")
                 return False
             else:
                 logger.warning(f"Input file {gzip_path} was empty, output {output_path} is also empty.")


        if remove_gz:
            try:
                os.remove(gzip_path)
                logger.info(f"Removed compressed file: {gzip_path}")
            except OSError as e:
                logger.warning(f"Could not remove compressed file {gzip_path}: {e}")

        return True
    except gzip.BadGzipFile:
        logger.error(f"Error extracting file: {gzip_path} is not a valid gzip file.")
        return False
    except Exception as e:
        logger.error(f"Error extracting file {gzip_path}: {e}", exc_info=True)
        return False

def download_vcfs(output_dir: str, phased_count: int = 1, unphased_count: int = 1, extract: bool = False, force_unphased: bool = False, specific_samples: Optional[List[str]] = None) -> Tuple[List[str], List[str]]:
    """
    Download specified number of phased and unphased VCF files.
    If specific_samples are provided and bcftools is available, extracts them into
    separate .vcf.gz files. Otherwise, downloads multi-sample files.

    Args:
        output_dir: Directory to save downloads.
        phased_count: Max number of phased VCFs (either multi-sample or per specified sample).
        unphased_count: Max number of unphased VCFs (either multi-sample or per specified sample).
        extract: If True, extract downloaded .gz files (DEPRECATED - use bcftools for sample extraction).
        force_unphased: Force download of fallback unphased VCFs if needed.
        specific_samples: List of specific sample IDs to extract.

    Returns:
        A tuple containing two lists: (paths_to_phased_vcfs, paths_to_unphased_vcfs).
        Paths will point to .vcf.gz files.
    """
    if extract:
        logger.warning("--extract is deprecated for VCF download. Sample extraction now uses bcftools and outputs .vcf.gz.")
        # We will ignore the extract flag for the main logic, as workflow expects .vcf.gz

    bcftools_available = check_bcftools_available()
    can_extract_samples = bcftools_available and specific_samples
    if specific_samples and not bcftools_available:
        logger.warning("bcftools not available - cannot extract specific samples. Will download multi-sample VCF files instead.")

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    downloaded_phased_vcfs = []
    downloaded_unphased_vcfs = []

    # --- Phased VCFs ---
    if phased_count > 0:
        logger.info(f"--- Processing Phased VCFs (Target: {phased_count}) ---")
        # Strategy: Download one main phased VCF, then extract samples if possible.
        # If not extracting samples, download multiple distinct phased VCFs.

        main_phased_vcf_remote = KNOWN_PHASED_VCFS[0]
        main_phased_filename = os.path.basename(main_phased_vcf_remote)
        main_phased_local_path = output_dir_path / f"phased_full_{main_phased_filename}"

        # Download the main multi-sample phased VCF if needed or if not extracting
        if not main_phased_local_path.exists() or not can_extract_samples:
             logger.info(f"Downloading main phased VCF: {main_phased_filename}")
             url = f"ftp://{FTP_SERVER}{main_phased_vcf_remote}"
             success = download_file_with_progress(url, str(main_phased_local_path), is_ftp=True)
             if not success:
                 logger.error(f"Failed to download main phased VCF: {main_phased_filename}. Skipping phased VCFs.")
                 # Skip phased processing if main download fails
                 phased_count = 0 # Set count to 0 to prevent further attempts

        if phased_count > 0 and main_phased_local_path.exists():
            if can_extract_samples:
                logger.info(f"Attempting to extract samples from {main_phased_local_path.name}")
                all_samples_in_vcf = get_samples_from_vcf(str(main_phased_local_path))
                if all_samples_in_vcf:
                    samples_to_extract = [s for s in specific_samples if s in all_samples_in_vcf]
                    logger.info(f"Found {len(samples_to_extract)} of the requested samples in the VCF: {', '.join(samples_to_extract)}")

                    if len(samples_to_extract) < len(specific_samples):
                        missing_req = [s for s in specific_samples if s not in samples_to_extract]
                        logger.warning(f"Requested samples not found in {main_phased_local_path.name}: {', '.join(missing_req)}")

                    # Limit extraction to phased_count
                    samples_to_extract = samples_to_extract[:phased_count]
                    logger.info(f"Will extract up to {phased_count} samples: {', '.join(samples_to_extract)}")

                    for i, sample_id in enumerate(samples_to_extract):
                        # Output filename includes sample ID and is compressed
                        sample_output_gz = output_dir_path / f"phased_{i+1}_{sample_id}_{main_phased_filename.replace('.vcf.gz', '.subset.vcf.gz')}"
                        if extract_sample_from_vcf(str(main_phased_local_path), str(sample_output_gz), sample_id):
                            downloaded_phased_vcfs.append(str(sample_output_gz))
                        else:
                            logger.warning(f"Failed to extract sample {sample_id} for phased VCF.")
                else:
                    logger.error(f"Could not get sample list from {main_phased_local_path.name}. Cannot extract samples.")
                    # Fallback: Add the main VCF if extraction failed but download succeeded
                    downloaded_phased_vcfs.append(str(main_phased_local_path))

            else: # Not extracting samples, download multiple files if needed
                logger.info("Adding multi-sample phased VCF(s).")
                downloaded_phased_vcfs.append(str(main_phased_local_path))
                # Download additional distinct phased VCFs if requested count > 1
                needed_more = phased_count - len(downloaded_phased_vcfs)
                for i in range(1, min(needed_more + 1, len(KNOWN_PHASED_VCFS))):
                    add_vcf_remote = KNOWN_PHASED_VCFS[i]
                    add_filename = os.path.basename(add_vcf_remote)
                    add_local_path = output_dir_path / f"phased_{i+1}_{add_filename}"
                    if not add_local_path.exists():
                         logger.info(f"Downloading additional phased VCF: {add_filename}")
                         url = f"ftp://{FTP_SERVER}{add_vcf_remote}"
                         if download_file_with_progress(url, str(add_local_path), is_ftp=True):
                             downloaded_phased_vcfs.append(str(add_local_path))
                         else:
                             logger.warning(f"Failed to download additional phased VCF: {add_filename}")
                    else:
                         logger.info(f"Additional phased VCF already exists: {add_local_path.name}")
                         downloaded_phased_vcfs.append(str(add_local_path))
                    # Ensure we don't exceed phased_count
                    if len(downloaded_phased_vcfs) >= phased_count:
                        break

    # --- Unphased VCFs ---
    if unphased_count > 0:
        logger.info(f"--- Processing Unphased VCFs (Target: {unphased_count}) ---")
        # Priority 1: Sample-specific VCFs if possible
        if can_extract_samples:
            logger.info("Attempting to download sample-specific unphased VCFs first.")
            sample_idx = 0
            for sample_id in specific_samples:
                if len(downloaded_unphased_vcfs) >= unphased_count:
                    break # Stop if we have enough

                if sample_id in SAMPLE_TO_UNPHASED_VCF:
                    # Try downloading one VCF per sample from their list
                    # Cycle through chromosomes 20, 21, 22 for variety if needed
                    vcf_options = SAMPLE_TO_UNPHASED_VCF[sample_id]
                    vcf_remote = vcf_options[sample_idx % len(vcf_options)]
                    filename = os.path.basename(vcf_remote)
                    local_path = output_dir_path / f"unphased_{sample_idx+1}_{sample_id}_{filename}"

                    if not local_path.exists():
                        logger.info(f"Downloading unphased VCF for {sample_id}: {filename}")
                        url = f"ftp://{FTP_SERVER}{vcf_remote}"
                        if download_file_with_progress(url, str(local_path), is_ftp=True):
                            downloaded_unphased_vcfs.append(str(local_path))
                            sample_idx += 1
                        else:
                            logger.warning(f"Failed to download unphased VCF for {sample_id}: {filename}")
                    else:
                        logger.info(f"Unphased VCF for {sample_id} already exists: {local_path.name}")
                        downloaded_unphased_vcfs.append(str(local_path))
                        sample_idx += 1
                else:
                    logger.warning(f"No known sample-specific unphased VCF path defined for {sample_id}")

        # Priority 2: Fallback to multi-sample unphased VCFs if needed
        unphased_needed = unphased_count - len(downloaded_unphased_vcfs)
        if unphased_needed > 0:
            logger.info(f"Need {unphased_needed} more unphased VCFs. Trying multi-sample fallback VCFs.")
            # Try downloading from fallback list
            for i in range(min(unphased_needed, len(KNOWN_UNPHASED_VCFS_FALLBACK))):
                vcf_remote = KNOWN_UNPHASED_VCFS_FALLBACK[i]
                filename = os.path.basename(vcf_remote)
                local_path = output_dir_path / f"unphased_fallback_{i+1}_{filename}"

                if not local_path.exists():
                    logger.info(f"Downloading fallback unphased VCF: {filename}")
                    url = f"ftp://{FTP_SERVER}{vcf_remote}"
                    if download_file_with_progress(url, str(local_path), is_ftp=True):
                        downloaded_unphased_vcfs.append(str(local_path))
                    else:
                        logger.warning(f"Failed to download fallback unphased VCF: {filename}")
                else:
                    logger.info(f"Fallback unphased VCF already exists: {local_path.name}")
                    downloaded_unphased_vcfs.append(str(local_path))

                if len(downloaded_unphased_vcfs) >= unphased_count:
                    break # Stop if we have enough

        # Final check if we still don't have enough (and force_unphased is not used much now)
        if len(downloaded_unphased_vcfs) < unphased_count:
             logger.warning(f"Could not obtain the target number ({unphased_count}) of unphased VCFs. Found {len(downloaded_unphased_vcfs)}.")


    return downloaded_phased_vcfs, downloaded_unphased_vcfs

def main():
    parser = argparse.ArgumentParser(description="Download VCF files (.vcf.gz) from 1000 Genomes. Optionally extracts samples using bcftools.")
    parser.add_argument('--output-dir', type=str, default='data/vcf',
                        help='Directory to save the downloaded VCF files (default: data/vcf)')
    parser.add_argument('--phased', type=int, default=1,
                        help='Max number of phased VCF files to obtain (per sample if specified, or multi-sample otherwise) (default: 1)')
    parser.add_argument('--unphased', type=int, default=1,
                        help='Max number of unphased VCF files to obtain (per sample if specified, or multi-sample otherwise) (default: 1)')
    parser.add_argument('--extract', action='store_true',
                        help='DEPRECATED: This flag is ignored. Sample extraction uses bcftools if --sample-ids is provided.')
    parser.add_argument('--force-unphased', action='store_true',
                        help='DEPRECATED: Fallback sources are now tried automatically if needed.')
    parser.add_argument('--sample-ids', type=str, nargs='+',
                        help='Specific sample IDs to extract using bcftools (requires bcftools). If not provided, multi-sample VCFs are downloaded.')
    args = parser.parse_args()

    if args.extract:
         logger.warning("--extract flag is deprecated and ignored.")
    if args.force_unphased:
         logger.warning("--force-unphased flag is deprecated and ignored.")


    phased_vcfs, unphased_vcfs = download_vcfs(
        args.output_dir,
        phased_count=args.phased,
        unphased_count=args.unphased,
        extract=False, # Ignored internally now
        force_unphased=False, # Ignored internally now
        specific_samples=args.sample_ids
    )

    logger.info("--- Download Summary ---")
    logger.info(f"Phased VCFs obtained ({len(phased_vcfs)}):")
    if phased_vcfs:
        for vcf in phased_vcfs:
            logger.info(f"  - {vcf}")
    else:
        logger.info("  (None)")

    logger.info(f"Unphased VCFs obtained ({len(unphased_vcfs)}):")
    if unphased_vcfs:
        for vcf in unphased_vcfs:
            logger.info(f"  - {vcf}")
    else:
        logger.info("  (None)")

    logger.info("------------------------")

if __name__ == "__main__":
    main()
