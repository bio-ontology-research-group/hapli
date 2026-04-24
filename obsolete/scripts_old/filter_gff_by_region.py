import argparse
import logging
import sys
from pathlib import Path

try:
    import pysam
except ImportError:
    print("Error: pysam is not installed. Please install it with 'pip install pysam'")
    sys.exit(1)


def setup_logging(verbose: bool = False):
    """Set up logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")


def get_sam_region(sam_path: Path) -> tuple[str, int, int]:
    """Extracts the genomic region covered by alignments in a SAM/BAM file."""
    logging.info(f"Reading SAM/BAM file to determine region: {sam_path}")
    samfile = pysam.AlignmentFile(sam_path, "r")

    chrom = None
    min_start = float('inf')
    max_end = float('-inf')
    count = 0

    for alignment in samfile.fetch(until_eof=True):
        if alignment.is_unmapped:
            continue

        if chrom is None:
            chrom = alignment.reference_name
        elif chrom != alignment.reference_name:
            # This script assumes a single contiguous region on one chromosome.
            logging.warning(
                f"Multiple chromosomes found. Using first one: {chrom}. "
                f"Ignoring alignment on {alignment.reference_name}"
            )
            continue

        min_start = min(min_start, alignment.reference_start)
        max_end = max(max_end, alignment.reference_end)
        count += 1

    if count == 0:
        raise ValueError("No aligned reads found in the SAM/BAM file.")

    logging.info(f"Determined region from {count} alignments: {chrom}:{min_start}-{max_end}")
    return chrom, int(min_start), int(max_end)


def filter_gff(gff_path: Path, output_path: Path, chrom: str, start: int, end: int):
    """
    Filters a GFF3 file to include only features that overlap the given region.
    This is a single-pass filter that streams the file and does not load it
    into memory or build a database. It only keeps features that physically
    lie within the specified coordinates.
    """
    logging.info(f"Filtering GFF file {gff_path} for region {chrom}:{start}-{end}")
    
    count = 0
    with open(gff_path, 'r') as in_f, open(output_path, 'w') as out_f:
        for line in in_f:
            if line.startswith('#'):
                out_f.write(line)
                continue

            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue

            line_chrom, _, _, line_start_str, line_end_str, _, _, _, _ = parts
            
            if line_chrom == chrom:
                try:
                    line_start = int(line_start_str)
                    line_end = int(line_end_str)
                except ValueError:
                    logging.warning(f"Could not parse coordinates in line: {line.strip()}")
                    continue
                
                # Check for overlap: region_A_end >= region_B_start AND region_A_start <= region_B_end
                if line_end >= start and line_start <= end:
                    out_f.write(line)
                    count += 1

    logging.info(f"Found and wrote {count} features within the region.")
    logging.info(f"Filtered GFF written to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter a GFF3 file to a specific genomic region defined by a SAM/BAM file. "
                    "Keeps only features that physically overlap the region."
    )
    parser.add_argument("--sam-file", required=True, type=Path, help="Input SAM or BAM file to define the region.")
    parser.add_argument("--gff-file", required=True, type=Path, help="Input GFF3 file to filter.")
    parser.add_argument("--output-file", required=True, type=Path, help="Path for the output filtered GFF3 file.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging.")
    args = parser.parse_args()

    setup_logging(args.verbose)

    try:
        chrom, start, end = get_sam_region(args.sam_file)
        filter_gff(args.gff_file, args.output_file, chrom, start, end)
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
