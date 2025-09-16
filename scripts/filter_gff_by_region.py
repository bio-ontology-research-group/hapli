import argparse
import logging
import sys
from collections import defaultdict
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
    Filters a GFF3 file to include only features overlapping the given region,
    including their parent and child features. This version avoids building a
    gffutils database for performance by performing two passes over the file.
    """
    logging.info(f"Filtering GFF file {gff_path} for region {chrom}:{start}-{end}")
    logging.info("Pass 1: Identifying overlapping features and building relationship map...")

    child_to_parents = defaultdict(list)
    parent_to_children = defaultdict(list)
    overlapping_feature_ids = set()

    def get_attrs(attr_str):
        attrs = {}
        for part in attr_str.split(';'):
            if '=' in part:
                key, val = part.split('=', 1)
                attrs[key] = val
        return attrs

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue

            line_chrom, _, _, line_start_str, line_end_str, _, _, _, attrs_str = parts
            
            attrs = get_attrs(attrs_str)
            feature_id = attrs.get('ID')
            if not feature_id:
                continue

            parent_id_val = attrs.get('Parent')
            if parent_id_val:
                parent_ids = parent_id_val.split(',')
                child_to_parents[feature_id].extend(parent_ids)
                for p_id in parent_ids:
                    parent_to_children[p_id].append(feature_id)

            line_start, line_end = int(line_start_str), int(line_end_str)
            if line_chrom == chrom and line_end >= start and line_start <= end:
                overlapping_feature_ids.add(feature_id)

    logging.info(f"Found {len(overlapping_feature_ids)} directly overlapping features.")
    logging.info("Expanding set to include all relatives...")

    features_to_keep = set()
    queue = list(overlapping_feature_ids)
    
    while queue:
        feature_id = queue.pop(0)
        if feature_id in features_to_keep:
            continue
        
        features_to_keep.add(feature_id)

        for p_id in child_to_parents.get(feature_id, []):
            if p_id not in features_to_keep:
                queue.append(p_id)
        
        for c_id in parent_to_children.get(feature_id, []):
            if c_id not in features_to_keep:
                queue.append(c_id)

    logging.info(f"Total features to keep (including relatives): {len(features_to_keep)}")

    logging.info(f"Pass 2: Writing filtered GFF to {output_path}")
    with open(output_path, "w") as out_f, open(gff_path) as in_f:
        for line in in_f:
            if line.startswith("#"):
                out_f.write(line)
                continue
            
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue
            
            attrs = get_attrs(parts[8])
            feature_id = attrs.get('ID')
            if feature_id and feature_id in features_to_keep:
                out_f.write(line)

    logging.info("Filtering complete.")


def main():
    parser = argparse.ArgumentParser(
        description="Filter a GFF3 file to a specific genomic region defined by a SAM/BAM file. "
                    "Keeps all features that overlap the region, plus their parent and child features."
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
