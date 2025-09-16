import argparse
import logging
import sys
from pathlib import Path

try:
    import gffutils
except ImportError:
    print("Error: gffutils is not installed. Please install it with 'pip install gffutils'")
    sys.exit(1)

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
    including their parent and child features.
    """
    logging.info(f"Filtering GFF file {gff_path} for region {chrom}:{start}-{end}")

    db_path = output_path.with_suffix(".db")
    if db_path.exists():
        db_path.unlink()

    logging.info("Creating temporary GFF database (this may take a moment)...")
    db = gffutils.create_db(str(gff_path), dbfn=str(db_path), force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True)

    features_to_keep = set()

    def add_relatives(feature, db, collection):
        """Recursively add a feature and all its parents and children to the collection."""
        collection.add(feature.id)
        for parent in db.parents(feature.id):
            if parent.id not in collection:
                add_relatives(parent, db, collection)
        for child in db.children(feature.id):
            if child.id not in collection:
                add_relatives(child, db, collection)

    logging.info("Identifying overlapping features and their relatives...")
    for feature in db.region(region=(chrom, start, end), completely_within=False):
        if feature.id not in features_to_keep:
            add_relatives(feature, db, features_to_keep)

    logging.info(f"Found {len(features_to_keep)} features (including relatives) to keep.")

    logging.info(f"Writing filtered GFF to {output_path}")
    with open(output_path, "w") as out_f:
        # Write original headers
        with open(gff_path) as in_f:
            for line in in_f:
                if line.startswith("##"):
                    out_f.write(line)
                else:
                    break
        
        # Write features from the database
        for feature_id in sorted(list(features_to_keep), key=lambda x: db[x].start):
            try:
                feature = db[feature_id]
                out_f.write(str(feature) + '\n')
            except gffutils.exceptions.FeatureNotFoundError:
                logging.warning(f"Could not find feature with ID '{feature_id}' in the database, skipping.")

    if db_path.exists():
        db_path.unlink()
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
