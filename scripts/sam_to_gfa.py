import argparse
import logging
from collections import defaultdict
from pathlib import Path

try:
    import pysam
except ImportError:
    print("Error: pysam is not installed. Please install it with 'pip install pysam'")
    exit(1)

def setup_logging(verbose: bool = False):
    """Set up logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")

class SamToGfaConverter:
    """Converts SAM/BAM alignments to a GFA graph."""

    def __init__(self, sam_path: Path, ref_fasta_path: Path, output_gfa_path: Path):
        self.sam_path = sam_path
        self.ref_fasta_path = ref_fasta_path
        self.output_gfa_path = output_gfa_path
        
        self.segments = {}  # seq -> id
        self.next_segment_id = 1
        self.links = set()
        self.paths = defaultdict(list)

    def _add_segment(self, seq: str) -> str:
        if not seq:
            return None
        if seq not in self.segments:
            seg_id = str(self.next_segment_id)
            self.segments[seq] = seg_id
            self.next_segment_id += 1
        return self.segments[seq]

    def convert(self):
        """
        Main conversion logic.
        Creates a simple GFA with a single path representing the reference
        sequence over the region covered by the SAM/BAM file.
        """
        logging.info("Starting SAM to GFA conversion.")
        
        try:
            samfile = pysam.AlignmentFile(self.sam_path)
            ref_fasta = pysam.FastaFile(str(self.ref_fasta_path))
        except (ValueError, FileNotFoundError) as e:
            logging.error(f"Error opening input files: {e}")
            return

        mapped_alignments = [aln for aln in samfile.fetch(until_eof=True) if not aln.is_unmapped]
        
        if not mapped_alignments:
            logging.warning("No mapped alignments found in the input file.")
            self.write_gfa()
            return

        chrom = mapped_alignments[0].reference_name
        ref_start = min(aln.reference_start for aln in mapped_alignments)
        ref_end = max(aln.reference_end for aln in mapped_alignments)
        
        logging.info(f"Determined region from SAM file: {chrom}:{ref_start}-{ref_end}")

        try:
            ref_seq = ref_fasta.fetch(chrom, ref_start, ref_end)
        except KeyError:
            logging.error(f"Chromosome '{chrom}' not found in reference FASTA file.")
            return
        
        if not ref_seq:
            logging.warning("Fetched empty reference sequence for the region. GFA will be empty.")
            self.write_gfa()
            return

        seg_id = self._add_segment(ref_seq)
        
        if seg_id:
            # Create a path named after the chromosome.
            self.paths[chrom] = [(seg_id, '+')]
            # Also create a generic 'reference' path for compatibility.
            self.paths['reference'] = [(seg_id, '+')]

        self.write_gfa()
        logging.info(f"GFA file written to {self.output_gfa_path}")

    def write_gfa(self):
        with open(self.output_gfa_path, 'w') as f:
            f.write("H\tVN:Z:1.0\n")
            
            sorted_segments = sorted(self.segments.items(), key=lambda item: int(item[1]))
            for seq, seg_id in sorted_segments:
                f.write(f"S\t{seg_id}\t{seq}\n")
            
            sorted_links = sorted(list(self.links))
            for from_seg, from_orient, to_seg, to_orient, overlap in sorted_links:
                f.write(f"L\t{from_seg}\t{from_orient}\t{to_seg}\t{to_orient}\t{overlap}\n")

            sorted_paths = sorted(self.paths.items())
            for path_name, segments in sorted_paths:
                path_str = "".join([f"{seg_id}{orient}" for seg_id, orient in segments])
                f.write(f"P\t{path_name}\t{path_str}\t*\n")

def main():
    parser = argparse.ArgumentParser(description="Convert SAM/BAM alignments to a GFA graph.")
    parser.add_argument("-i", "--input-sam", required=True, type=Path,
                        help="Input SAM or BAM file.")
    parser.add_argument("-r", "--reference-fasta", required=True, type=Path,
                        help="Reference FASTA file.")
    parser.add_argument("-o", "--output-gfa", required=True, type=Path,
                        help="Output GFA file path.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose logging.")
    args = parser.parse_args()

    setup_logging(args.verbose)
    
    converter = SamToGfaConverter(
        sam_path=args.input_sam,
        ref_fasta_path=args.reference_fasta,
        output_gfa_path=args.output_gfa
    )
    converter.convert()

if __name__ == "__main__":
    main()
