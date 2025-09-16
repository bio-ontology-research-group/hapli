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
        """Main conversion logic."""
        logging.info("Starting SAM to GFA conversion.")
        
        samfile = pysam.AlignmentFile(self.sam_path)
        ref_fasta = pysam.FastaFile(str(self.ref_fasta_path))

        alignments = list(samfile.fetch())
        if not alignments:
            logging.warning("No alignments found.")
            return

        chrom = alignments[0].reference_name
        ref_start = min(aln.reference_start for aln in alignments if not aln.is_unmapped)
        ref_end = max(aln.reference_end for aln in alignments if not aln.is_unmapped)
        
        logging.info(f"Processing region {chrom}:{ref_start}-{ref_end}")

        breakpoints = self._collect_breakpoints(alignments, ref_fasta, chrom, ref_start, ref_end)
        sorted_bps = sorted(list(breakpoints))

        self._build_haplotype_paths(alignments, ref_fasta, chrom, sorted_bps)
        self._build_reference_path(ref_fasta, chrom, sorted_bps)

        self.write_gfa()
        logging.info(f"GFA file written to {self.output_gfa_path}")

    def _collect_breakpoints(self, alignments, ref_fasta, chrom, ref_start, ref_end):
        breakpoints = {ref_start, ref_end}
        for aln in alignments:
            if aln.is_unmapped:
                continue
            
            # Indel breakpoints from CIGAR
            ref_pos = aln.reference_start
            for op, length in aln.cigartuples:
                if op in (0, 7, 8):  # M, X, =
                    ref_pos += length
                elif op == 1:  # Insertion
                    breakpoints.add(ref_pos)
                elif op == 2:  # Deletion
                    breakpoints.add(ref_pos)
                    breakpoints.add(ref_pos + length)
                    ref_pos += length
                elif op == 3: # N
                    ref_pos += length

            # SNV breakpoints
            for q_idx, r_idx in aln.get_aligned_pairs():
                if q_idx is not None and r_idx is not None:
                    if ref_start <= r_idx < ref_end:
                        q_base = aln.query_sequence[q_idx]
                        r_base = ref_fasta.fetch(chrom, r_idx, r_idx + 1).upper()
                        if q_base != r_base:
                            breakpoints.add(r_idx)
                            breakpoints.add(r_idx + 1)
        return breakpoints

    def _build_haplotype_paths(self, alignments, ref_fasta, chrom, sorted_bps):
        for aln in alignments:
            if aln.is_unmapped:
                continue

            path_name = aln.query_name
            current_path_segs = []
            
            ref_pos = aln.reference_start
            q_pos = 0

            for op, length in aln.cigartuples:
                if op in (0, 7, 8):  # Match/mismatch
                    op_ref_end = ref_pos + length
                    
                    bps_in_match = [bp for bp in sorted_bps if ref_pos < bp < op_ref_end]
                    block_bps = [ref_pos] + bps_in_match + [op_ref_end]

                    for i in range(len(block_bps) - 1):
                        start, end = block_bps[i], block_bps[i+1]
                        q_start_offset = start - ref_pos
                        q_end_offset = end - ref_pos
                        q_sub_seq = aln.query_sequence[q_pos + q_start_offset : q_pos + q_end_offset]
                        
                        seg_id = self._add_segment(q_sub_seq)
                        if seg_id:
                            current_path_segs.append((seg_id, '+'))
                    
                    ref_pos = op_ref_end
                    q_pos += length
                
                elif op == 1:  # Insertion
                    ins_seq = aln.query_sequence[q_pos:q_pos + length]
                    seg_id = self._add_segment(ins_seq)
                    if seg_id:
                        current_path_segs.append((seg_id, '+'))
                    q_pos += length

                elif op == 2:  # Deletion
                    ref_pos += length

                elif op == 4:  # Soft clip
                    q_pos += length
            
            self.paths[path_name] = current_path_segs
            for i in range(len(current_path_segs) - 1):
                from_seg, from_orient = current_path_segs[i]
                to_seg, to_orient = current_path_segs[i+1]
                self.links.add((from_seg, from_orient, to_seg, to_orient, '0M'))

    def _build_reference_path(self, ref_fasta, chrom, sorted_bps):
        ref_path_segs = []
        for i in range(len(sorted_bps) - 1):
            start, end = sorted_bps[i], sorted_bps[i+1]
            if start < end:
                seq = ref_fasta.fetch(chrom, start, end)
                seg_id = self._add_segment(seq)
                if seg_id:
                    ref_path_segs.append((seg_id, '+'))
        
        self.paths['reference'] = ref_path_segs
        for i in range(len(ref_path_segs) - 1):
            from_seg, from_orient = ref_path_segs[i]
            to_seg, to_orient = ref_path_segs[i+1]
            self.links.add((from_seg, from_orient, to_seg, to_orient, '0M'))

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
