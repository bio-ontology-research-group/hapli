import argparse
import logging
from collections import defaultdict
from pathlib import Path
import re

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
    """Converts SAM/BAM alignments to a GFA graph with proper haplotype handling."""

    def __init__(self, sam_path: Path, ref_fasta_path: Path, output_gfa_path: Path, 
                 min_overlap: int = 100):
        self.sam_path = sam_path
        self.ref_fasta_path = ref_fasta_path
        self.output_gfa_path = output_gfa_path
        self.min_overlap = min_overlap
        
        self.segments = {}  # seq -> id
        self.next_segment_id = 1
        self.links = set()
        self.paths = defaultdict(list)
        self.haplotype_groups = defaultdict(list)  # sample -> [hp1_path, hp2_path]
        
    def _add_segment(self, seq: str) -> str:
        if not seq:
            return None
        if seq not in self.segments:
            seg_id = str(self.next_segment_id)
            self.segments[seq] = seg_id
            self.next_segment_id += 1
        return self.segments[seq]

    def _parse_haplotype_name(self, query_name: str):
        """
        Parse sample and haplotype from query name.
        Expected format: SAMPLE.hp1 or SAMPLE.hp2
        Returns: (sample_name, haplotype_number)
        """
        match = re.match(r'(.+)\.(hp[12])$', query_name)
        if match:
            return match.group(1), match.group(2)
        return query_name, None

    def _create_variant_graph(self, alignments_by_position, ref_seq, chrom):
        """
        Create a variant graph from overlapping alignments.
        This creates nodes at variation points and links between them.
        """
        # Sort alignments by start position
        sorted_alignments = sorted(alignments_by_position, key=lambda x: x.reference_start)
        
        if not sorted_alignments:
            return
            
        # Find common coordinate range
        min_start = min(aln.reference_start for aln in sorted_alignments)
        max_end = max(aln.reference_end for aln in sorted_alignments)
        
        # Extract reference sequence for this region
        ref_region = ref_seq[min_start:max_end] if len(ref_seq) > max_end else ref_seq[min_start:]
        
        # Create reference segment
        ref_seg_id = self._add_segment(ref_region)
        
        # Create reference path
        self.paths[f"{chrom}_reference"] = [(ref_seg_id, '+')]
        
        # Process each haplotype alignment
        for aln in sorted_alignments:
            sample_name, haplotype = self._parse_haplotype_name(aln.query_name)
            
            # Get the aligned sequence
            hap_seq = aln.query_sequence
            
            # Create segment for this haplotype sequence
            hap_seg_id = self._add_segment(hap_seq)
            
            # Create path for this haplotype
            path_name = f"{sample_name}.{haplotype}" if haplotype else aln.query_name
            self.paths[path_name] = [(hap_seg_id, '+')]
            
            # Track haplotype groupings
            if haplotype:
                self.haplotype_groups[sample_name].append(path_name)
            
            # Optionally add links between overlapping segments
            # This creates edges in the variation graph
            if ref_seg_id != hap_seg_id:
                # Add link from reference to variant and back
                # Using 0M overlap for now - adjust based on actual overlap
                self.links.add((ref_seg_id, '+', hap_seg_id, '+', '0M'))

    def convert(self):
        """
        Main conversion logic for haplotype-aware GFA creation.
        """
        logging.info("Starting SAM to GFA conversion with haplotype handling.")
        
        try:
            samfile = pysam.AlignmentFile(self.sam_path)
            ref_fasta = pysam.FastaFile(str(self.ref_fasta_path))
        except (ValueError, FileNotFoundError) as e:
            logging.error(f"Error opening input files: {e}")
            return

        # Get all mapped alignments
        mapped_alignments = [aln for aln in samfile.fetch(until_eof=True) if not aln.is_unmapped]
        
        if not mapped_alignments:
            logging.warning("No mapped alignments found in the input file.")
            self.write_gfa()
            return

        # Group alignments by chromosome and position for variant graph construction
        alignments_by_region = defaultdict(list)
        
        for aln in mapped_alignments:
            # Group by chromosome and approximate region (for overlapping alignments)
            region_key = (aln.reference_name, aln.reference_start // 10000)
            alignments_by_region[region_key].append(aln)
            
        # Process each region
        for (chrom, region_start), alignments in alignments_by_region.items():
            logging.info(f"Processing chromosome '{chrom}' region {region_start * 10000}")
            
            try:
                # Get reference sequence
                ref_seq = ref_fasta.fetch(chrom)
            except KeyError:
                logging.error(f"Chromosome '{chrom}' not found in reference FASTA.")
                continue
                
            # Create variant graph for this region
            self._create_variant_graph(alignments, ref_seq, chrom)

        # Log haplotype pairing information
        for sample, haplotypes in self.haplotype_groups.items():
            logging.info(f"Sample '{sample}' has haplotypes: {', '.join(sorted(haplotypes))}")

        self.write_gfa()
        logging.info(f"GFA file written to {self.output_gfa_path}")

    def write_gfa(self):
        """Write GFA with proper formatting and haplotype grouping metadata."""
        with open(self.output_gfa_path, 'w') as f:
            # Header
            f.write("H\tVN:Z:1.0\n")
            
            # Write sample grouping information as comments
            for sample, haplotypes in sorted(self.haplotype_groups.items()):
                f.write(f"# Sample: {sample} Haplotypes: {','.join(sorted(haplotypes))}\n")
            
            # Segments (S lines)
            sorted_segments = sorted(self.segments.items(), key=lambda item: int(item[1]))
            for seq, seg_id in sorted_segments:
                # Truncate very long sequences in output for readability
                seq_display = seq if len(seq) <= 1000 else f"{seq[:500]}...{seq[-500:]}"
                f.write(f"S\t{seg_id}\t{seq}\n")
            
            # Links (L lines) 
            sorted_links = sorted(list(self.links))
            for from_seg, from_orient, to_seg, to_orient, overlap in sorted_links:
                f.write(f"L\t{from_seg}\t{from_orient}\t{to_seg}\t{to_orient}\t{overlap}\n")

            # Paths (P lines) - grouped by sample
            # First write reference paths
            ref_paths = {k: v for k, v in self.paths.items() if 'reference' in k.lower()}
            for path_name, segments in sorted(ref_paths.items()):
                path_str = "".join([f"{seg_id}{orient}" for seg_id, orient in segments])
                cigars = "*"  # Can be computed if needed
                f.write(f"P\t{path_name}\t{path_str}\t{cigars}\n")
            
            # Then write haplotype paths grouped by sample
            processed = set()
            for sample in sorted(self.haplotype_groups.keys()):
                for hap_path in sorted(self.haplotype_groups[sample]):
                    if hap_path in self.paths:
                        segments = self.paths[hap_path]
                        path_str = "".join([f"{seg_id}{orient}" for seg_id, orient in segments])
                        cigars = "*"
                        f.write(f"P\t{hap_path}\t{path_str}\t{cigars}\n")
                        processed.add(hap_path)
            
            # Write any remaining paths not in haplotype groups
            for path_name, segments in sorted(self.paths.items()):
                if path_name not in processed and 'reference' not in path_name.lower():
                    path_str = "".join([f"{seg_id}{orient}" for seg_id, orient in segments])
                    cigars = "*"
                    f.write(f"P\t{path_name}\t{path_str}\t{cigars}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Convert SAM/BAM alignments to a GFA graph with haplotype handling."
    )
    parser.add_argument("-i", "--input-sam", required=True, type=Path,
                        help="Input SAM or BAM file.")
    parser.add_argument("-r", "--reference-fasta", required=True, type=Path,
                        help="Reference FASTA file (should be region-extracted).")
    parser.add_argument("-o", "--output-gfa", required=True, type=Path,
                        help="Output GFA file path.")
    parser.add_argument("--min-overlap", type=int, default=100,
                        help="Minimum overlap for creating links (default: 100).")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose logging.")
    args = parser.parse_args()

    setup_logging(args.verbose)
    
    converter = SamToGfaConverter(
        sam_path=args.input_sam,
        ref_fasta_path=args.reference_fasta,
        output_gfa_path=args.output_gfa,
        min_overlap=args.min_overlap
    )
    converter.convert()

if __name__ == "__main__":
    main()
    
