#!/usr/bin/env python3

import argparse
import os
import random
import sys
import subprocess # For optional validation call
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF # To write GFF3
# import gfapy # Removed gfapy import

# --- Default Settings ---
DEFAULT_REF_LENGTH = 1000
DEFAULT_NUM_FEATURES = 5
DEFAULT_PAIRED = False
DEFAULT_OUTPUT_DIR = "data"
DEFAULT_BASENAME = "example"
REF_SEQ_ID = "ref_chr1"
SEGMENT_LEN = 100 # Base length for GFA segments

# --- Helper Functions ---

def generate_random_sequence(length: int) -> str:
    """Generates a random DNA sequence."""
    return "".join(random.choice("ACGT") for _ in range(length))

def create_gff_feature(seq_id, start, end, strand, feature_type, feature_id, qualifiers=None):
    """Creates a Bio.SeqFeature object."""
    if qualifiers is None:
        qualifiers = {}
    qualifiers["ID"] = [feature_id]
    # Ensure 0-based start, exclusive end for FeatureLocation
    # GFF3 is 1-based, inclusive. Convert GFF3 coords before calling this.
    location = FeatureLocation(start, end, strand=strand)
    feature = SeqFeature(location, type=feature_type, id=feature_id, qualifiers=qualifiers)
    return feature

# --- Main Generation Function ---

def generate_data(ref_length: int, num_features: int, paired_haplotypes: bool, output_dir: str, basename: str):
    """Generates FASTA, GFF3, and GFA2 files."""

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    fasta_path = os.path.join(output_dir, f"{basename}.fasta")
    gff_path = os.path.join(output_dir, f"{basename}.gff3")
    gfa_path = os.path.join(output_dir, f"{basename}.gfa")

    print(f"Generating data with basename '{basename}' in '{output_dir}':")
    print(f"  Reference Length: {ref_length}")
    print(f"  Num Features: {num_features}")
    print(f"  Paired Haplotypes: {paired_haplotypes}")

    # 1. Generate FASTA
    print(f"  Generating FASTA: {fasta_path}")
    ref_seq_str = generate_random_sequence(ref_length)
    ref_record = SeqRecord(Seq(ref_seq_str), id=REF_SEQ_ID, description="Generated reference sequence")
    with open(fasta_path, "w") as f_out:
        f_out.write(ref_record.format("fasta"))

    # 2. Generate GFF3
    print(f"  Generating GFF3: {gff_path}")
    gff_records = []
    gff_features = []
    min_feature_gap = 10
    min_feature_len = 50
    max_feature_len = max(min_feature_len + 1, ref_length // (num_features + 1) - min_feature_gap)

    current_pos = 0
    for i in range(num_features):
        gene_id = f"gene{i+1}"
        mrna_id = f"mrna{i+1}"

        # Place gene start, ensuring space
        # Use 0-based coords internally, convert to 1-based for GFF output if needed by BCBio.GFF
        gene_start_0based = random.randint(current_pos, max(current_pos, ref_length - min_feature_len * 2 - min_feature_gap * 3))
        gene_len = random.randint(min_feature_len, max_feature_len)
        gene_end_0based = min(ref_length, gene_start_0based + gene_len)
        if gene_end_0based <= gene_start_0based: continue # Skip if too close to end

        strand = random.choice([1, -1])

        # Create Gene feature (using 0-based for FeatureLocation)
        gene_qual = {"source": ["generate_script"]}
        gene_feat = create_gff_feature(REF_SEQ_ID, gene_start_0based, gene_end_0based, strand, "gene", gene_id, gene_qual)

        # Create mRNA feature within gene
        mrna_qual = {"Parent": [gene_id], "source": ["generate_script"]}
        mrna_feat = create_gff_feature(REF_SEQ_ID, gene_start_0based, gene_end_0based, strand, "mRNA", mrna_id, mrna_qual)

        # Create Exon features within mRNA
        exon1_id = f"{mrna_id}_exon1"
        exon2_id = f"{mrna_id}_exon2"
        exon_len = max(10, (gene_end_0based - gene_start_0based) // 3)
        exon1_start_0based = gene_start_0based
        exon1_end_0based = min(gene_end_0based, exon1_start_0based + exon_len)
        exon2_start_0based = max(exon1_end_0based + 5, gene_end_0based - exon_len)
        exon2_end_0based = gene_end_0based

        if exon1_end_0based <= exon1_start_0based or exon2_end_0based <= exon2_start_0based or exon2_start_0based >= exon2_end_0based: continue # Skip if exons invalid

        exon_qual = {"Parent": [mrna_id], "source": ["generate_script"]}
        exon1_feat = create_gff_feature(REF_SEQ_ID, exon1_start_0based, exon1_end_0based, strand, "exon", exon1_id, exon_qual.copy())
        exon2_feat = create_gff_feature(REF_SEQ_ID, exon2_start_0based, exon2_end_0based, strand, "exon", exon2_id, exon_qual.copy())

        gff_features.extend([gene_feat, mrna_feat, exon1_feat, exon2_feat])

        current_pos = gene_end_0based + min_feature_gap # Update position for next feature

    # Create a SeqRecord for GFF writing
    gff_seq_record = SeqRecord(Seq(""), id=REF_SEQ_ID) # Sequence itself not needed here
    gff_seq_record.features = gff_features

    with open(gff_path, "w") as gff_out:
        # BCBio.GFF handles the conversion to 1-based GFF3 format during writing
        GFF.write([gff_seq_record], gff_out)


    # 3. Generate GFA2 Manually as Strings
    print(f"  Generating GFA2: {gfa_path}")
    gfa_lines = []
    gfa_lines.append("H\tVN:Z:2.0") # GFA2 Header

    # Create segments based on reference, introducing a variation point
    num_segments = (ref_length + SEGMENT_LEN - 1) // SEGMENT_LEN
    segments = [] # List of segment IDs
    segment_lengths = {} # Dict mapping seg_id -> length
    variation_point_idx = num_segments // 2 # Introduce variation around the middle
    alt_segment_id = None
    alt_segment_len = 0

    for i in range(num_segments):
        seg_id = f"s{i+1}"
        seg_start_0based = i * SEGMENT_LEN
        seg_end_0based = min(ref_length, (i + 1) * SEGMENT_LEN)
        seg_len = seg_end_0based - seg_start_0based
        if seg_len <= 0: continue

        seg_seq = ref_seq_str[seg_start_0based:seg_end_0based]
        segment_lengths[seg_id] = seg_len
        segments.append(seg_id)

        # Add reference segment line
        gfa_lines.append(f"S\t{seg_id}\t{seg_len}\t{seg_seq}")
        # Add Fragment line linking segment to reference
        # F <sid> <external> <sbeg> <send> <fbeg> <fend> <alignment>
        # Using 0-based, exclusive-end coordinates for external reference (sbeg, send)
        # Using 0-based, exclusive-end coordinates for segment (fbeg, fend)
        gfa_lines.append(f"F\t{seg_id}\t{REF_SEQ_ID}+\t{seg_start_0based}\t{seg_end_0based}\t0\t{seg_len}\t*")

        # Create alternative segment at variation point
        if i == variation_point_idx:
            alt_segment_id = f"s{i+1}_alt"
            alt_segment_len = seg_len # Use same length for simplicity
            alt_seq = generate_random_sequence(alt_segment_len)
            # Ensure alt_seq is different from original
            while alt_seq == seg_seq and alt_segment_len > 0:
                alt_seq = generate_random_sequence(alt_segment_len)

            # Add alternative segment line
            gfa_lines.append(f"S\t{alt_segment_id}\t{alt_segment_len}\t{alt_seq}")
            # Add Fragment line for alt segment (maps to same reference region)
            gfa_lines.append(f"F\t{alt_segment_id}\t{REF_SEQ_ID}+\t{seg_start_0based}\t{seg_end_0based}\t0\t{alt_segment_len}\t*")


    # Create Edges linking consecutive reference segments
    for i in range(len(segments) - 1):
        from_seg = segments[i]
        to_seg = segments[i+1]
        from_seg_len = segment_lengths[from_seg]
        # E <eid> <sid1> <sid2> <beg1> <end1> <beg2> <end2> <alignment>
        # End-to-start link: Use segment length for end position on sid1, 0 for start on sid2. Removed '$'.
        gfa_lines.append(f"E\t*\t{from_seg}+\t{to_seg}+\t{from_seg_len}\t{from_seg_len}\t0\t0\t*")

        # Add edges for the variation
        if i == variation_point_idx - 1 and alt_segment_id: # Link segment before variation point to alt segment
             gfa_lines.append(f"E\t*\t{from_seg}+\t{alt_segment_id}+\t{from_seg_len}\t{from_seg_len}\t0\t0\t*")
        if i == variation_point_idx and alt_segment_id: # Link alt segment to segment after variation point
             from_alt_seg = alt_segment_id
             # Need length of the alt segment
             gfa_lines.append(f"E\t*\t{from_alt_seg}+\t{to_seg}+\t{alt_segment_len}\t{alt_segment_len}\t0\t0\t*")


    # Create Ordered Groups (Paths)
    ref_path_elements = [f"{s}+" for s in segments]
    gfa_lines.append(f"O\t{REF_SEQ_ID}\t{' '.join(ref_path_elements)}") # Reference path

    # Sample paths
    sample1_path = list(ref_path_elements) # Sample 1 follows reference
    sample2_path = list(ref_path_elements) # Sample 2 uses variation
    if alt_segment_id and variation_point_idx < len(sample2_path):
        sample2_path[variation_point_idx] = f"{alt_segment_id}+"

    if paired_haplotypes:
        gfa_lines.append(f"O\tsample1_h1\t{' '.join(sample1_path)}")
        gfa_lines.append(f"O\tsample1_h2\t{' '.join(sample1_path)}") # Both follow ref
        gfa_lines.append(f"O\tsample2_h1\t{' '.join(sample1_path)}") # H1 follows ref
        gfa_lines.append(f"O\tsample2_h2\t{' '.join(sample2_path)}") # H2 uses variation
    else:
        gfa_lines.append(f"O\tsample1\t{' '.join(sample1_path)}")
        gfa_lines.append(f"O\tsample2\t{' '.join(sample2_path)}")


    # Write GFA file
    try:
        with open(gfa_path, "w") as f_out:
            for line in gfa_lines:
                f_out.write(line + "\n")
    except Exception as e:
         print(f"Error writing GFA file: {e}")
         sys.exit(1) # Exit if writing fails

    print("Data generation complete.")


# --- Main Execution Block ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate consistent FASTA, GFF3, and GFA2 example files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--ref-length", type=int, default=DEFAULT_REF_LENGTH,
        help="Approximate length of the reference sequence."
    )
    parser.add_argument(
        "--num-features", type=int, default=DEFAULT_NUM_FEATURES,
        help="Number of gene features to generate in GFF3."
    )
    parser.add_argument(
        "--paired-haplotypes", action='store_true', default=DEFAULT_PAIRED,
        help="Generate paired haplotypes (e.g., sample1_h1, sample1_h2) instead of single sample paths."
    )
    parser.add_argument(
        "--output-dir", type=str, default=DEFAULT_OUTPUT_DIR,
        help="Directory to save the generated files."
    )
    parser.add_argument(
        "--basename", type=str, default=DEFAULT_BASENAME,
        help="Basename for the output files (e.g., 'example' -> example.fasta, example.gff3, example.gfa)."
    )

    args = parser.parse_args()

    # Basic validation
    if args.ref_length < 100:
        print("Error: Reference length must be at least 100.", file=sys.stderr)
        sys.exit(1)
    if args.num_features < 1:
        print("Error: Number of features must be at least 1.", file=sys.stderr)
        sys.exit(1)
    if SEGMENT_LEN * args.num_features > args.ref_length:
         print(f"Warning: Reference length ({args.ref_length}) might be too small for {args.num_features} features and segment length {SEGMENT_LEN}.", file=sys.stderr)


    generate_data(
        ref_length=args.ref_length,
        num_features=args.num_features,
        paired_haplotypes=args.paired_haplotypes,
        output_dir=args.output_dir,
        basename=args.basename
    )

    # Optional: Validate generated files immediately after generation
    try:
        gfa_output_path = os.path.join(args.output_dir, f"{args.basename}.gfa")
        validator_script = os.path.join(os.path.dirname(__file__), "validate_gfa2.py")
        if os.path.exists(validator_script) and os.path.exists(gfa_output_path):
            print("\nRunning GFA2 validation script...")
            result = subprocess.run(
                [sys.executable, validator_script, gfa_output_path],
                capture_output=True, text=True, check=False, encoding='utf-8'
            )
            if result.returncode == 0:
                print("GFA2 Validation successful.")
            else:
                print("GFA2 Validation failed:")
                print("--- Stderr ---")
                print(result.stderr)
                print("--- Stdout ---")
                print(result.stdout)
                # Decide if failure here should stop the script
                # sys.exit(1)
        else:
            print("\nSkipping GFA2 validation (script or GFA file not found).")
    except Exception as e:
        print(f"\nCould not run validation script: {e}")
