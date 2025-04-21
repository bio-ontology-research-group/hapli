#!/usr/bin/env python3

import pygfa
import os

# Define the output path relative to this script's location
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_GFA_PATH = os.path.join(BASE_DIR, "data", "example.gfa")

print(f"Generating GFA2 file at: {OUTPUT_GFA_PATH}")

# Create a GFA graph object
gfa = pygfa.graph.GraphicalFragmentAssembly()

# Add Header (pygfa adds VN:Z:2.0 automatically for GFA2 context)
# If specific header tags are needed, add them like:
# gfa.add_header("VN", "Z", "2.0") # pygfa might handle this automatically
# gfa.add_header("TS", "i", "5000") # Example optional tag

# Add Segments
gfa.add_segment("1", sequence="A") # Length is inferred from sequence
gfa.add_segment("2", sequence="T")
gfa.add_segment("3", sequence="C")
gfa.add_segment("4", sequence="G")
gfa.add_segment("5", sequence="X") # Alternative segment

# Add Edges (Links in GFA1)
# pygfa uses Edge(edge_id, from_segment, from_orient, to_segment, to_orient,
#                 begin1, end1, begin2, end2, alignment)
# Using '*' for edge_id
# Assuming simple end-to-start links (overlap '*')
# pygfa's Edge coordinates seem to be 1-based end position for first segment,
# and 0-based start position for second segment, based on how it writes GFA2.
# Let's use the segment length for end1 and 0 for beg2.
# The alignment '*' represents the simple overlap.
gfa.add_edge(pygfa.types.Edge("*", "1", "+", "2", "+", gfa.segments["1"].slen, gfa.segments["1"].slen, 0, 0, "*"))
gfa.add_edge(pygfa.types.Edge("*", "2", "+", "3", "+", gfa.segments["2"].slen, gfa.segments["2"].slen, 0, 0, "*"))
gfa.add_edge(pygfa.types.Edge("*", "3", "+", "4", "+", gfa.segments["3"].slen, gfa.segments["3"].slen, 0, 0, "*"))
gfa.add_edge(pygfa.types.Edge("*", "2", "+", "5", "+", gfa.segments["2"].slen, gfa.segments["2"].slen, 0, 0, "*")) # Variation link

# Add Ordered Groups (Paths in GFA1)
# pygfa uses OrderedGroup(group_id, elements) where elements are segment references
# Elements need to be specified with orientation, e.g., ('1', '+')
gfa.add_ordered_group("sample1", [("1", "+"), ("2", "+"), ("3", "+"), ("4", "+")])
gfa.add_ordered_group("sample2", [("1", "+"), ("2", "+"), ("5", "+"), ("4", "+")])

# Write the GFA file
try:
    # Ensure the data directory exists
    os.makedirs(os.path.dirname(OUTPUT_GFA_PATH), exist_ok=True)

    with open(OUTPUT_GFA_PATH, "w") as f:
        # Use GFA2 context for writing
        f.write(str(gfa)) # pygfa's __str__ should produce GFA2 format
    print("Successfully generated GFA2 file.")

except Exception as e:
    print(f"Error writing GFA file: {e}")
    exit(1)

# Optional: Validate the generated file immediately
try:
    import subprocess
    import sys
    validator_script = os.path.join(BASE_DIR, "scripts", "validate_gfa2.py")
    if os.path.exists(validator_script):
        print("Running validation script...")
        result = subprocess.run(
            [sys.executable, validator_script, OUTPUT_GFA_PATH],
            capture_output=True, text=True, check=False, encoding='utf-8'
        )
        if result.returncode == 0:
            print("Validation successful.")
        else:
            print("Validation failed:")
            print("--- Stderr ---")
            print(result.stderr)
            print("--- Stdout ---")
            print(result.stdout)
            exit(1) # Exit if validation fails after generation
    else:
        print("Validation script not found, skipping validation.")
except Exception as e:
    print(f"Could not run validation script: {e}")

