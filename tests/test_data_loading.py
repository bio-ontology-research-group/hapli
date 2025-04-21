import os
import pytest

# Define the paths relative to the project root
DATA_DIR = "data"
GFA_FILE = os.path.join(DATA_DIR, "example.gfa")
GFF_FILE = os.path.join(DATA_DIR, "example.gff3")
FASTA_FILE = os.path.join(DATA_DIR, "reference.fasta")

def test_data_files_exist():
    """Check if the basic test data files exist."""
    assert os.path.exists(GFA_FILE), f"GFA file not found: {GFA_FILE}"
    assert os.path.exists(GFF_FILE), f"GFF3 file not found: {GFF_FILE}"
    assert os.path.exists(FASTA_FILE), f"FASTA file not found: {FASTA_FILE}"

def test_can_read_data_files():
    """Attempt to open and read a small part of each data file."""
    try:
        with open(GFA_FILE, 'r') as f:
            f.readline()
    except Exception as e:
        pytest.fail(f"Could not read GFA file {GFA_FILE}: {e}")

    try:
        with open(GFF_FILE, 'r') as f:
            f.readline()
    except Exception as e:
        pytest.fail(f"Could not read GFF3 file {GFF_FILE}: {e}")

    try:
        with open(FASTA_FILE, 'r') as f:
            f.readline()
    except Exception as e:
        pytest.fail(f"Could not read FASTA file {FASTA_FILE}: {e}")

# Future tests could involve actual parsing using relevant libraries
# def test_parse_gfa():
#     # Add code to parse GFA_FILE using a GFA parsing library
#     pass
#
# def test_parse_gff():
#     # Add code to parse GFF_FILE using a GFF parsing library (e.g., Biopython)
#     pass
#
# def test_parse_fasta():
#     # Add code to parse FASTA_FILE using a FASTA parsing library (e.g., Biopython)
#     pass

