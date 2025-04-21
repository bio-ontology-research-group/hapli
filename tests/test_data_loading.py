import os
import pytest

# Define the paths relative to the project root
# Assumes tests are run from the project root directory
DATA_DIR = "data"
GFA_FILE = os.path.join(DATA_DIR, "example.gfa")
GFF_FILE = os.path.join(DATA_DIR, "example.gff3")
FASTA_FILE = os.path.join(DATA_DIR, "reference.fasta")

# List of files to check for convenience
TEST_FILES = [GFA_FILE, GFF_FILE, FASTA_FILE]

# Note: To see the print statements during test execution (even on pass),
# run pytest with the '-s' flag: pytest -s

def test_data_files_exist():
    """Check if the basic test data files exist and report their size."""
    print("\n--- Checking File Existence and Size ---")
    all_files_found = True
    for file_path in TEST_FILES:
        print(f"Checking: {file_path}")
        if os.path.exists(file_path):
            file_size = os.path.getsize(file_path)
            print(f"  -> Found. Size: {file_size} bytes.")
            assert True # Technically redundant, but confirms check
        else:
            print(f"  -> *** NOT FOUND ***")
            all_files_found = False
            # Fail assertion immediately if a file is missing
            pytest.fail(f"Required data file not found: {file_path}")

    if not all_files_found:
         pytest.fail("One or more data files were not found.")
    print("--- File Existence Check PASSED ---")


def test_read_and_stat_data_files():
    """Attempt to open and read each data file, reporting line counts."""
    print("\n--- Checking File Readability and Stats ---")
    stats = {}
    all_files_readable = True

    for file_path in TEST_FILES:
        print(f"Reading: {file_path}")
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
                line_count = len(lines)
                stats[file_path] = {'lines': line_count}
                print(f"  -> Read successful. Lines: {line_count}")
        except Exception as e:
            print(f"  -> *** FAILED TO READ ***: {e}")
            all_files_readable = False
            # Fail assertion immediately if a file cannot be read
            pytest.fail(f"Could not read data file {file_path}: {e}")

    if not all_files_readable:
        pytest.fail("One or more data files could not be read.")

    print("--- File Readability Check PASSED ---")
    print("Collected File Statistics:", stats)


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
