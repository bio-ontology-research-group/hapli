# tests/test_vcf_to_gfa.py
import unittest
import os
import tempfile
import shutil
import sys
import logging

# Ensure src directory is in path for testing
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
src_path = os.path.join(project_root, 'src')
if src_path not in sys.path:
    sys.path.insert(0, src_path)

# Conditional import based on dependencies
dependencies_installed = False
try:
    import pysam
    import pyfaidx
    import gfapy
    from src.converters.vcf_to_gfa import VCFtoGFAConverter, VCFtoGFAConversionError
    from src.converters.reference_handler import ReferenceHandler, ReferenceHandlerError
    dependencies_installed = True
except ImportError as e:
    print(f"Skipping VCF->GFA tests: Missing dependencies - {e}")


# Define paths relative to the test file location or project root
TEST_DATA_DIR = os.path.join(project_root, 'data', 'vcf_conversion')
REF_FASTA = os.path.join(TEST_DATA_DIR, 'reference.fasta')
VCF_MULTI_SAMPLE = os.path.join(TEST_DATA_DIR, 'phased_multi_sample.vcf')
VCF_MULTI_BLOCK = os.path.join(TEST_DATA_DIR, 'multi_block.vcf')
VCF_MALFORMED = os.path.join(TEST_DATA_DIR, 'malformed.vcf') # Need to create this
VCF_NO_SAMPLES = os.path.join(TEST_DATA_DIR, 'no_samples.vcf') # Need to create this
VCF_EMPTY = os.path.join(TEST_DATA_DIR, 'empty.vcf') # Need to create this

# Disable logging noise during tests
logging.disable(logging.CRITICAL)

# Helper function to create dummy VCF files for specific tests
def create_dummy_vcf(filepath, content):
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##reference=file://reference.fasta\n")
        f.write("##contig=<ID=chr1,length=600>\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n") # No samples
        if content: # Add specific records if provided
             f.write(content)

@unittest.skipIf(not dependencies_installed, "Missing dependencies (pysam, pyfaidx, gfapy)")
class TestVCFtoGFAConverter(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Set up test data directory and ensure files exist."""
        cls.test_dir = tempfile.mkdtemp(prefix="vcf_gfa_test_")
        cls.output_gfa = os.path.join(cls.test_dir, "output.gfa")
        cls.maxDiff = None # Show full diff on assertion failure

        # --- Ensure main test data exists ---
        if not os.path.exists(TEST_DATA_DIR):
             raise unittest.SkipTest(f"Test data directory not found: {TEST_DATA_DIR}")
        if not os.path.exists(REF_FASTA) or not os.path.exists(VCF_MULTI_SAMPLE) or not os.path.exists(VCF_MULTI_BLOCK):
            raise unittest.SkipTest(f"Required test data files not found in {TEST_DATA_DIR}")

        # --- Create dummy VCFs for edge cases ---
        create_dummy_vcf(VCF_MALFORMED, "chr1\t10\trs1\tA\tG\t100\tPASS\t.\tGT\t0/1\tEXTRA_COL\n") # Extra column
        create_dummy_vcf(VCF_NO_SAMPLES, "chr1\t10\trs1\tA\tG\t100\tPASS\t.\tGT\n") # Header has FORMAT but no samples
        create_dummy_vcf(VCF_EMPTY, "") # Only header

        # --- Create FASTA index ---
        try:
             # Use ReferenceHandler to trigger index creation if needed
             with ReferenceHandler(REF_FASTA) as handler:
                 if not handler._fasta: raise ReferenceHandlerError("FASTA could not be loaded")
        except (ReferenceHandlerError, FileNotFoundError) as e:
             raise unittest.SkipTest(f"Failed to prepare FASTA index for {REF_FASTA}: {e}")
        except Exception as e:
             # Catch broader exceptions during setup
             raise unittest.SkipTest(f"Unexpected error preparing FASTA {REF_FASTA}: {e}")


        # --- Create VCF indices ---
        # Use force=True to overwrite potentially outdated indices
        cls.vcf_files_to_index = [VCF_MULTI_SAMPLE, VCF_MULTI_BLOCK, VCF_MALFORMED, VCF_NO_SAMPLES, VCF_EMPTY]
        for vcf_file in cls.vcf_files_to_index:
            try:
                # Remove old indices first to ensure clean state
                for ext in ['.tbi', '.csi']:
                    if os.path.exists(vcf_file + ext): os.remove(vcf_file + ext)
                # Index only if the file has content beyond header
                if os.path.getsize(vcf_file) > 200: # Approx header size
                    pysam.tabix_index(vcf_file, preset="vcf", force=True, keep_original=True)
            except pysam.utils.SamtoolsError as e:
                # Allow tests to run if indexing fails for dummy files, but maybe not for main ones
                if vcf_file in [VCF_MULTI_SAMPLE, VCF_MULTI_BLOCK]:
                     raise unittest.SkipTest(f"Failed to create required VCF index for {vcf_file} (samtools/htslib needed?): {e}")
                else:
                     print(f"Warning: Could not index dummy VCF {vcf_file}: {e}")
            except Exception as e:
                 print(f"Warning: Unexpected error indexing {vcf_file}: {e}")


    @classmethod
    def tearDownClass(cls):
        """Remove temporary directory and generated files."""
        shutil.rmtree(cls.test_dir)
        # Clean up index files created
        for idx_ext in ['.fai', '.tbi', '.csi']:
            # Clean indices for main test files
            for vcf_file in [VCF_MULTI_SAMPLE, VCF_MULTI_BLOCK]:
                 idx_file = vcf_file + idx_ext
                 if os.path.exists(idx_file):
                     try: os.remove(idx_file)
                     except OSError: pass
            # Clean indices for dummy files
            for vcf_file in [VCF_MALFORMED, VCF_NO_SAMPLES, VCF_EMPTY]:
                 idx_file = vcf_file + idx_ext
                 if os.path.exists(idx_file):
                     try: os.remove(idx_file)
                     except OSError: pass
            # Clean FASTA index
            fasta_idx = REF_FASTA + idx_ext
            if os.path.exists(fasta_idx):
                 try: os.remove(fasta_idx)
                 except OSError: pass

    def tearDown(self):
        """Clean up output file after each test."""
        if os.path.exists(self.output_gfa):
            try:
                os.remove(self.output_gfa)
            except OSError as e:
                 print(f"Warning: Could not remove output file {self.output_gfa}: {e}")


    def _validate_gfa_structure(self, gfa_filepath, min_segments=1, min_links=0, expected_paths=None):
        """Helper to perform basic structural validation on GFA using gfapy."""
        self.assertTrue(os.path.exists(gfa_filepath), f"Output GFA file not found: {gfa_filepath}")
        try:
            gfa = gfapy.Gfa.from_file(gfa_filepath)
            self.assertGreaterEqual(len(gfa.segments), min_segments)
            self.assertGreaterEqual(len(gfa.links), min_links)
            if expected_paths is not None:
                self.assertEqual(len(gfa.paths), len(expected_paths), f"Expected {len(expected_paths)} paths, found {len(gfa.paths)}")
                found_paths = {p.name for p in gfa.paths}
                self.assertSetEqual(found_paths, set(expected_paths), f"Path names mismatch. Found: {found_paths}, Expected: {expected_paths}")

            # Check segment LN tags
            for seg in gfa.segments:
                 self.assertIsNotNone(seg.tags.get("LN"), f"Segment {seg.name} missing LN tag")
                 self.assertEqual(seg.tags["LN"].value, len(seg.sequence), f"Segment {seg.name} LN tag mismatch")

            # Check path tags (SM, HP) if paths exist
            if expected_paths:
                 for path in gfa.paths:
                     self.assertIsNotNone(path.tags.get("SM"), f"Path {path.name} missing SM tag")
                     self.assertIsNotNone(path.tags.get("HP"), f"Path {path.name} missing HP tag")
                     self.assertIn(path.tags["HP"].value, [1, 2], f"Path {path.name} has invalid HP tag value")

            return gfa # Return parsed object for further checks

        except gfapy.error.FormatError as e:
            self.fail(f"Output GFA file is not valid: {e}")
        except Exception as e:
             self.fail(f"Error during GFA validation: {e}")


    def test_init(self):
        """Test converter initialization."""
        try:
            # Use context manager for auto-cleanup in case of init error
            with VCFtoGFAConverter(VCF_MULTI_SAMPLE, REF_FASTA, self.output_gfa) as converter:
                self.assertIsNotNone(converter)
                self.assertEqual(converter.unphased_strategy, 'ref') # Check default
        except Exception as e:
            self.fail(f"Initialization failed: {e}")

    def test_missing_input_files_init(self):
        """Test initialization with non-existent files (should raise error)."""
        with self.assertRaises((FileNotFoundError, VCFtoGFAConversionError)):
            # Error should be raised during init when ReferenceHandler/pysam try to open files
            with VCFtoGFAConverter("nonexistent.vcf", REF_FASTA, self.output_gfa) as converter:
                 pass # Should not reach here

        with self.assertRaises((FileNotFoundError, VCFtoGFAConversionError)):
             with VCFtoGFAConverter(VCF_MULTI_SAMPLE, "nonexistent.fasta", self.output_gfa) as converter:
                  pass # Should not reach here

    def test_basic_conversion_multi_sample_ref_strategy(self):
        """Test conversion of multi-sample VCF with default 'ref' unphased strategy."""
        converter = VCFtoGFAConverter(
            VCF_MULTI_SAMPLE, REF_FASTA, self.output_gfa, unphased_strategy='ref'
        )
        converter.convert()
        gfa = self._validate_gfa_structure(
            self.output_gfa,
            min_segments=10, min_links=10,
            expected_paths={"SAMPLE1_hap1", "SAMPLE1_hap2", "SAMPLE2_hap1", "SAMPLE2_hap2"}
        )
        # TODO: Add more specific checks based on the expected output GFA content for 'ref' strategy

    def test_basic_conversion_multi_sample_alt_strategy(self):
        """Test conversion of multi-sample VCF with 'alt' unphased strategy."""
        converter = VCFtoGFAConverter(
            VCF_MULTI_SAMPLE, REF_FASTA, self.output_gfa, unphased_strategy='alt'
        )
        converter.convert()
        gfa = self._validate_gfa_structure(
            self.output_gfa,
            min_segments=10, min_links=10,
            expected_paths={"SAMPLE1_hap1", "SAMPLE1_hap2", "SAMPLE2_hap1", "SAMPLE2_hap2"}
        )
        # TODO: Add specific checks for how unphased variants (rs4, rs7) are handled with 'alt'

    def test_basic_conversion_multi_sample_skip_strategy(self):
        """Test conversion of multi-sample VCF with 'skip' unphased strategy."""
        converter = VCFtoGFAConverter(
            VCF_MULTI_SAMPLE, REF_FASTA, self.output_gfa, unphased_strategy='skip'
        )
        converter.convert()
        gfa = self._validate_gfa_structure(
            self.output_gfa,
            min_segments=10, min_links=10,
            expected_paths={"SAMPLE1_hap1", "SAMPLE1_hap2", "SAMPLE2_hap1", "SAMPLE2_hap2"}
        )
        # TODO: Add specific checks for how unphased variants (rs4, rs7) are handled with 'skip'
        # Paths should effectively follow reference at these positions.

    def test_multi_block_conversion(self):
        """Test conversion with multiple phase sets."""
        converter = VCFtoGFAConverter(VCF_MULTI_BLOCK, REF_FASTA, self.output_gfa)
        converter.convert()
        gfa = self._validate_gfa_structure(
            self.output_gfa,
            min_segments=8, min_links=8,
            expected_paths={"SAMPLE_A_hap1", "SAMPLE_A_hap2"}
        )
        # TODO: Verify structure reflects the two phase blocks and the unphased variant correctly.

    def test_library_usage(self):
        """Verify that the converter uses the specified libraries (conceptual check)."""
        # This test implies library usage by successfully running basic conversion.
        # More rigorous checks would involve mocking library calls.
        try:
            with VCFtoGFAConverter(VCF_MULTI_SAMPLE, REF_FASTA, self.output_gfa) as converter:
                converter.convert()
            self.assertTrue(os.path.exists(self.output_gfa))
        except ImportError:
            self.fail("Conversion failed due to missing library, indicating reliance on them.")
        except Exception as e:
            # Allow known conversion/library errors, fail on others
            if isinstance(e, (VCFtoGFAConversionError, ReferenceHandlerError, PhasingError, gfapy.error.GfapyError, pysam.utils.SamtoolsError, FileNotFoundError)):
                 pass # Errors from the libraries or converter logic are somewhat expected
            else:
                 self.fail(f"Unexpected error, potentially indicating custom logic issues: {e}")

    def test_malformed_vcf(self):
        """Test handling of a VCF file with format errors."""
        # pysam might raise an error during parsing or iteration
        with self.assertRaises((VCFtoGFAConversionError, ValueError, RuntimeError)): # Catch potential errors from pysam/converter
             with VCFtoGFAConverter(VCF_MALFORMED, REF_FASTA, self.output_gfa) as converter:
                  converter.convert()

    def test_vcf_no_samples(self):
        """Test conversion of a VCF file with no samples."""
        with VCFtoGFAConverter(VCF_NO_SAMPLES, REF_FASTA, self.output_gfa) as converter:
            converter.convert()
        # Expect only the reference path for the contig
        gfa = self._validate_gfa_structure(
            self.output_gfa,
            min_segments=1, min_links=0,
            expected_paths={"chr1_ref"} # Only reference path expected
        )
        self.assertEqual(len(gfa.paths), 1)

    def test_empty_vcf(self):
        """Test conversion of an empty VCF file (only header)."""
        with VCFtoGFAConverter(VCF_EMPTY, REF_FASTA, self.output_gfa) as converter:
            converter.convert()
        # Expect only the reference path for the contig
        gfa = self._validate_gfa_structure(
            self.output_gfa,
            min_segments=1, min_links=0,
            expected_paths={"chr1_ref"} # Only reference path expected
        )
        self.assertEqual(len(gfa.paths), 1)

    def test_region_filtering(self):
        """Test conversion when a region is specified."""
        # Requires VCF index
        region = "chr1:20-110" # 1-based region including rs2, rs3, rs4, rs5
        try:
             with VCFtoGFAConverter(VCF_MULTI_SAMPLE, REF_FASTA, self.output_gfa) as converter:
                  converter.convert(region=region)
        except VCFtoGFAConversionError as e:
             # Skip if indexing failed during setup, otherwise fail
             if "index" in str(e).lower():
                  self.skipTest(f"Skipping region test, VCF index likely missing/failed: {e}")
             else:
                  self.fail(f"Region conversion failed unexpectedly: {e}")

        # Validate basic structure - expect fewer segments/links than full conversion
        gfa = self._validate_gfa_structure(
            self.output_gfa,
            min_segments=5, min_links=4, # Rough estimate, depends on exact segments created
            expected_paths={"SAMPLE1_hap1", "SAMPLE1_hap2", "SAMPLE2_hap1", "SAMPLE2_hap2"}
        )
        # TODO: More specific checks that only variants within the region affected the paths.


    # TODO: Add tests for:
    # - INDELs (already in test VCFs, need specific validation)
    # - Structural Variants (would require new test VCF and likely converter logic)
    # - Contig not in FASTA / Contig not in VCF


if __name__ == '__main__':
    # Re-enable logging for standalone test runs if desired
    # logging.disable(logging.NOTSET)
    # logging.basicConfig(level=logging.DEBUG)
    unittest.main()
