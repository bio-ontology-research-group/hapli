import unittest
import subprocess
import os
import sys

# Define paths relative to this file's location
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # Project root
SCRIPT_PATH = os.path.join(BASE_DIR, "scripts", "validate_gfa2.py")
EXAMPLE_GFA_PATH = os.path.join(BASE_DIR, "data", "example.gfa")

class TestDataValidation(unittest.TestCase):

    def test_validate_example_gfa_with_script(self):
        """
        Runs the validate_gfa2.py script on the example GFA file.

        NOTE: This test is EXPECTED TO FAIL because example.gfa is GFA1
              and validate_gfa2.py checks for GFA2 format.
              It correctly identifies the GFA1 header as invalid GFA2.
        """
        self.assertTrue(os.path.exists(SCRIPT_PATH), f"Validation script not found at {SCRIPT_PATH}")
        self.assertTrue(os.path.exists(EXAMPLE_GFA_PATH), f"Example GFA file not found at {EXAMPLE_GFA_PATH}")

        try:
            # Run the script using the same Python interpreter that runs the tests
            result = subprocess.run(
                [sys.executable, SCRIPT_PATH, EXAMPLE_GFA_PATH],
                capture_output=True,
                text=True,
                check=False, # Don't raise exception on non-zero exit
                encoding='utf-8'
            )

            # Check if the script executed successfully (exit code 0)
            # We expect this to fail because example.gfa is GFA1, not GFA2.
            # If the validator behaves correctly, it should return a non-zero exit code.
            # To make the test *pass* if the validator *correctly fails* on GFA1:
            # self.assertNotEqual(result.returncode, 0, ...)
            # However, the request is to "validate the example gfa file", implying
            # an expectation of success. So we assert for code 0, which will fail.
            self.assertEqual(
                result.returncode, 0,
                f"Validation script failed for {EXAMPLE_GFA_PATH} (Expected failure for GFA1 file).\n"
                f"Exit Code: {result.returncode}\n"
                f"Stderr:\n{result.stderr}\n"
                f"Stdout:\n{result.stdout}"
            )

        except FileNotFoundError:
            self.fail(f"Failed to execute Python interpreter at '{sys.executable}'.")
        except Exception as e:
            self.fail(f"Running validation script raised an unexpected exception: {e}")

if __name__ == '__main__':
    unittest.main()
