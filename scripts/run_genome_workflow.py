#!/usr/bin/env python3
"""
Main script to run the genome analysis workflow (downloading data, converting VCF to GFA).
Handles state management for resuming failed workflows.
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import time
import shutil # Import shutil for which
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from tqdm import tqdm # Import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('genome_workflow.log') # Log file name changed
    ]
)
logger = logging.getLogger(__name__)

# --- Prerequisite Check ---
def check_prerequisites(tools: List[str]) -> bool:
    """Check if required external tools are installed using shutil.which."""
    missing_tools = []
    for tool in tools:
        tool_path = shutil.which(tool)
        if tool_path:
            logger.info(f"Found {tool} at {tool_path}")
        else:
            logger.warning(f"{tool} not found in PATH.")
            missing_tools.append(tool)

    if missing_tools:
        logger.error(f"Missing required tools: {', '.join(missing_tools)}")
        logger.error("Please install these tools and ensure they are in your PATH accessible to this script.")
        # Log the PATH known to the script for debugging
        script_path_env = os.environ.get('PATH', 'PATH environment variable not set.')
        logger.error(f"Script's PATH: {script_path_env}")
        # Add specific installation instructions if possible
        if 'tabix' in missing_tools or 'bgzip' in missing_tools:
             logger.error("  'tabix' and 'bgzip' are part of HTSlib (https://github.com/samtools/htslib)")
        if 'samtools' in missing_tools:
             logger.error("  'samtools' can be installed via package managers or from https://github.com/samtools/samtools")
        if 'vg' in missing_tools:
             logger.error("  'vg' can be installed from https://github.com/vgteam/vg")
        if 'bcftools' in missing_tools:
             logger.error("  'bcftools' can be installed via package managers or from https://github.com/samtools/bcftools")
        return False
    return True


# --- Workflow Step Base Classes ---
class WorkflowStep:
    """Base class for workflow steps."""
    def __init__(self, name: str, description: str):
        self.name = name
        self.description = description
        self.start_time = None
        self.end_time = None
        self.success = None
        self.skipped = False
        self.outputs = [] # Expected output files/directories

    def should_run(self, state_file: Dict) -> bool:
        """Determine if this step should run based on state file and output existence."""
        # Check if step was already completed successfully
        step_state = state_file.get(self.name, {}) # Get state for this specific step
        if step_state.get('success', False):
            # Verify outputs still exist
            outputs_from_state = step_state.get('outputs', [])
            # Use the outputs defined in the step instance for checking existence
            outputs_to_check = self.get_expected_outputs()
            if not outputs_to_check: # If no outputs defined, assume success means skippable
                 logger.info(f"Skipping {self.name} - already completed (no specific outputs to check)")
                 self.skipped = True
                 return False
            if all(os.path.exists(out) for out in outputs_to_check):
                logger.info(f"Skipping {self.name} - already completed and outputs exist")
                self.skipped = True
                return False
            else:
                 missing = [out for out in outputs_to_check if not os.path.exists(out)]
                 logger.warning(f"Running {self.name} again - previously completed but expected outputs are missing: {missing}")

        return True

    def run(self, context: Dict[str, Any]) -> bool:
        """Run the step. Should be implemented by subclasses. Context can pass data between steps."""
        raise NotImplementedError("Subclasses must implement run()")

    def get_expected_outputs(self) -> List[str]:
        """Return the list of expected output files/directories for this step."""
        # Default implementation returns the static list. Subclasses can override.
        return self.outputs

    def get_state(self) -> Dict:
        """Get the state of this step for the state file."""
        # Use get_expected_outputs to ensure the saved state reflects the actual check
        start_iso = datetime.fromtimestamp(self.start_time).isoformat() if self.start_time else None
        end_iso = datetime.fromtimestamp(self.end_time).isoformat() if self.end_time else None
        return {
            'success': self.success,
            'start_time': start_iso,
            'end_time': end_iso,
            'outputs': self.get_expected_outputs(),
            'skipped': self.skipped,
            'status': 'completed' if self.success else ('skipped' if self.skipped else 'failed') # Add explicit status
        }

    def execute(self, state_file: Dict, context: Dict[str, Any]) -> bool:
        """Execute the step and record timing."""
        if not self.should_run(state_file):
            # Even if skipped, update context if needed (e.g., file paths)
            self.update_context(context)
            return True

        logger.info(f"--- Running Step: {self.name} ---")
        logger.info(f"Description: {self.description}")
        self.start_time = time.time()
        self.end_time = None # Reset end time
        self.success = None # Reset success status
        self.skipped = False # Reset skipped status

        try:
            self.success = self.run(context)
            # Update context after successful run
            if self.success:
                self.update_context(context)
        except Exception as e:
            logger.error(f"Error in step {self.name}: {e}", exc_info=True) # Add traceback
            self.success = False

        self.end_time = time.time()
        duration = self.end_time - self.start_time if self.start_time and self.end_time else 0

        if self.success:
            logger.info(f"Step {self.name} completed successfully in {duration:.2f} seconds")
        else:
            logger.error(f"Step {self.name} failed after {duration:.2f} seconds")
        logger.info(f"--- Finished Step: {self.name} ---")

        return self.success

    def update_context(self, context: Dict[str, Any]):
        """Update the shared context dictionary after the step runs (or is skipped)."""
        # Default implementation does nothing. Subclasses can override.
        pass


class CommandStep(WorkflowStep):
    """A workflow step that runs an external command."""
    def __init__(self, name: str, description: str, command: List[str], outputs: List[str]):
        super().__init__(name, description)
        self.command = command
        self.outputs = outputs # Static list of expected outputs

    def run(self, context: Dict[str, Any]) -> bool:
        """Run the command."""
        try:
            # Replace placeholders in command if needed (example)
            # formatted_command = [str(arg).format(**context) for arg in self.command]
            formatted_command = self.command # Assuming command is fully formed for now
            logger.info(f"Running command: {' '.join(formatted_command)}")
            result = subprocess.run(formatted_command, check=True, text=True, capture_output=True)
            logger.debug(f"Command stdout: {result.stdout}")
            if result.stderr:
                 logger.debug(f"Command stderr: {result.stderr}") # Log stderr as debug

            # Verify outputs were created
            missing_outputs = [out for out in self.get_expected_outputs() if not os.path.exists(out)]
            if missing_outputs:
                logger.error(f"Missing expected outputs after command execution: {', '.join(missing_outputs)}")
                return False

            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed with exit code {e.returncode}: {' '.join(e.cmd)}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"Error running command {self.name}: {e}", exc_info=True)
            return False

class PythonScriptStep(WorkflowStep):
    """A workflow step that runs a Python script."""
    def __init__(self, name: str, description: str, script_path: str, args: List[str], outputs: List[str]):
        super().__init__(name, description)
        self.script_path = script_path
        self.args = args
        self.outputs = outputs # Static list of expected outputs

    def run(self, context: Dict[str, Any]) -> bool:
        """Run the Python script."""
        try:
            # Replace placeholders in args if needed
            formatted_args = [str(arg).format(**context) for arg in self.args]
            # Filter out empty strings that might result from conditional formatting
            formatted_args = [arg for arg in formatted_args if arg]

            cmd = [sys.executable, self.script_path] + formatted_args
            logger.info(f"Running Python script: {' '.join(cmd)}")

            result = subprocess.run(cmd, check=True, text=True, capture_output=True)
            logger.debug(f"Script stdout: {result.stdout}")
            if result.stderr:
                 logger.debug(f"Script stderr: {result.stderr}") # Log stderr as debug

            # Verify outputs were created
            missing_outputs = [out for out in self.get_expected_outputs() if not os.path.exists(out)]
            if missing_outputs:
                # If output is a directory, check if it was created, even if empty
                dir_outputs = [out for out in missing_outputs if out in self.outputs and Path(out).suffix == '']
                file_outputs = [out for out in missing_outputs if out not in dir_outputs]

                all_missing = True
                for do in dir_outputs:
                    if os.path.isdir(do):
                         logger.warning(f"Directory output '{do}' exists but was flagged missing (might be empty). Considering it present.")
                         missing_outputs.remove(do) # Remove from missing list
                    else:
                         logger.error(f"Expected output directory not found: {do}")

                if file_outputs: # Only error if file outputs are missing
                     logger.error(f"Missing expected file outputs after script execution: {', '.join(file_outputs)}")
                     return False
                elif not missing_outputs: # All were directories that now exist
                     return True
                else: # Some other error state?
                     logger.error(f"Script finished but failed output validation.")
                     return False

            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Script failed with exit code {e.returncode}: {' '.join(e.cmd)}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            return False
        except KeyError as e:
             logger.error(f"Missing context variable '{e}' needed for script arguments.")
             return False
        except Exception as e:
            logger.error(f"Error running script {self.name}: {e}", exc_info=True)
            return False

class IndexVcfStep(WorkflowStep):
    """Indexes VCF files using tabix."""
    def __init__(self, name: str, description: str, vcf_context_key: str):
        super().__init__(name, description)
        self.vcf_context_key = vcf_context_key # Key in context holding list of VCF paths
        self._outputs = [] # Dynamically determined

    def run(self, context: Dict[str, Any]) -> bool:
        """Find .vcf.gz files and index them."""
        vcf_files = context.get(self.vcf_context_key)
        if not vcf_files:
            logger.warning(f"No VCF files found in context key '{self.vcf_context_key}'. Skipping indexing.")
            # If no files, is it an error? Or just nothing to do? Assume success.
            return True

        if not isinstance(vcf_files, list):
             vcf_files = [vcf_files] # Handle single file case

        logger.info(f"Found {len(vcf_files)} VCF files to potentially index from context key '{self.vcf_context_key}'")
        self._outputs = [] # Reset outputs for this run
        all_indexed = True
        processed_vcf_paths = [] # Keep track of the final .vcf.gz paths

        for vcf_file_str in vcf_files:
            vcf_file = Path(vcf_file_str)
            if not vcf_file.exists():
                logger.warning(f"VCF file not found: {vcf_file}. Skipping indexing.")
                continue # Skip this file, but don't fail the whole step yet

            # Check if it's bgzipped
            is_gzipped = vcf_file.name.endswith(".gz")
            target_vcf = vcf_file
            bgzip_success = True

            if not is_gzipped:
                logger.warning(f"VCF file {vcf_file} is not gzipped. Attempting to bgzip.")
                target_vcf = vcf_file.with_suffix(vcf_file.suffix + '.gz')
                # Check if bgzipped version already exists from a previous run
                if target_vcf.exists():
                     logger.info(f"Using existing bgzipped file: {target_vcf}")
                else:
                    try:
                        # Use shutil.which to find bgzip
                        bgzip_exe = shutil.which('bgzip')
                        if not bgzip_exe:
                             logger.error("bgzip command not found in PATH.")
                             all_indexed = False
                             break # Fail the step if bgzip isn't available and needed

                        bgzip_cmd = [bgzip_exe, '-f', str(vcf_file)] # -f to force overwrite if original exists? Be careful.
                        logger.info(f"Running: {' '.join(bgzip_cmd)}")
                        subprocess.run(bgzip_cmd, check=True, capture_output=True, text=True)
                        logger.info(f"Successfully created bgzipped file: {target_vcf}")
                        # Original file is removed by bgzip by default
                    except subprocess.CalledProcessError as e:
                        logger.error(f"bgzip failed for {vcf_file}: {e.stderr}")
                        bgzip_success = False
                        all_indexed = False
                        continue # Skip indexing this file
                    except Exception as e:
                         logger.error(f"Error running bgzip for {vcf_file}: {e}")
                         bgzip_success = False
                         all_indexed = False
                         continue

            if not bgzip_success:
                 continue

            # Index the (potentially newly) bgzipped file
            index_file = Path(str(target_vcf) + ".tbi")
            self._outputs.append(str(index_file)) # Add expected index file to outputs

            # Check if index already exists AND if it's older than the VCF
            reindex_needed = True
            if index_file.exists():
                try:
                    vcf_mtime = target_vcf.stat().st_mtime
                    idx_mtime = index_file.stat().st_mtime
                    if idx_mtime >= vcf_mtime:
                        logger.info(f"Index file {index_file.name} is up-to-date. Skipping indexing for {target_vcf.name}.")
                        reindex_needed = False
                        processed_vcf_paths.append(str(target_vcf)) # Add to list even if skipped
                    else:
                        logger.warning(f"Index file {index_file.name} is older than VCF file {target_vcf.name}. Re-indexing.")
                        # Remove the old index before re-running tabix
                        index_file.unlink()
                except OSError as e:
                    logger.warning(f"Could not check/remove old index file {index_file.name}: {e}. Will attempt re-indexing.")

            if reindex_needed:
                try:
                    # Use shutil.which to find tabix
                    tabix_exe = shutil.which('tabix')
                    if not tabix_exe:
                         logger.error("tabix command not found in PATH.")
                         all_indexed = False
                         break # Fail the step if tabix isn't available

                    tabix_cmd = [tabix_exe, '-p', 'vcf', str(target_vcf)]
                    logger.info(f"Running: {' '.join(tabix_cmd)}")
                    result = subprocess.run(tabix_cmd, check=True, capture_output=True, text=True)
                    logger.debug(f"tabix stdout: {result.stdout}")
                    if result.stderr:
                        logger.debug(f"tabix stderr: {result.stderr}")
                    if not index_file.exists():
                         logger.error(f"tabix command completed but index file {index_file} was not created.")
                         all_indexed = False
                    else:
                         logger.info(f"Successfully indexed {target_vcf.name}")
                         processed_vcf_paths.append(str(target_vcf)) # Add successfully indexed file
                except subprocess.CalledProcessError as e:
                    logger.error(f"tabix failed for {target_vcf.name} with exit code {e.returncode}: {' '.join(e.cmd)}")
                    logger.error(f"STDOUT: {e.stdout}")
                    logger.error(f"STDERR: {e.stderr}")
                    all_indexed = False
                except Exception as e:
                    logger.error(f"Error running tabix for {target_vcf.name}: {e}", exc_info=True)
                    all_indexed = False

        # Update context with the list of VCF files that are now confirmed to be gzipped
        # Only include files that were successfully processed (indexed or skipped because up-to-date)
        context[self.vcf_context_key] = processed_vcf_paths
        logger.debug(f"Updated context key '{self.vcf_context_key}' with processed VCF paths: {processed_vcf_paths}")

        return all_indexed

    def get_expected_outputs(self) -> List[str]:
        """Dynamically determine the expected .tbi files."""
        # Outputs are determined during run based on input VCFs
        return self._outputs

    def update_context(self, context: Dict[str, Any]):
        """Updates the context with the paths to the indexed VCF files (now guaranteed .gz)."""
        # The run method already updates the context list in place
        pass


class WorkflowManager:
    """Manages the execution of a sequence of workflow steps."""
    def __init__(self, state_file_path: str):
        self.state_file_path = state_file_path
        self.steps: List[WorkflowStep] = []
        self.state: Dict[str, Any] = {}
        self.context: Dict[str, Any] = {} # Shared context
        self.interrupted = False

        # Register signal handlers - Removed for simplicity, can be added back if needed
        # signal.signal(signal.SIGINT, self._handle_interrupt)
        # signal.signal(signal.SIGTERM, self._handle_interrupt)

        # Load state if it exists
        self._load_state()

    def _handle_interrupt(self, sig, frame):
        """Handle interrupt signals."""
        if not self.interrupted: # Prevent multiple messages
            logger.warning(f"Received signal {sig}. Workflow will attempt to stop gracefully after the current step.")
            self.interrupted = True

    def _load_state(self):
        """Load the state file if it exists."""
        if os.path.exists(self.state_file_path):
            try:
                with open(self.state_file_path, 'r') as f:
                    self.state = json.load(f)
                logger.info(f"Loaded workflow state from {self.state_file_path}")
            except json.JSONDecodeError as e:
                logger.error(f"Error decoding state file {self.state_file_path}: {e}. Starting fresh.")
                self.state = {"steps": {}} # Initialize with empty steps dict
            except Exception as e:
                logger.error(f"Error loading state file: {e}. Starting fresh.")
                self.state = {"steps": {}} # Initialize with empty steps dict
        else:
            logger.info("No existing state file found. Starting from the beginning.")
            self.state = {"steps": {}} # Initialize with empty steps dict

    def _save_state(self):
        """Save the current state to the state file."""
        try:
            state_dir = os.path.dirname(self.state_file_path)
            if state_dir and not os.path.exists(state_dir):
                os.makedirs(state_dir, exist_ok=True)

            # Use the original path provided, ensuring directory exists
            path_to_save = self.state_file_path
            with open(path_to_save, 'w') as f:
                json.dump(self.state, f, indent=2, sort_keys=True)
            logger.debug(f"Saved workflow state to {path_to_save}") # Debug level
        except Exception as e:
            logger.error(f"Error saving state file: {e}")

    def add_step(self, step: WorkflowStep):
        """Add a step to the workflow."""
        self.steps.append(step)

    def run(self) -> bool:
        """Run all steps in the workflow."""
        total_steps = len(self.steps)
        completed_steps = 0
        overall_success = True

        logger.info(f"Starting workflow with {total_steps} steps. State file: {self.state_file_path}")
        logger.info(f"Shared context initialized as: {self.context}")

        # Create a progress bar for the entire workflow
        with tqdm(total=total_steps, desc="Workflow Progress", unit="step", disable=None) as pbar:
            for i, step in enumerate(self.steps):
                pbar.set_description(f"Workflow Progress ({step.name})")
                if self.interrupted:
                    logger.warning("Workflow interrupted by user. Stopping.")
                    overall_success = False
                    break

                # Pass the state specific to this step to execute method
                step_state_data = self.state.get("steps", {}).get(step.name, {})
                step_success = step.execute(step_state_data, self.context)

                # Update state with step results immediately after execution
                if "steps" not in self.state:
                     self.state["steps"] = {}
                self.state["steps"][step.name] = step.get_state()
                self._save_state() # Save state after each step

                # Update progress bar only if the step wasn't skipped initially
                if not step.skipped:
                     completed_steps += 1 # Increment completed only if run/failed
                pbar.update(1)


                if not step_success and not step.skipped:
                    overall_success = False
                    logger.error(f"Step {step.name} failed. Stopping workflow.")
                    break
                elif step.skipped:
                     logger.debug(f"Step {step.name} was skipped, continuing.")
                else:
                     logger.debug(f"Step {step.name} succeeded, continuing.")


        # Final status logging
        executed_count = sum(1 for step in self.steps if step.start_time is not None and not step.skipped)
        skipped_count = sum(1 for step in self.steps if step.skipped)
        not_run_count = total_steps - executed_count - skipped_count

        if overall_success and not self.interrupted:
            logger.info(f"Workflow completed successfully.")
        elif self.interrupted:
            logger.warning(f"Workflow was interrupted.")
        else:
            logger.error(f"Workflow failed.")

        logger.info(f"Summary: {executed_count} steps executed, {skipped_count} skipped, {not_run_count} not run.")

        self._generate_report(total_steps)
        return overall_success

    def _generate_report(self, total_steps: int):
        """Generate a summary report of the workflow execution."""
        report_path = "workflow_report.txt" # Keep simple name for now
        executed_steps = 0
        successful_steps = 0
        failed_steps = 0
        skipped_steps = 0

        try:
            with open(report_path, 'w') as f:
                f.write("=== Genome Workflow Execution Report ===\n") # Updated title
                f.write(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"State File: {self.state_file_path}\n")

                total_duration = 0.0

                f.write("\n--- Step Details ---\n")
                step_states = self.state.get("steps", {})
                for step in self.steps:
                    state = step_states.get(step.name, {}) # Get state from loaded/saved state
                    status = "□ Not Run"
                    duration_str = ""

                    if state:
                        # Use the explicit status saved in the state if available
                        step_status = state.get('status', 'pending')

                        if step_status == "skipped":
                            status = "⏭ Skipped"
                            skipped_steps += 1
                        elif step_status == "completed":
                            status = "✓ Success"
                            successful_steps += 1
                            executed_steps += 1 # Count successful as executed
                        elif step_status == "failed":
                            status = "✗ Failed"
                            failed_steps += 1
                            executed_steps += 1 # Count failed as executed
                        elif step_status == "running": # Should not happen in final report, but handle
                             status = "? Running (Error?)"
                             executed_steps += 1
                        else: # Pending or other
                             status = "□ Not Run"


                        # Calculate duration from state file times if available
                        start_t = state.get('start_time')
                        end_t = state.get('end_time')
                        if start_t and end_t and step_status != "skipped":
                             try:
                                 # Parse ISO format times
                                 start_dt = datetime.fromisoformat(start_t)
                                 end_dt = datetime.fromisoformat(end_t)
                                 secs = (end_dt - start_dt).total_seconds()
                                 total_duration += secs
                                 mins, secs = divmod(secs, 60)
                                 duration_str = f" ({int(mins)}m {secs:.1f}s)"
                             except (ValueError, TypeError):
                                 duration_str = " (duration error)"
                    else:
                         # Step wasn't even reached in the state
                         pass


                    f.write(f"\n{status}: {step.name}{duration_str}\n")
                    f.write(f"  Desc: {step.description}\n")

                    outputs = state.get('outputs', step.get_expected_outputs()) # Show expected if not run
                    if outputs:
                        f.write("  Outputs:\n")
                        max_outputs_to_show = 5
                        for i, output in enumerate(outputs):
                            if i >= max_outputs_to_show:
                                f.write(f"   ... and {len(outputs) - max_outputs_to_show} more\n")
                                break
                            exists_str = ""
                            if state: # Only check existence if step was attempted/skipped
                                exists_str = "✓" if os.path.exists(output) else "✗"
                            f.write(f"   - {exists_str} {output}\n")

                f.write("\n--- Summary ---\n")
                f.write(f"Total Steps: {total_steps}\n")
                f.write(f"Successful:  {successful_steps}\n")
                f.write(f"Failed:      {failed_steps}\n")
                f.write(f"Skipped:     {skipped_steps}\n")
                f.write(f"Not Run:     {total_steps - (successful_steps + failed_steps + skipped_steps)}\n")

                mins, secs = divmod(total_duration, 60)
                hours, mins = divmod(mins, 60)
                f.write(f"\nTotal Execution Time (for run steps): {int(hours)}h {int(mins)}m {secs:.1f}s\n")

            logger.info(f"Workflow report generated at {report_path}")
        except Exception as e:
            logger.error(f"Error generating report: {e}", exc_info=True)


# Helper Step Classes for Dynamic VCF handling
class DownloadVcfStep(PythonScriptStep):
    """Custom step to download VCFs and update context."""
    def __init__(self, name, description, script_path, args, outputs, context_key_phased, context_key_unphased):
        super().__init__(name, description, script_path, args, outputs)
        self.context_key_phased = context_key_phased
        self.context_key_unphased = context_key_unphased

    def update_context(self, context: Dict[str, Any]):
        """Update context with paths to downloaded VCF files."""
        # This relies on the download script having predictable output names or
        # ideally, the download script should output a manifest file that this
        # step reads to update the context accurately.
        # For now, using glob based on expected patterns from download_1000g_vcf.py
        output_dir = Path(context['vcf_dir'])
        phased_files = list(output_dir.glob("phased_*.vcf.gz"))
        unphased_files = list(output_dir.glob("unphased_*.vcf.gz"))

        context[self.context_key_phased] = [str(f) for f in phased_files if f.exists()]
        context[self.context_key_unphased] = [str(f) for f in unphased_files if f.exists()]
        # Combine them for the indexing step
        context['all_vcf_files'] = context[self.context_key_phased] + context[self.context_key_unphased]

        logger.info(f"Context updated: {self.context_key_phased} = {context[self.context_key_phased]}")
        logger.info(f"Context updated: {self.context_key_unphased} = {context[self.context_key_unphased]}")
        logger.info(f"Context updated: all_vcf_files = {context['all_vcf_files']}")


class SeparateVcfToGfaStep(PythonScriptStep):
    """Runs VCF to GFA conversion in 'separate' mode."""
    def __init__(self, name, description, script_path, args_template, outputs_dir_context_key, vcf_context_key):
        super().__init__(name, description, script_path, [], []) # No static args/outputs initially
        self.args_template_base = args_template
        self.outputs_dir_context_key = outputs_dir_context_key
        self.vcf_context_key = vcf_context_key
        self._outputs = [] # Dynamically determined output files

    def run(self, context: Dict[str, Any]) -> bool:
        """Constructs args and runs the script for separate conversion."""
        vcf_files = context.get(self.vcf_context_key)
        output_dir = context.get(self.outputs_dir_context_key)
        if not vcf_files:
             logger.warning(f"Skipping {self.name}: No VCF files found in context key '{self.vcf_context_key}'.")
             return True # Success, nothing to do
        if not output_dir:
            logger.error(f"Error in {self.name}: Output directory context key '{self.outputs_dir_context_key}' not found.")
            return False

        if not isinstance(vcf_files, list):
            vcf_files = [vcf_files]

        # Build the command line arguments
        try:
            # Format base arguments
            formatted_args = [str(arg).format(**context) for arg in self.args_template_base]
            # Filter out empty strings from conditional args
            formatted_args = [arg for arg in formatted_args if arg]

            # Add VCF file arguments
            for vcf_file in vcf_files:
                formatted_args.extend(['--vcf', str(vcf_file)])
            # Add mode
            formatted_args.extend(['--mode', 'separate'])
            # Add output directory
            formatted_args.extend(['--output', str(output_dir)])

        except KeyError as e:
            logger.error(f"Missing context variable '{e}' needed for script arguments for {self.name}.")
            return False

        # Set the formatted args for the parent class run method
        self.args = formatted_args # Override the instance's args for this run

        # Determine expected outputs before running
        self._outputs = []
        output_dir_path = Path(output_dir)
        for vcf_path_str in vcf_files:
            vcf_path = Path(vcf_path_str)
            base_name = vcf_path.name
            if base_name.endswith(".vcf.gz"):
                base_name = base_name[:-len(".vcf.gz")]
            elif base_name.endswith(".vcf"):
                 base_name = base_name[:-len(".vcf")]
            output_gfa_path = output_dir_path / f"{base_name}.gfa"
            self._outputs.append(str(output_gfa_path))

        # Call the parent run method
        success = super().run(context)
        return success

    def get_expected_outputs(self) -> List[str]:
        """Returns the list of GFA files expected to be created."""
        # Outputs are determined during run based on input VCFs and output dir
        return self._outputs

    def update_context(self, context: Dict[str, Any]):
        """Adds the list of generated GFA files to the context."""
        # Add the list of generated GFA files if they exist
        # Note: The script itself should handle resume logic for individual files
        context['separate_gfa_files'] = [str(o) for o in self._outputs if Path(o).exists()]


class JointVcfToGfaStep(PythonScriptStep):
    """Runs VCF to GFA conversion in 'joint' mode."""
    def __init__(self, name, description, script_path, args_template, output_file_context_key, vcf_context_key):
        super().__init__(name, description, script_path, [], []) # Output is dynamic
        self.args_template_base = args_template
        self.output_file_context_key = output_file_context_key
        self.vcf_context_key = vcf_context_key
        self._output_file = None # Dynamically determined

    def run(self, context: Dict[str, Any]) -> bool:
        """Constructs args and runs the script for joint conversion."""
        vcf_files = context.get(self.vcf_context_key)
        output_file = context.get(self.output_file_context_key)
        if not vcf_files:
             logger.warning(f"Skipping {self.name}: No VCF files found in context key '{self.vcf_context_key}'.")
             return True # Success, nothing to do
        if not output_file:
            logger.error(f"Error in {self.name}: Output file context key '{self.output_file_context_key}' not found.")
            return False

        if not isinstance(vcf_files, list):
            vcf_files = [vcf_files]

        self._output_file = str(output_file) # Store for get_expected_outputs

        # Build the command line arguments
        try:
            # Format base arguments
            formatted_args = [str(arg).format(**context) for arg in self.args_template_base]
            # Filter out empty strings from conditional args
            formatted_args = [arg for arg in formatted_args if arg]

            # Add VCF file arguments
            for vcf_file in vcf_files:
                formatted_args.extend(['--vcf', str(vcf_file)])
            # Add mode
            formatted_args.extend(['--mode', 'joint'])
            # Add output file
            formatted_args.extend(['--output', str(output_file)])

        except KeyError as e:
            logger.error(f"Missing context variable '{e}' needed for script arguments for {self.name}.")
            return False

        # Set the formatted args for the parent class run method
        self.args = formatted_args # Override the instance's args for this run

        # Call the parent run method
        success = super().run(context)
        return success

    def get_expected_outputs(self) -> List[str]:
        """Returns the single joint GFA file path."""
        return [self._output_file] if self._output_file else []

    def update_context(self, context: Dict[str, Any]):
        """Adds the path to the generated joint GFA file to the context."""
        # Add the joint GFA file path if it exists
        if self._output_file and Path(self._output_file).exists():
             context['joint_gfa_file'] = self._output_file


def create_workflow(args):
    """Create the workflow with all necessary steps."""
    # Define base directories relative to a root project directory if possible
    # Or use absolute paths based on args
    base_output_dir = Path(args.output_dir).resolve()
    ref_dir = base_output_dir / "reference_hs37d5"
    annot_dir = base_output_dir / "annotation_grch37"
    vcf_dir = base_output_dir / "vcf"
    graph_dir = base_output_dir / "graphs_hs37d5"
    state_file = base_output_dir / "workflow_state_genome.json" # Changed state file name

    # --- Workflow Manager ---
    workflow = WorkflowManager(state_file_path=str(state_file))

    # --- Initial Context ---
    # These can be overridden by loaded state or step outputs
    workflow.context = {
        "base_output_dir": str(base_output_dir),
        "ref_dir": str(ref_dir),
        "annot_dir": str(annot_dir),
        "vcf_dir": str(vcf_dir),
        "graph_dir": str(graph_dir),
        "reference_fasta_gz": str(ref_dir / "hs37d5.fa.gz"), # Path to downloaded gz
        "reference_fasta": str(ref_dir / "hs37d5.fa"),      # Path to extracted fa
        "reference_fai": str(ref_dir / "hs37d5.fa.fai"),    # Path to index
        # Updated GFF3 filenames
        "annotation_file_gz": str(annot_dir / "Homo_sapiens.GRCh37.87.chromosome.22.gff3.gz"),
        "annotation_file": str(annot_dir / "Homo_sapiens.GRCh37.87.chromosome.22.gff3"),
        "vcf_files_phased": [], # To be populated by download step
        "vcf_files_unphased": [], # To be populated by download step
        "all_vcf_files": [], # Combined list
        "region": args.region, # Optional region for VCF->GFA
        "vg_mem_gb": args.vg_mem,
        "vg_threads": args.vg_threads,
        "skip_existing": args.skip_existing, # Pass skip flag to download steps
        "num_phased": args.num_phased, # Pass counts to download step
        "num_unphased": args.num_unphased,
        "sample_ids": args.sample_ids # Pass sample IDs
    }

    # Create necessary directories early
    for d_key in ["ref_dir", "annot_dir", "vcf_dir", "graph_dir"]:
        Path(workflow.context[d_key]).mkdir(parents=True, exist_ok=True)

    # --- Step 1: Download Reference Genome (hs37d5) ---
    # Outputs now include the index file which is generated by the script
    workflow.add_step(PythonScriptStep(
        name="download_reference",
        description="Download hs37d5 reference genome FASTA and generate index.",
        script_path="scripts/download_hs37d5.py",
        args=[
            "--output-dir", "{ref_dir}",
            "--extract", # Always extract for downstream use
            "--skip-existing" if args.skip_existing else "" # Conditionally add skip flag
        ],
        outputs=[workflow.context["reference_fasta"], workflow.context["reference_fai"]],
    ))

    # --- Step 2: Download Annotation (GRCh37 GFF3 - Chr22) ---
    # Updated description and outputs list
    workflow.add_step(PythonScriptStep(
        name="download_annotation",
        description="Download GRCh37 GFF3 annotation file (Chr 22).",
        script_path="scripts/download_grch37_gff3.py",
        args=[
            "--output-dir", "{annot_dir}",
            "--extract", # Always extract for downstream use
            "--skip-existing" if args.skip_existing else ""
        ],
        outputs=[workflow.context["annotation_file"]], # Expect only the extracted file
    ))

    # --- Step 3: Download VCF files (Phased and Unphased) ---
    vcf_download_args = [
            "--output-dir", "{vcf_dir}",
            "--phased", "{num_phased}",
            "--unphased", "{num_unphased}",
            # "--extract", # Ignored by download script now
            "--skip-existing" if args.skip_existing else ""
        ]
    # Add sample IDs if provided
    if args.sample_ids:
         vcf_download_args.append("--sample-ids")
         vcf_download_args.extend(args.sample_ids)

    workflow.add_step(DownloadVcfStep(
        name="download_vcf",
        description="Download VCF files from 1000 Genomes Project (as .vcf.gz)",
        script_path="scripts/download_1000g_vcf.py",
        args=vcf_download_args,
        outputs=[workflow.context["vcf_dir"]], # Output is the directory, check existence
        context_key_phased="vcf_files_phased",
        context_key_unphased="vcf_files_unphased"
    ))

    # --- Step 4: Index VCF files ---
    # This step now uses the 'all_vcf_files' key populated by DownloadVcfStep
    index_step = IndexVcfStep(
        name="index_vcf",
        description="Index downloaded VCF files using tabix (and bgzip if needed).",
        vcf_context_key="all_vcf_files" # Index all downloaded VCFs
        # Outputs determined dynamically by the step
    )
    workflow.add_step(index_step)


    # --- Step 5: Convert VCFs to GFA (separately) ---
    separate_graph_dir = graph_dir / "separate"
    workflow.context["separate_graph_dir"] = str(separate_graph_dir) # Add to context
    Path(separate_graph_dir).mkdir(parents=True, exist_ok=True) # Ensure dir exists

    # Define the base template for arguments, VCFs and output are added dynamically
    vcf_to_gfa_base_args = [
         "--reference", "{reference_fasta}",
         "--memory", "{vg_mem_gb}",
         "--threads", "{vg_threads}",
         # Conditionally add region
         "--region" if args.region else "",
         "{region}" if args.region else "",
         # Conditionally add skip flag (assuming vcf_to_gfa_converter supports it)
         # "--skip-existing" if args.skip_existing else ""
    ]

    workflow.add_step(SeparateVcfToGfaStep(
        name="vcf_to_gfa_separate",
        description="Convert VCF files to GFA (separately)",
        script_path="scripts/vcf_to_gfa_converter.py", # ASSUMES THIS SCRIPT EXISTS
        args_template=vcf_to_gfa_base_args,
        outputs_dir_context_key="separate_graph_dir",
        vcf_context_key="all_vcf_files" # Use the combined, indexed list
    ))

    # --- Step 6: Convert VCFs to GFA (jointly) ---
    joint_gfa_path = graph_dir / "joint_graph.gfa"
    workflow.context["joint_gfa_file_path"] = str(joint_gfa_path) # Add to context
    Path(joint_gfa_path).parent.mkdir(parents=True, exist_ok=True) # Ensure dir exists

    workflow.add_step(JointVcfToGfaStep(
        name="vcf_to_gfa_joint",
        description="Convert VCF files to GFA (jointly)",
        script_path="scripts/vcf_to_gfa_converter.py", # ASSUMES THIS SCRIPT EXISTS
        args_template=vcf_to_gfa_base_args,
        output_file_context_key="joint_gfa_file_path",
        vcf_context_key="all_vcf_files" # Use the combined, indexed list
    ))

    return workflow

def main():
    parser = argparse.ArgumentParser(description="Run the Genome Analysis Workflow (hs37d5/GRCh37).")
    parser.add_argument('--output-dir', '-o', default="workflow_output_hs37d5",
                        help='Base directory for all workflow outputs and state.')
    parser.add_argument('--num-phased', type=int, default=1,
                        help='Max number of phased VCF files/samples to download/use.')
    parser.add_argument('--num-unphased', type=int, default=1,
                        help='Max number of unphased VCF files/samples to download/use.')
    parser.add_argument('--sample-ids', type=str, nargs='+', default=None,
                        help='Specific sample IDs to download/extract (e.g., HG01383 NA12878). Overrides defaults in download script.')
    parser.add_argument('--region', type=str, default=None, # Changed default to None
                        help='Restrict VCF conversion to a specific region (e.g., "22", "1:1M-2M"). Recommended if VCFs span multiple contigs.')
    # Increased default memory for vg construct
    parser.add_argument('--vg-mem', type=int, default=32,
                        help='Memory limit (GB) for vg construct.')
    parser.add_argument('--vg-threads', type=int, default=os.cpu_count() or 8,
                        help='Number of threads for vg construct.')
    parser.add_argument('--skip-existing', action='store_true',
                        help='Skip download/processing steps if final output files already exist.')
    # Add arguments for specific VCF files if needed, instead of downloading

    args = parser.parse_args()

    # --- Prerequisite Check (Overall) ---
    # Check tools needed by any step upfront
    # Added samtools here
    all_prereqs = ['bcftools', 'tabix', 'bgzip', 'samtools', 'vg']
    if not check_prerequisites(all_prereqs):
        sys.exit(1)

    # Handle region argument
    if args.region:
         logger.info(f"Running VCF conversion restricted to region: {args.region}")
    else:
         logger.info("Running VCF conversion without region restriction (ensure VCFs/reference match contigs or vg may fail).")


    # --- Create and Run Workflow ---
    try:
        # Ensure the converter script exists before starting
        converter_script = Path("scripts/vcf_to_gfa_converter.py")
        if not converter_script.exists():
             # Check if the placeholder error message is in the file
             try:
                 with open(converter_script, 'r') as f:
                     content = f.read()
                     if "Variant/reference sequence mismatch" in content:
                          logger.critical(f"Required script '{converter_script}' contains only an error message. Please restore its code (e.g., from git hash 7cf2f5c).")
                     else:
                          logger.critical(f"Required script '{converter_script}' not found. Please restore it.")
             except Exception:
                 logger.critical(f"Required script '{converter_script}' not found or cannot be read. Please restore it.")
             sys.exit(2)

        workflow = create_workflow(args)
        success = workflow.run()
        sys.exit(0 if success else 1)
    except Exception as e:
        logger.critical(f"An unexpected critical error occurred in the workflow setup or execution: {e}", exc_info=True)
        sys.exit(2)


if __name__ == "__main__":
    main()
