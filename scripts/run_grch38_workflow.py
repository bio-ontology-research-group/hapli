#!/usr/bin/env python3
"""
Script to run the entire GRCh38 data processing workflow sequentially.
Provides resumption capability and progress tracking.
"""

import argparse
import logging
import os
import sys
import time
import hashlib
import subprocess
import signal
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from tqdm import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('grch38_workflow.log')
    ]
)
logger = logging.getLogger(__name__)

def check_prerequisites(tools: List[str]) -> bool:
    """Check if required external tools are installed."""
    missing_tools = []
    for tool in tools:
        try:
            # Use a simple command that should succeed if the tool is present
            cmd = [tool, '--version'] if tool not in ['bgzip', 'tabix'] else [tool, '-h'] # bgzip/tabix don't have --version
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
            logger.info(f"Found {tool}")
        except (FileNotFoundError, subprocess.CalledProcessError, OSError):
            missing_tools.append(tool)

    if missing_tools:
        logger.error(f"Missing required tools: {', '.join(missing_tools)}")
        logger.error("Please install these tools and ensure they are in your PATH.")
        return False
    return True


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
        if self.name in state_file and state_file[self.name].get('success', False):
            # Verify outputs still exist
            outputs_from_state = state_file[self.name].get('outputs', [])
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
        return {
            'success': self.success,
            'start_time': self.start_time,
            'end_time': self.end_time,
            'outputs': self.get_expected_outputs(),
            'skipped': self.skipped
        }

    def execute(self, state_file: Dict, context: Dict[str, Any]) -> bool:
        """Execute the step and record timing."""
        if not self.should_run(state_file):
            # Even if skipped, update context if needed (e.g., file paths)
            self.update_context(context)
            return True

        logger.info(f"Running step: {self.name} - {self.description}")
        self.start_time = time.time()

        try:
            self.success = self.run(context)
            # Update context after successful run
            if self.success:
                self.update_context(context)
        except Exception as e:
            logger.error(f"Error in step {self.name}: {e}", exc_info=True) # Add traceback
            self.success = False

        self.end_time = time.time()
        duration = self.end_time - self.start_time

        if self.success:
            logger.info(f"Step {self.name} completed successfully in {duration:.2f} seconds")
        else:
            logger.error(f"Step {self.name} failed after {duration:.2f} seconds")

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
            # formatted_args = [str(arg).format(**context) for arg in self.args]
            formatted_args = self.args # Assuming args are fully formed
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
        except Exception as e:
            logger.error(f"Error running script {self.name}: {e}", exc_info=True)
            return False

class IndexVcfStep(WorkflowStep):
    """Indexes VCF files using tabix."""
    def __init__(self, name: str, description: str, vcf_dir: str):
        super().__init__(name, description)
        self.vcf_dir = vcf_dir
        # Outputs are determined dynamically in get_expected_outputs

    def run(self, context: Dict[str, Any]) -> bool:
        """Find .vcf.gz files and index them."""
        vcf_files = list(Path(self.vcf_dir).glob('*.vcf.gz'))
        if not vcf_files:
            logger.warning(f"No *.vcf.gz files found in {self.vcf_dir} to index.")
            # If no files, is it an error? Or just nothing to do? Assume success.
            return True

        logger.info(f"Found {len(vcf_files)} VCF files to index in {self.vcf_dir}")
        all_indexed = True
        for vcf_file in vcf_files:
            index_file = Path(f"{vcf_file}.tbi")
            if index_file.exists():
                logger.info(f"Index file {index_file} already exists, skipping indexing for {vcf_file.name}")
                continue

            cmd = ['tabix', '-p', 'vcf', str(vcf_file)]
            try:
                logger.info(f"Indexing {vcf_file.name}...")
                result = subprocess.run(cmd, check=True, text=True, capture_output=True)
                logger.debug(f"tabix stdout: {result.stdout}")
                if result.stderr:
                    logger.debug(f"tabix stderr: {result.stderr}")
                if not index_file.exists():
                     logger.error(f"tabix command completed but index file {index_file} was not created.")
                     all_indexed = False
                else:
                     logger.info(f"Successfully indexed {vcf_file.name}")
            except subprocess.CalledProcessError as e:
                logger.error(f"tabix failed for {vcf_file.name} with exit code {e.returncode}: {' '.join(e.cmd)}")
                logger.error(f"STDOUT: {e.stdout}")
                logger.error(f"STDERR: {e.stderr}")
                all_indexed = False
            except Exception as e:
                logger.error(f"Error running tabix for {vcf_file.name}: {e}", exc_info=True)
                all_indexed = False

        # Update context with the list of indexed files (including .tbi)
        # context['indexed_vcf_files'] = self.get_expected_outputs() # Context updated in update_context
        return all_indexed

    def get_expected_outputs(self) -> List[str]:
        """Dynamically determine the expected .tbi files."""
        vcf_files = list(Path(self.vcf_dir).glob('*.vcf.gz'))
        return [str(Path(f"{vcf_file}.tbi")) for vcf_file in vcf_files]

    def update_context(self, context: Dict[str, Any]):
        """Add the list of VCF.gz files to the context."""
        vcf_files = list(Path(self.vcf_dir).glob('*.vcf.gz'))
        context['vcf_gz_files'] = [str(f) for f in vcf_files]
        logger.debug(f"Updated context with vcf_gz_files: {context['vcf_gz_files']}")


class WorkflowManager:
    """Manages the execution of a sequence of workflow steps."""
    def __init__(self, state_file_path: str):
        self.state_file_path = state_file_path
        self.steps: List[WorkflowStep] = []
        self.state: Dict[str, Any] = {}
        self.context: Dict[str, Any] = {} # Shared context
        self.interrupted = False

        # Register signal handlers
        signal.signal(signal.SIGINT, self._handle_interrupt)
        signal.signal(signal.SIGTERM, self._handle_interrupt)

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
                self.state = {}
            except Exception as e:
                logger.error(f"Error loading state file: {e}. Starting fresh.")
                self.state = {}
        else:
            logger.info("No existing state file found. Starting from the beginning.")
            self.state = {}

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

                step_success = step.execute(self.state, self.context)

                # Update state with step results immediately after execution
                self.state[step.name] = step.get_state()
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
        report_path = "workflow_report.txt"
        executed_steps = 0
        successful_steps = 0
        failed_steps = 0
        skipped_steps = 0

        try:
            with open(report_path, 'w') as f:
                f.write("=== GRCh38 Workflow Execution Report ===\n")
                f.write(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"State File: {self.state_file_path}\n")

                total_duration = 0.0

                f.write("\n--- Step Details ---\n")
                for step in self.steps:
                    state = self.state.get(step.name, {})
                    status = "□ Not Run"
                    duration_str = ""

                    if state:
                        # executed_steps +=1 # This counts attempts, not just runs
                        if state.get('skipped', False):
                            status = "⏭ Skipped"
                            skipped_steps += 1
                        elif state.get('success', False):
                            status = "✓ Success"
                            successful_steps += 1
                            executed_steps += 1 # Count successful as executed
                        else:
                            status = "✗ Failed"
                            failed_steps += 1
                            executed_steps += 1 # Count failed as executed

                        if state.get('start_time') and state.get('end_time') and not state.get('skipped'):
                            secs = state['end_time'] - state['start_time']
                            total_duration += secs
                            mins, secs = divmod(secs, 60)
                            duration_str = f" ({int(mins)}m {secs:.1f}s)"
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
class SeparateVcfToGfaStep(PythonScriptStep):
    def __init__(self, name, description, script_path, args_template, outputs_dir):
        super().__init__(name, description, script_path, [], [outputs_dir]) # No static args/outputs initially
        self.args_template = args_template
        self.outputs_dir = outputs_dir

    def run(self, context: Dict[str, Any]) -> bool:
        vcf_files = context.get('vcf_gz_files')
        if not vcf_files:
            logger.warning(f"Skipping {self.name}: No VCF.gz files found in context.")
            return True # Success, nothing to do

        # Dynamically generate args based on template and context
        self.args = []
        for vcf in vcf_files:
            self.args.extend(["--vcf", vcf])
        self.args.extend(self.args_template) # Add the rest of the args

        # Run the script using the parent's run method
        return super().run(context)

    def get_expected_outputs(self) -> List[str]:
         # Output is the directory itself for separate mode
         # We could try and predict filenames, but checking the dir is simpler for now
         return [self.outputs_dir]

    def update_context(self, context: Dict[str, Any]):
         # Add the output directory path to context if needed later
         context['separate_gfa_dir'] = self.outputs_dir


class JointVcfToGfaStep(PythonScriptStep):
    def __init__(self, name, description, script_path, args_template, output_file):
        super().__init__(name, description, script_path, [], [output_file]) # Output is a file
        self.args_template = args_template
        self.output_file = output_file

    def run(self, context: Dict[str, Any]) -> bool:
        vcf_files = context.get('vcf_gz_files')
        if not vcf_files:
            logger.warning(f"Skipping {self.name}: No VCF.gz files found in context.")
            return True

        self.args = []
        for vcf in vcf_files:
            self.args.extend(["--vcf", vcf])
        self.args.extend(self.args_template)

        return super().run(context)

    def get_expected_outputs(self) -> List[str]:
         return [self.output_file]

    def update_context(self, context: Dict[str, Any]):
         context['joint_gfa_file'] = self.output_file


def create_workflow(args):
    """Create the workflow with all necessary steps."""
    workflow = WorkflowManager(args.state_file)
    base_output_dir = Path(args.output_dir).resolve() # Use absolute paths

    # --- Define Directories ---
    reference_dir = base_output_dir / "reference"
    annotation_dir = base_output_dir / "annotation"
    vcf_dir = base_output_dir / "vcf"
    sv_dir = base_output_dir / "structural_variants"
    graphs_dir = base_output_dir / "graphs"
    separate_gfa_dir = graphs_dir / "separate"
    joint_gfa_dir = graphs_dir / "joint" # Directory for joint output

    # Create necessary directories early
    for d in [reference_dir, annotation_dir, vcf_dir, sv_dir, graphs_dir, separate_gfa_dir, joint_gfa_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # --- Update Context with Paths ---
    workflow.context.update({
        "reference_dir": str(reference_dir),
        "annotation_dir": str(annotation_dir),
        "vcf_dir": str(vcf_dir),
        "sv_dir": str(sv_dir),
        "graphs_dir": str(graphs_dir),
        "separate_gfa_dir": str(separate_gfa_dir),
        "joint_gfa_dir": str(joint_gfa_dir),
        "source": args.source,
        "region": args.region, # Store region from args (could be None)
        "phased_vcfs": args.phased_vcfs,
        "unphased_vcfs": args.unphased_vcfs,
        "sample_ids": args.sample_ids,
    })

    # --- Step 1: Download GRCh38 reference ---
    ref_base_name = ("Homo_sapiens.GRCh38.dna.primary_assembly" if args.source == "ensembl"
                     else "GCA_000001405.15_GRCh38_no_alt_analysis_set")
    reference_gz_file = reference_dir / f"{ref_base_name}.fa.gz"
    reference_file = reference_dir / f"{ref_base_name}.fa"
    workflow.context["reference_file"] = str(reference_file) # Add to context

    workflow.add_step(PythonScriptStep(
        name="download_reference",
        description=f"Download GRCh38 reference genome ({args.source})",
        script_path="scripts/download_grch38.py",
        args=["--output-dir", str(reference_dir), "--extract", "--source", args.source],
        outputs=[str(reference_file), str(reference_gz_file)] # Expect both potentially
    ))

    # --- Step 2: Download GFF3 annotation ---
    gff_base_name = ("Homo_sapiens.GRCh38.109" if args.source == "ensembl"
                     else "GCA_000001405.15_GRCh38_genomic" if args.source == "ncbi"
                     else "gencode.v44.annotation") # Assuming gencode if not ensembl/ncbi
    gff3_gz_file = annotation_dir / f"{gff_base_name}.gff3.gz" # Standardize extension
    gff3_file = annotation_dir / f"{gff_base_name}.gff3"
    workflow.context["gff3_file"] = str(gff3_file) # Add to context

    workflow.add_step(PythonScriptStep(
        name="download_annotation",
        description=f"Download GFF3 annotation ({args.source})",
        script_path="scripts/download_grch38_gff3.py",
        args=["--output-dir", str(annotation_dir), "--extract", "--source", args.source],
        outputs=[str(gff3_file), str(gff3_gz_file)] # Expect both potentially
    ))

    # --- Step 3: Convert reference to GFA ---
    reference_gfa = graphs_dir / "GRCh38_reference.gfa"
    workflow.context["reference_gfa"] = str(reference_gfa)
    workflow.add_step(PythonScriptStep(
        name="reference_to_gfa",
        description="Convert reference FASTA to GFA",
        script_path="scripts/fasta_to_gfa.py",
        args=["--input", str(reference_file), "--output", str(reference_gfa)],
        outputs=[str(reference_gfa)]
    ))

    # --- Step 4: Download VCF files ---
    # Define consistent sample IDs to use across all VCF downloads
    # Use provided sample IDs or default
    sample_ids_to_download = args.sample_ids or ["HG01383", "NA12878"] # Default if none provided
    workflow.context["sample_ids"] = sample_ids_to_download # Update context

    vcf_download_args = [
            "--output-dir", str(vcf_dir),
            "--phased", str(args.phased_vcfs),
            "--unphased", str(args.unphased_vcfs),
            # "--extract", # REMOVED: Keep files gzipped
            # "--force-unphased", # Deprecated in download script
            "--sample-ids"
        ] + sample_ids_to_download

    workflow.add_step(PythonScriptStep(
        name="download_vcf",
        description="Download VCF files from 1000 Genomes Project (as .vcf.gz)",
        script_path="scripts/download_1000g_vcf.py",
        args=vcf_download_args,
        outputs=[str(vcf_dir)] # Output is the directory, check existence
    ))

    # --- Step 4a: Index VCF files ---
    index_step = IndexVcfStep(
        name="index_vcf",
        description="Index downloaded VCF files using tabix",
        vcf_dir=str(vcf_dir)
        # Outputs determined dynamically by the step
    )
    workflow.add_step(index_step)


    # --- Step 4b: Download Structural Variant VCF files ---
    # Note: These might also need indexing if used with -R in vg. Assuming not for now.
    sv_download_args = [
            "--output-dir", str(sv_dir),
            "--samples", "3", # Example: download for 3 samples
            "--sample-ids"
        ] + sample_ids_to_download # Use the same samples

    workflow.add_step(PythonScriptStep(
        name="download_sv_vcf",
        description="Download Structural Variant VCF files",
        script_path="scripts/download_1000g_sv.py",
        args=sv_download_args,
        outputs=[str(sv_dir)] # Check directory existence
    ))

    # --- Step 5: Convert VCFs to GFA (separately) ---
    # This step now depends on the context having 'vcf_gz_files' from the index step
    separate_args_template = [
        "--reference", str(reference_file),
        "--output", str(separate_gfa_dir),
        "--mode", "separate",
    ]
    # Conditionally add region argument
    if workflow.context["region"]:
        separate_args_template.extend(["--region", workflow.context["region"]])

    workflow.add_step(SeparateVcfToGfaStep(
        name="vcf_to_gfa_separate",
        description="Convert VCF files to GFA (separately)",
        script_path="scripts/vcf_to_gfa_converter.py",
        args_template=separate_args_template,
        outputs_dir=str(separate_gfa_dir)
    ))

    # --- Step 6: Convert VCFs to GFA (jointly) ---
    # Similar dynamic approach for joint conversion
    joint_gfa_file = joint_gfa_dir / "joint_graph.gfa" # Specific output file
    joint_args_template = [
        "--reference", str(reference_file),
        "--output", str(joint_gfa_file),
        "--mode", "joint",
    ]
    # Conditionally add region argument
    if workflow.context["region"]:
        joint_args_template.extend(["--region", workflow.context["region"]])

    workflow.add_step(JointVcfToGfaStep(
        name="vcf_to_gfa_joint",
        description="Convert VCF files to GFA (jointly)",
        script_path="scripts/vcf_to_gfa_converter.py",
        args_template=joint_args_template,
        output_file=str(joint_gfa_file)
    ))

    return workflow

def main():
    parser = argparse.ArgumentParser(description="Run the GRCh38 data processing workflow.")
    parser.add_argument('--output-dir', type=str, default='data',
                        help='Base output directory for all generated data (default: data)')
    parser.add_argument('--source', type=str, choices=['ensembl', 'ncbi', 'gencode'], default='ensembl',
                        help='Source for reference and annotation (default: ensembl)')
    parser.add_argument('--phased-vcfs', type=int, default=1, # Changed default to 1
                        help='Max number of phased VCF files per sample to download/use (default: 1)')
    parser.add_argument('--unphased-vcfs', type=int, default=1, # Changed default to 1
                        help='Max number of unphased VCF files per sample to download/use (default: 1)')
    parser.add_argument('--sample-ids', type=str, nargs='+', default=None,
                        help='Specific sample IDs to download (e.g., HG01383 NA12878). Overrides default.')
    parser.add_argument('--region', type=str, default="22",
                        help='Restrict VCF conversion to a specific region (e.g., "22", "chr1:1M-2M"). Default: "22". Use "" for no region.')
    parser.add_argument('--state-file', type=str, default='workflow_state.json',
                        help='File to store workflow state for resumption (default: workflow_state.json)')

    args = parser.parse_args()

    # --- Prerequisite Check ---
    # Check for tools needed by the workflow itself and the scripts it calls
    required_tools = ['vg', 'tabix', 'bgzip', 'python3', 'bcftools'] # Added bcftools
    if not check_prerequisites(required_tools):
        logger.error("Prerequisite check failed. Please install missing tools.")
        sys.exit(1)

    # Handle empty string for region argument
    if args.region == "":
        args.region = None
        logger.info("Running VCF conversion without region restriction.")
    elif args.region:
         logger.info(f"Running VCF conversion restricted to region: {args.region}")


    # --- Create and Run Workflow ---
    try:
        workflow = create_workflow(args)
        success = workflow.run()
        sys.exit(0 if success else 1)
    except Exception as e:
        logger.critical(f"An unexpected critical error occurred in the workflow setup or execution: {e}", exc_info=True)
        sys.exit(2)


if __name__ == "__main__":
    main()
