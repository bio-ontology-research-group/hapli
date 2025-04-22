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
from typing import Dict, List, Optional, Tuple
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

class WorkflowStep:
    """Base class for workflow steps."""
    def __init__(self, name: str, description: str):
        self.name = name
        self.description = description
        self.start_time = None
        self.end_time = None
        self.success = None
        self.skipped = False
        self.outputs = []
    
    def should_run(self, state_file: Dict) -> bool:
        """Determine if this step should run based on state file and output existence."""
        # Check if step was already completed successfully
        if self.name in state_file and state_file[self.name].get('success', False):
            # Verify outputs still exist
            outputs = state_file[self.name].get('outputs', [])
            if all(os.path.exists(out) for out in outputs):
                logger.info(f"Skipping {self.name} - already completed")
                self.skipped = True
                return False
        return True
    
    def run(self) -> bool:
        """Run the step. Should be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement run()")
    
    def get_state(self) -> Dict:
        """Get the state of this step for the state file."""
        return {
            'success': self.success,
            'start_time': self.start_time,
            'end_time': self.end_time,
            'outputs': self.outputs,
            'skipped': self.skipped
        }
    
    def execute(self, state_file: Dict) -> bool:
        """Execute the step and record timing."""
        if not self.should_run(state_file):
            return True
        
        logger.info(f"Running step: {self.name} - {self.description}")
        self.start_time = time.time()
        
        try:
            self.success = self.run()
        except Exception as e:
            logger.error(f"Error in step {self.name}: {e}")
            self.success = False
        
        self.end_time = time.time()
        duration = self.end_time - self.start_time
        
        if self.success:
            logger.info(f"Step {self.name} completed successfully in {duration:.2f} seconds")
        else:
            logger.error(f"Step {self.name} failed after {duration:.2f} seconds")
        
        return self.success

class CommandStep(WorkflowStep):
    """A workflow step that runs a command."""
    def __init__(self, name: str, description: str, command: List[str], outputs: List[str]):
        super().__init__(name, description)
        self.command = command
        self.outputs = outputs
    
    def run(self) -> bool:
        """Run the command."""
        try:
            logger.info(f"Running command: {' '.join(self.command)}")
            result = subprocess.run(self.command, check=True, text=True, capture_output=True)
            logger.debug(f"Command output: {result.stdout}")
            
            # Verify outputs were created
            missing_outputs = [out for out in self.outputs if not os.path.exists(out)]
            if missing_outputs:
                logger.error(f"Missing expected outputs: {', '.join(missing_outputs)}")
                return False
                
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed with exit code {e.returncode}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"Error running command: {e}")
            return False

class PythonScriptStep(WorkflowStep):
    """A workflow step that runs a Python script."""
    def __init__(self, name: str, description: str, script_path: str, args: List[str], outputs: List[str]):
        super().__init__(name, description)
        self.script_path = script_path
        self.args = args
        self.outputs = outputs
    
    def run(self) -> bool:
        """Run the Python script."""
        try:
            cmd = [sys.executable, self.script_path] + self.args
            logger.info(f"Running Python script: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, check=True, text=True, capture_output=True)
            logger.debug(f"Script output: {result.stdout}")
            
            # Verify outputs were created
            missing_outputs = [out for out in self.outputs if not os.path.exists(out)]
            if missing_outputs:
                logger.error(f"Missing expected outputs: {', '.join(missing_outputs)}")
                return False
                
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Script failed with exit code {e.returncode}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"Error running script: {e}")
            return False

class WorkflowManager:
    """Manages the execution of a sequence of workflow steps."""
    def __init__(self, state_file_path: str):
        self.state_file_path = state_file_path
        self.steps = []
        self.state = {}
        self.interrupted = False
        
        # Register signal handlers
        signal.signal(signal.SIGINT, self._handle_interrupt)
        signal.signal(signal.SIGTERM, self._handle_interrupt)
        
        # Load state if it exists
        self._load_state()
    
    def _handle_interrupt(self, sig, frame):
        """Handle interrupt signals."""
        logger.warning("Received interrupt signal. Workflow will stop after current step.")
        self.interrupted = True
    
    def _load_state(self):
        """Load the state file if it exists."""
        if os.path.exists(self.state_file_path):
            try:
                with open(self.state_file_path, 'r') as f:
                    self.state = json.load(f)
                logger.info(f"Loaded workflow state from {self.state_file_path}")
            except Exception as e:
                logger.error(f"Error loading state file: {e}")
                self.state = {}
        else:
            logger.info("No existing state file found. Starting from the beginning.")
            self.state = {}
    
    def _save_state(self):
        """Save the current state to the state file."""
        try:
            # Handle case where state_file_path might be just a filename without directory
            if os.path.dirname(self.state_file_path) == '':
                # No directory specified, save in current directory
                path_to_save = os.path.join(os.getcwd(), self.state_file_path)
            else:
                path_to_save = self.state_file_path
                # Create directory if it doesn't exist
                os.makedirs(os.path.dirname(path_to_save), exist_ok=True)
            
            with open(path_to_save, 'w') as f:
                json.dump(self.state, f, indent=2)
            logger.info(f"Saved workflow state to {path_to_save}")
        except Exception as e:
            logger.error(f"Error saving state file: {e}")
    
    def add_step(self, step: WorkflowStep):
        """Add a step to the workflow."""
        self.steps.append(step)
    
    def run(self) -> bool:
        """Run all steps in the workflow."""
        total_steps = len(self.steps)
        completed_steps = 0
        success = True
        
        logger.info(f"Starting workflow with {total_steps} steps")
        
        # Create a progress bar for the entire workflow
        with tqdm(total=total_steps, desc="Workflow Progress", unit="step") as pbar:
            for step in self.steps:
                if self.interrupted:
                    logger.warning("Workflow interrupted by user. Stopping.")
                    break
                
                step_success = step.execute(self.state)
                
                # Update state with step results
                self.state[step.name] = step.get_state()
                self._save_state()
                
                completed_steps += 1
                pbar.update(1)
                
                if not step_success:
                    success = False
                    logger.error(f"Step {step.name} failed. Stopping workflow.")
                    break
        
        if success and not self.interrupted:
            logger.info(f"Workflow completed successfully. {completed_steps}/{total_steps} steps executed.")
        elif self.interrupted:
            logger.warning(f"Workflow was interrupted. {completed_steps}/{total_steps} steps executed.")
            success = False
        else:
            logger.error(f"Workflow failed. {completed_steps}/{total_steps} steps executed.")
        
        self._generate_report(completed_steps, total_steps)
        return success
    
    def _generate_report(self, completed_steps: int, total_steps: int):
        """Generate a summary report of the workflow execution."""
        report_path = "workflow_report.txt"
        
        try:
            with open(report_path, 'w') as f:
                f.write("=== GRCh38 Workflow Execution Report ===\n\n")
                f.write(f"Steps completed: {completed_steps}/{total_steps}\n")
                f.write("\nStep Details:\n")
                
                for step in self.steps:
                    state = self.state.get(step.name, {})
                    status = "✓ Success" if state.get('success', False) else "✗ Failed" if state else "□ Not Run"
                    
                    if state.get('skipped', False):
                        status = "⏭ Skipped"
                    
                    duration = ""
                    if state.get('start_time') and state.get('end_time'):
                        secs = state['end_time'] - state['start_time']
                        mins, secs = divmod(secs, 60)
                        duration = f" ({int(mins)}m {int(secs)}s)"
                    
                    f.write(f"\n{status}: {step.name}{duration}\n")
                    f.write(f"  Description: {step.description}\n")
                    
                    if 'outputs' in state and state['outputs']:
                        f.write("  Outputs:\n")
                        for output in state['outputs']:
                            exists = "✓" if os.path.exists(output) else "✗"
                            f.write(f"   - {exists} {output}\n")
            
            logger.info(f"Workflow report generated at {report_path}")
        except Exception as e:
            logger.error(f"Error generating report: {e}")

def create_workflow(args):
    """Create the workflow with all necessary steps."""
    workflow = WorkflowManager(args.state_file)
    
    # Create necessary directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, "reference"), exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, "annotation"), exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, "vcf"), exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, "graphs"), exist_ok=True)
    
    # Step 1: Download GRCh38 reference
    reference_dir = os.path.join(args.output_dir, "reference")
    reference_gz_file = os.path.join(reference_dir, 
                                   "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" 
                                   if args.source == "ensembl" else
                                   "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz")
    
    reference_file = os.path.join(reference_dir, 
                                 "Homo_sapiens.GRCh38.dna.primary_assembly.fa" 
                                 if args.source == "ensembl" else
                                 "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
    
    workflow.add_step(PythonScriptStep(
        name="download_reference",
        description="Download GRCh38 reference genome",
        script_path="scripts/download_grch38.py",
        args=["--output-dir", reference_dir, "--extract", "--source", args.source],
        outputs=[reference_file, reference_gz_file]  # Add both files as potential outputs
    ))
    
    # Step 2: Download GFF3 annotation
    annotation_dir = os.path.join(args.output_dir, "annotation")
    # Define both compressed and extracted filenames
    gff3_gz_file = os.path.join(annotation_dir, 
                              "Homo_sapiens.GRCh38.109.gff3.gz" 
                              if args.source == "ensembl" else
                              "GCA_000001405.15_GRCh38_genomic.gff.gz"
                              if args.source == "ncbi" else
                              "gencode.v44.annotation.gff3.gz")
    
    gff3_file = os.path.join(annotation_dir,
                            "Homo_sapiens.GRCh38.109.gff3" 
                            if args.source == "ensembl" else
                            "GCA_000001405.15_GRCh38_genomic.gff"
                            if args.source == "ncbi" else
                            "gencode.v44.annotation.gff3")
    
    workflow.add_step(PythonScriptStep(
        name="download_annotation",
        description="Download GFF3 annotation",
        script_path="scripts/download_grch38_gff3.py",
        args=["--output-dir", annotation_dir, "--extract", "--source", args.source],
        outputs=[gff3_file, gff3_gz_file]  # Add both files as potential outputs
    ))
    
    # Step 3: Convert reference to GFA
    reference_gfa = os.path.join(args.output_dir, "graphs", "GRCh38_reference.gfa")
    workflow.add_step(PythonScriptStep(
        name="reference_to_gfa",
        description="Convert reference FASTA to GFA",
        script_path="scripts/fasta_to_gfa.py",
        args=["--input", reference_file, "--output", reference_gfa],
        outputs=[reference_gfa]
    ))
    
    # Define consistent sample IDs to use across all VCF downloads
    sample_ids = ["HG01383", "NA12878", "NA20530"]  # Specific samples for all downloads
    
    # Step 4: Download VCF files
    vcf_dir = os.path.join(args.output_dir, "vcf")
    # We can't predict the exact filenames, so we'll just check the directory
    workflow.add_step(PythonScriptStep(
        name="download_vcf",
        description="Download VCF files from 1000 Genomes Project",
        script_path="scripts/download_1000g_vcf.py",
        args=[
            "--output-dir", vcf_dir,
            "--phased", str(args.phased_vcfs),
            "--unphased", str(args.unphased_vcfs),
            "--extract",
            "--force-unphased",  # Force download of unphased VCFs
            "--sample-ids"
        ] + sample_ids,
        outputs=[vcf_dir]  # Just check if the directory exists
    ))
    
    # Step 4b: Download Structural Variant VCF files
    sv_dir = os.path.join(args.output_dir, "structural_variants")
    os.makedirs(sv_dir, exist_ok=True)
    workflow.add_step(PythonScriptStep(
        name="download_sv_vcf",
        description="Download Structural Variant VCF files",
        script_path="scripts/download_1000g_sv.py",
        args=[
            "--output-dir", sv_dir,
            "--samples", "3",  # Download for 3 samples
            "--sample-ids"
        ] + sample_ids,  # Use the same samples as regular VCFs
        outputs=[sv_dir]  # Just check if the directory exists
    ))
    
    # Step 5: Convert VCFs to GFA (separately)
    graphs_dir = os.path.join(args.output_dir, "graphs")
    separate_gfa_dir = os.path.join(graphs_dir, "separate")
    os.makedirs(separate_gfa_dir, exist_ok=True)
    
    # We'll determine VCF files dynamically
    vcf_files = []
    if os.path.exists(vcf_dir):
        vcf_files = [os.path.join(vcf_dir, f) for f in os.listdir(vcf_dir) 
                    if f.endswith('.vcf') and not f.endswith('.gz')]
    
    if vcf_files:
        vcf_args = []
        for vcf in vcf_files[:2]:  # Use only first two VCFs for demonstration
            vcf_args.extend(["--vcf", vcf])
        
        workflow.add_step(PythonScriptStep(
            name="vcf_to_gfa_separate",
            description="Convert VCF files to GFA (separately)",
            script_path="scripts/vcf_to_gfa_converter.py",
            args=[
                *vcf_args,
                "--reference", reference_file,
                "--output", separate_gfa_dir,
                "--mode", "separate",
                "--region", args.region or "22"  # Default to chromosome 22
            ],
            outputs=[separate_gfa_dir]
        ))
        
        # Step 6: Convert VCFs to GFA (jointly)
        joint_gfa = os.path.join(graphs_dir, "joint_output.gfa")
        workflow.add_step(PythonScriptStep(
            name="vcf_to_gfa_joint",
            description="Convert VCF files to GFA (jointly)",
            script_path="scripts/vcf_to_gfa_converter.py",
            args=[
                *vcf_args,
                "--reference", reference_file,
                "--output", joint_gfa,
                "--mode", "joint",
                "--region", args.region or "22"  # Default to chromosome 22
            ],
            outputs=[joint_gfa]
        ))
    
    return workflow

def main():
    parser = argparse.ArgumentParser(description="Run the entire GRCh38 data processing workflow")
    parser.add_argument('--output-dir', type=str, default='data',
                        help='Base output directory for all data')
    parser.add_argument('--source', type=str, choices=['ensembl', 'ncbi'], default='ensembl',
                        help='Source for reference and annotation (ensembl or ncbi)')
    parser.add_argument('--phased-vcfs', type=int, default=2,
                        help='Number of phased VCF files to download (max 3)')
    parser.add_argument('--unphased-vcfs', type=int, default=3,
                        help='Number of unphased VCF files to download (max 3)')
    parser.add_argument('--region', type=str, default=None,
                        help='Restrict VCF conversion to region (e.g., "22")')
    parser.add_argument('--state-file', type=str, default='workflow_state.json',
                        help='File to store workflow state for resumption')
    args = parser.parse_args()
    
    # Create and run the workflow
    workflow = create_workflow(args)
    success = workflow.run()
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
