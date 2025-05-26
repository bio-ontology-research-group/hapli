#!/usr/bin/env python3
"""
Pangenome construction wrapper script.

This script provides a command-line interface for building pangenomes
from haplotype FASTA files using Snakemake and Cactus tools.
"""

import argparse
import os
import sys
import glob
import yaml
import subprocess
from pathlib import Path
from typing import List, Dict, Any


def find_haplotype_files(sample_dirs: List[str], patterns: List[str] = None) -> List[str]:
    """
    Find haplotype FASTA files in the given sample directories.
    
    Args:
        sample_dirs: List of directories to search
        patterns: List of filename patterns to match
    
    Returns:
        List of absolute paths to haplotype FASTA files
    """
    if patterns is None:
        patterns = ["*_hap*.fasta", "*.fasta"]
    
    haplotype_files = []
    
    for sample_dir in sample_dirs:
        sample_path = Path(sample_dir).resolve()
        if not sample_path.exists():
            print(f"Warning: Sample directory does not exist: {sample_path}")
            continue
            
        print(f"Searching in: {sample_path}")
        
        for pattern in patterns:
            files = list(sample_path.glob(pattern))
            if files:
                print(f"  Found {len(files)} files matching '{pattern}':")
                for f in sorted(files):
                    print(f"    {f}")
                    haplotype_files.append(str(f))
                break
    
    return sorted(list(set(haplotype_files)))  # Remove duplicates and sort


def find_common_mount_point(paths: List[str]) -> str:
    """
    Find a common parent directory that contains all the given paths.
    This will be used as the Docker mount point.
    """
    if not paths:
        return os.getcwd()
    
    # Convert to Path objects and get absolute paths
    path_objects = [Path(p).resolve() for p in paths]
    
    # Find common parent
    common_parts = path_objects[0].parts
    for path_obj in path_objects[1:]:
        # Find common prefix
        new_common = []
        for i, (a, b) in enumerate(zip(common_parts, path_obj.parts)):
            if a == b:
                new_common.append(a)
            else:
                break
        common_parts = tuple(new_common)
    
    if common_parts:
        return str(Path(*common_parts))
    else:
        # Fallback to root
        return "/"


def create_config(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Create the Snakemake config dictionary from command line arguments.
    """
    # Find haplotype files
    if args.sample_dirs:
        haplotype_files = find_haplotype_files(args.sample_dirs, args.patterns)
    elif args.haplotype_files:
        haplotype_files = [str(Path(f).resolve()) for f in args.haplotype_files]
    else:
        raise ValueError("Must provide either --sample-dirs or --haplotype-files")
    
    if not haplotype_files:
        raise ValueError("No haplotype files found!")
    
    print(f"Found {len(haplotype_files)} haplotype files")
    
    # Resolve paths
    reference = str(Path(args.reference).resolve())
    output_dir = str(Path(args.output).resolve())
    
    # Find common mount point for Docker
    all_paths = [reference, output_dir] + haplotype_files
    mount_dir = find_common_mount_point(all_paths)
    
    print(f"Reference: {reference}")
    print(f"Output directory: {output_dir}")
    print(f"Docker mount point: {mount_dir}")
    
    # Verify reference exists
    if not Path(reference).exists():
        raise FileNotFoundError(f"Reference file not found: {reference}")
    
    # Create config
    config = {
        "reference": reference,
        "haplotype_fastas": haplotype_files,
        "output_dir": output_dir,
        "mount_dir": mount_dir,
        "threads": args.threads,
        "cactus_docker": args.cactus_docker,
        "vg_docker": args.vg_docker
    }
    
    return config


def write_config_file(config: Dict[str, Any], config_path: str) -> None:
    """Write the config to a YAML file."""
    with open(config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)


def run_snakemake(workflow_dir: str, config_path: str, args: argparse.Namespace) -> int:
    """
    Run Snakemake with the generated config.
    
    Returns:
        Exit code from Snakemake
    """
    cmd = [
        "snakemake",
        "--snakefile", os.path.join(workflow_dir, "Snakefile"),
        "--configfile", config_path,
        "--cores", str(args.cores),
    ]
    
    if args.dry_run:
        cmd.append("--dry-run")
    
    if args.verbose:
        cmd.extend(["--verbose"])
    
    if args.force:
        cmd.extend(["--forceall"])
    
    # Add any additional snakemake arguments
    if args.snakemake_args:
        cmd.extend(args.snakemake_args)
    
    print(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, cwd=workflow_dir)
        return result.returncode
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        return 1
    except Exception as e:
        print(f"Error running Snakemake: {e}")
        return 1


def main():
    parser = argparse.ArgumentParser(
        description="Build pangenomes from haplotype FASTA files using Cactus",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Build pangenome from sample directories
  %(prog)s -r reference.fa -s sample_001 sample_002 sample_003 -o pangenome_output

  # Build pangenome from specific files
  %(prog)s -r reference.fa -f sample1_hap1.fa sample1_hap2.fa sample2_hap1.fa -o output

  # Dry run to see what would be executed
  %(prog)s -r reference.fa -s sample_* -o output --dry-run

  # Use more cores and custom Docker images
  %(prog)s -r reference.fa -s sample_* -o output --cores 8 --cactus-docker custom/cactus:latest
        """
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-s", "--sample-dirs",
        nargs="+",
        help="Directories containing haplotype FASTA files"
    )
    input_group.add_argument(
        "-f", "--haplotype-files",
        nargs="+",
        help="Specific haplotype FASTA files to include"
    )
    
    parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Reference genome FASTA file"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory for pangenome files"
    )
    
    # Search patterns
    parser.add_argument(
        "-p", "--patterns",
        nargs="+",
        default=["*_hap*.fasta", "*.fasta"],
        help="Filename patterns to search for in sample directories (default: *_hap*.fasta *.fasta)"
    )
    
    # Execution options
    parser.add_argument(
        "--cores",
        type=int,
        default=4,
        help="Number of CPU cores for Snakemake (default: 4)"
    )
    
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for individual tools (default: 4)"
    )
    
    parser.add_argument(
        "--cactus-docker",
        default="quay.io/comparative-genomics-toolkit/cactus:v2.9.8",
        help="Docker image for Cactus tools (default: quay.io/comparative-genomics-toolkit/cactus:v2.9.8)"
    )
    
    parser.add_argument(
        "--vg-docker",
        default="quay.io/vgteam/vg:v1.65.0",
        help="Docker image for vg tools (default: quay.io/vgteam/vg:v1.65.0)"
    )
    
    # Snakemake options
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be executed without running"
    )
    
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-execution of all rules"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose output"
    )
    
    parser.add_argument(
        "--snakemake-args",
        nargs=argparse.REMAINDER,
        help="Additional arguments to pass to Snakemake"
    )
    
    args = parser.parse_args()
    
    try:
        # Create config
        config = create_config(args)
        
        # Find workflow directory
        script_dir = Path(__file__).parent.resolve()
        workflow_dir = script_dir.parent / "workflows" / "pangenome"
        
        if not workflow_dir.exists():
            raise FileNotFoundError(f"Workflow directory not found: {workflow_dir}")
        
        # Write config file
        config_path = workflow_dir / "config.yaml"
        write_config_file(config, str(config_path))
        
        print(f"Config written to: {config_path}")
        
        # Run Snakemake
        exit_code = run_snakemake(str(workflow_dir), str(config_path), args)
        
        if exit_code == 0:
            print(f"\nPangenome construction completed successfully!")
            print(f"Output files in: {config['output_dir']}")
        else:
            print(f"\nPangenome construction failed with exit code: {exit_code}")
        
        sys.exit(exit_code)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
