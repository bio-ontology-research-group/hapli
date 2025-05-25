#!/usr/bin/env python3
"""
Variant calling wrapper script for pangenomes.

This script provides a command-line interface for calling variants
from pangenome graphs using vg tools to generate phased VCF files.
"""

import argparse
import os
import sys
import yaml
import subprocess
from pathlib import Path
from typing import List, Dict, Any


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


def extract_sample_names_from_gfa(gfa_path: str) -> List[str]:
    """
    Extract sample names from GFA file by looking at path names.
    """
    sample_names = set()
    
    with open(gfa_path, 'r') as f:
        for line in f:
            if line.startswith('P\t'):  # Path line
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    path_name = parts[1]
                    # Extract sample name from path like "sample_000_hap1"
                    if '_hap' in path_name:
                        sample_base = path_name.split('_hap')[0]
                        sample_names.add(sample_base)
    
    return sorted(list(sample_names))


def create_config(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Create the Snakemake config dictionary from command line arguments.
    """
    # Resolve paths
    pangenome_vg = str(Path(args.pangenome_vg).resolve())
    pangenome_gfa = str(Path(args.pangenome_gfa).resolve())
    reference = str(Path(args.reference).resolve())
    output_dir = str(Path(args.output).resolve())
    
    # Verify input files exist
    if not Path(pangenome_vg).exists():
        raise FileNotFoundError(f"Pangenome VG file not found: {pangenome_vg}")
    if not Path(pangenome_gfa).exists():
        raise FileNotFoundError(f"Pangenome GFA file not found: {pangenome_gfa}")
    if not Path(reference).exists():
        raise FileNotFoundError(f"Reference file not found: {reference}")
    
    # Extract sample names from GFA if not provided
    if args.samples:
        sample_names = args.samples
    else:
        sample_names = extract_sample_names_from_gfa(pangenome_gfa)
        if not sample_names:
            raise ValueError("No sample names found in GFA file. Please specify --samples explicitly.")
    
    print(f"Sample names: {sample_names}")
    
    # Find common mount point for Docker
    all_paths = [pangenome_vg, pangenome_gfa, reference, output_dir]
    mount_dir = find_common_mount_point(all_paths)
    
    print(f"Pangenome VG: {pangenome_vg}")
    print(f"Pangenome GFA: {pangenome_gfa}")
    print(f"Reference: {reference}")
    print(f"Output directory: {output_dir}")
    print(f"Docker mount point: {mount_dir}")
    
    # Create config
    config = {
        "pangenome_vg": pangenome_vg,
        "pangenome_gfa": pangenome_gfa,
        "reference": reference,
        "output_dir": output_dir,
        "mount_dir": mount_dir,
        "sample_names": sample_names,
        "threads": args.threads,
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
        description="Call variants from pangenome graphs to generate phased VCF files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Call variants for all samples in pangenome
  %(prog)s -v pangenome.vg -g pangenome.gfa -r reference.fa -o vcf_output

  # Call variants for specific samples
  %(prog)s -v pangenome.vg -g pangenome.gfa -r reference.fa -o vcf_output --samples sample_001 sample_002

  # Dry run to see what would be executed
  %(prog)s -v pangenome.vg -g pangenome.gfa -r reference.fa -o vcf_output --dry-run

  # Use more cores
  %(prog)s -v pangenome.vg -g pangenome.gfa -r reference.fa -o vcf_output --cores 8
        """
    )
    
    # Input files
    parser.add_argument(
        "-v", "--pangenome-vg",
        required=True,
        help="Pangenome VG file"
    )
    
    parser.add_argument(
        "-g", "--pangenome-gfa",
        required=True,
        help="Pangenome GFA file"
    )
    
    parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Reference genome FASTA file"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory for VCF files"
    )
    
    # Sample selection
    parser.add_argument(
        "--samples",
        nargs="+",
        help="Sample names to call variants for (if not specified, extract from GFA)"
    )
    
    # Execution options
    parser.add_argument(
        "--cores",
        type=int,
        default=8, # Changed default from 4 to 8
        help="Number of CPU cores for Snakemake (default: 8)"
    )
    
    parser.add_argument(
        "--threads",
        type=int,
        default=8, # Changed default from 4 to 8
        help="Number of threads for individual tools (default: 8)"
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
        workflow_dir = script_dir.parent / "workflows" / "variant_calling"
        
        if not workflow_dir.exists():
            raise FileNotFoundError(f"Workflow directory not found: {workflow_dir}")
        
        # Write config file
        config_path = workflow_dir / "config.yaml"
        write_config_file(config, str(config_path))
        
        print(f"Config written to: {config_path}")
        
        # Run Snakemake
        exit_code = run_snakemake(str(workflow_dir), str(config_path), args)
        
        if exit_code == 0:
            print(f"\nVariant calling completed successfully!")
            print(f"Output files in: {config['output_dir']}")
        else:
            print(f"\nVariant calling failed with exit code: {exit_code}")
        
        sys.exit(exit_code)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
