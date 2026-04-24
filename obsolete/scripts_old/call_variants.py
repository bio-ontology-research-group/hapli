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
from typing import List, Dict, Any, Set


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


def get_reference_chromosomes(reference_path: str) -> Set[str]:
    """
    Extract chromosome names from reference FASTA file.
    """
    chromosomes = set()
    
    with open(reference_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract chromosome name (everything after > until first space)
                chrom_name = line[1:].split()[0]
                chromosomes.add(chrom_name)
    
    return chromosomes


def get_pangenome_paths_from_vg(vg_path: str, vg_docker: str, mount_dir: str) -> Set[str]:
    """
    Extract path names from VG file using vg paths command.
    """
    paths = set()
    
    try:
        # Get relative path for Docker
        rel_vg = str(Path(vg_path).relative_to(mount_dir))
        
        # Run vg paths command
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{mount_dir}:/data",
            "-w", "/data",
            vg_docker,
            "vg", "paths", "-v", rel_vg, "-L"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        for line in result.stdout.strip().split('\n'):
            if line.strip():
                paths.add(line.strip())
                
    except subprocess.CalledProcessError as e:
        print(f"Warning: Could not extract paths from VG file: {e}")
        print(f"stderr: {e.stderr}")
    
    return paths


def get_reference_paths_from_pangenome(vg_path: str, vg_docker: str, mount_dir: str) -> Set[str]:
    """
    Extract reference path names from VG file.
    Reference paths typically follow pattern: chr*#0#chr*
    """
    all_paths = get_pangenome_paths_from_vg(vg_path, vg_docker, mount_dir)
    ref_paths = set()
    
    for path in all_paths:
        # Reference paths don't contain sample names and follow chr*#0#chr* pattern
        if not path.startswith('sample_') and '#0#chr' in path:
            ref_paths.add(path)
    
    return ref_paths


def extract_sample_names_from_vg(vg_path: str, vg_docker: str, mount_dir: str) -> List[str]:
    """
    Extract sample names from VG file by looking at path names.
    """
    all_paths = get_pangenome_paths_from_vg(vg_path, vg_docker, mount_dir)
    sample_names = set()
    
    for path in all_paths:
        # Extract sample name from path like "sample_000_hap1#0#chr1#0"
        if path.startswith('sample_') and '_hap' in path:
            # Extract everything before _hap
            sample_base = path.split('_hap')[0]
            sample_names.add(sample_base)
    
    return sorted(list(sample_names))


def check_chromosome_coverage(reference_path: str, vg_path: str, vg_docker: str, mount_dir: str) -> None:
    """
    Check if all chromosomes from reference are present in pangenome.
    """
    print("\n=== CHROMOSOME COVERAGE ANALYSIS ===")
    
    # Get chromosomes from reference
    ref_chroms = get_reference_chromosomes(reference_path)
    print(f"Reference chromosomes ({len(ref_chroms)}): {sorted(ref_chroms)}")
    
    # Get all paths from pangenome
    all_paths = get_pangenome_paths_from_vg(vg_path, vg_docker, mount_dir)
    print(f"All pangenome paths ({len(all_paths)}): {sorted(all_paths)}")
    
    # Get reference paths from pangenome
    ref_paths = get_reference_paths_from_pangenome(vg_path, vg_docker, mount_dir)
    print(f"Reference paths in pangenome ({len(ref_paths)}): {sorted(ref_paths)}")
    
    # Extract chromosome names from reference paths
    pangenome_chroms = set()
    for path in ref_paths:
        # Extract chromosome from path like "chr1#0#chr1"
        if '#0#' in path:
            chrom = path.split('#0#')[-1]
            pangenome_chroms.add(chrom)
    
    print(f"Chromosomes in pangenome ({len(pangenome_chroms)}): {sorted(pangenome_chroms)}")
    
    # Check coverage
    missing_chroms = ref_chroms - pangenome_chroms
    extra_chroms = pangenome_chroms - ref_chroms
    
    if missing_chroms:
        print(f"\n❌ MISSING chromosomes in pangenome: {sorted(missing_chroms)}")
    else:
        print(f"\n✅ All reference chromosomes found in pangenome")
    
    if extra_chroms:
        print(f"ℹ️  Extra chromosomes in pangenome: {sorted(extra_chroms)}")
    
    print("=" * 40)


def create_config(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Create the Snakemake config dictionary from command line arguments.
    """
    # Resolve paths
    pangenome_vg = str(Path(args.pangenome_vg).resolve())
    reference = str(Path(args.reference).resolve())
    output_dir = str(Path(args.output).resolve())
    
    # Verify input files exist
    if not Path(pangenome_vg).exists():
        raise FileNotFoundError(f"Pangenome VG file not found: {pangenome_vg}")
    if not Path(reference).exists():
        raise FileNotFoundError(f"Reference file not found: {reference}")
    
    # Find common mount point for Docker
    all_paths = [pangenome_vg, reference, output_dir]
    mount_dir = find_common_mount_point(all_paths)
    
    # Check chromosome coverage
    check_chromosome_coverage(reference, pangenome_vg, args.vg_docker, mount_dir)
    
    # Extract sample names from VG if not provided
    if args.samples:
        sample_names = args.samples
    else:
        sample_names = extract_sample_names_from_vg(pangenome_vg, args.vg_docker, mount_dir)
        if not sample_names:
            raise ValueError("No sample names found in VG file. Please specify --samples explicitly.")
    
    print(f"Sample names: {sample_names}")
    
    print(f"Pangenome VG: {pangenome_vg}")
    print(f"Reference: {reference}")
    print(f"Output directory: {output_dir}")
    print(f"Docker mount point: {mount_dir}")
    
    # Create config
    config = {
        "pangenome_vg": pangenome_vg,
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
  %(prog)s -v pangenome.vg -r reference.fa -o vcf_output

  # Call variants for specific samples
  %(prog)s -v pangenome.vg -r reference.fa -o vcf_output --samples sample_001 sample_002

  # Dry run to see what would be executed
  %(prog)s -v pangenome.vg -r reference.fa -o vcf_output --dry-run

  # Use more cores
  %(prog)s -v pangenome.vg -r reference.fa -o vcf_output --cores 8

  # Check chromosome coverage without running variant calling
  %(prog)s -v pangenome.vg -r reference.fa -o vcf_output --check-only
        """
    )
    
    # Input files
    parser.add_argument(
        "-v", "--pangenome-vg",
        required=True,
        help="Pangenome VG file"
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
        help="Sample names to call variants for (if not specified, extract from VG)"
    )
    
    # Execution options
    parser.add_argument(
        "--cores",
        type=int,
        default=8,
        help="Number of CPU cores for Snakemake (default: 8)"
    )
    
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads for individual tools (default: 8)"
    )
    
    parser.add_argument(
        "--vg-docker",
        default="quay.io/vgteam/vg:v1.65.0",
        help="Docker image for vg tools (default: quay.io/vgteam/vg:v1.65.0)"
    )
    
    # Diagnostic options
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Only check chromosome coverage, don't run variant calling"
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
        
        # If only checking, exit here
        if args.check_only:
            print("\nChromosome coverage check completed.")
            sys.exit(0)
        
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
