#!/usr/bin/env python3
"""
Master script for generating complete test datasets.

This script orchestrates the entire test data generation pipeline:
1. Generate reference genome
2. Generate variant samples
3. Build pangenome
4. Call variants to generate VCF files

Usage:
    python scripts/generate_test_data.py --preset small
    python scripts/generate_test_data.py --preset medium --output-dir data/medium_test
    python scripts/generate_test_data.py --custom --num-samples 5 --variants-per-haplotype 100
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import Dict, Any
import yaml


def setup_logging(verbose: bool = False) -> None:
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def run_command(cmd: list, description: str, cwd: Path = None) -> None:
    """Run a command and handle errors."""
    logging.info(f"Running: {description}")
    logging.debug(f"Command: {' '.join(map(str, cmd))}")
    
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            check=True,
            capture_output=True,
            text=True
        )
        
        if result.stdout:
            logging.debug(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            logging.debug(f"STDERR:\n{result.stderr}")
            
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {description}")
        logging.error(f"Exit code: {e.returncode}")
        logging.error(f"STDOUT:\n{e.stdout}")
        logging.error(f"STDERR:\n{e.stderr}")
        raise


def get_preset_config(preset: str) -> Dict[str, Any]:
    """Get configuration for predefined presets."""
    presets = {
        "tiny": {
            "description": "Tiny dataset for quick testing",
            "num_chromosomes": 2,
            "base_length": 10000,
            "num_samples": 2,
            "variants_per_haplotype": 10,
            "homozygous_fraction": 0.2,
            "variant_types": ["snv", "insertion", "deletion"],
            "gc_content": 0.42,
            "cores": 2
        },
        "small": {
            "description": "Small dataset for development and testing",
            "num_chromosomes": 3,
            "base_length": 50000,
            "num_samples": 3,
            "variants_per_haplotype": 25,
            "homozygous_fraction": 0.3,
            "variant_types": ["snv", "insertion", "deletion"],
            "gc_content": 0.42,
            "cores": 4
        },
        "medium": {
            "description": "Medium dataset for validation",
            "num_chromosomes": 5,
            "base_length": 100000,
            "num_samples": 5,
            "variants_per_haplotype": 50,
            "homozygous_fraction": 0.25,
            "variant_types": ["snv", "insertion", "deletion"],
            "gc_content": 0.42,
            "cores": 8
        },
        "large": {
            "description": "Large dataset for performance testing",
            "num_chromosomes": 10,
            "base_length": 200000,
            "num_samples": 10,
            "variants_per_haplotype": 100,
            "homozygous_fraction": 0.2,
            "variant_types": ["snv", "insertion", "deletion"],
            "gc_content": 0.42,
            "cores": 16
        }
    }
    
    if preset not in presets:
        raise ValueError(f"Unknown preset: {preset}. Available presets: {list(presets.keys())}")
    
    return presets[preset]


def generate_reference_genome(config: Dict[str, Any], output_dir: Path) -> Path:
    """Generate reference genome."""
    reference_path = output_dir / "reference.fa"
    
    script_path = Path(__file__).parent / "test_data" / "generate_reference_genome.py"
    
    cmd = [
        sys.executable, str(script_path),
        "--output", str(reference_path),
        "--num-chromosomes", str(config["num_chromosomes"]),
        "--chr-length", str(config["base_length"]),
        "--gc-content", str(config["gc_content"]),
        "--add-cpg-islands",
        "--add-telomeres",
        "--seed", str(config.get("seed", 42))
    ]
    
    run_command(cmd, "Generating reference genome")
    
    return reference_path


def generate_variant_samples(config: Dict[str, Any], reference_path: Path, output_dir: Path) -> Path:
    """Generate variant samples."""
    samples_dir = output_dir / "samples"
    
    script_path = Path(__file__).parent / "test_data" / "generate_variants.py"
    
    cmd = [
        sys.executable, str(script_path),
        "--reference", str(reference_path),
        "--output-dir", str(samples_dir),
        "--num-samples", str(config["num_samples"]),
        "--variants-per-haplotype", str(config["variants_per_haplotype"]),
        "--variant-types", ",".join(config["variant_types"]),
        "--homozygous-fraction", str(config["homozygous_fraction"]),
        "--seed", str(config.get("seed", 42) + 1)  # Different seed for variants
    ]
    
    run_command(cmd, "Generating variant samples")
    
    return samples_dir


def build_pangenome(config: Dict[str, Any], reference_path: Path, samples_dir: Path, output_dir: Path) -> Path:
    """Build pangenome graph."""
    pangenome_dir = output_dir / "pangenome"
    
    script_path = Path(__file__).parent / "build_pangenome.py"
    
    # Find all sample directories
    sample_dirs = [str(d) for d in samples_dir.iterdir() if d.is_dir() and d.name.startswith("sample_")]
    
    cmd = [
        sys.executable, str(script_path),
        "--reference", str(reference_path),
        "--sample-dirs"] + sample_dirs + [
        "--output", str(pangenome_dir),
        "--cores", str(config["cores"])
    ]
    
    # Add Docker image options if specified
    if "cactus_docker" in config:
        cmd.extend(["--cactus-docker", config["cactus_docker"]])
    if "vg_docker" in config:
        cmd.extend(["--vg-docker", config["vg_docker"]])
    
    run_command(cmd, "Building pangenome graph")
    
    return pangenome_dir


def call_variants(config: Dict[str, Any], reference_path: Path, pangenome_dir: Path, output_dir: Path) -> Path:
    """Call variants from pangenome."""
    vcf_dir = output_dir / "vcf_output"
    
    script_path = Path(__file__).parent / "call_variants.py"
    
    pangenome_vg = pangenome_dir / "pangenome.vg"
    
    cmd = [
        sys.executable, str(script_path),
        "--pangenome-vg", str(pangenome_vg),
        "--reference", str(reference_path),
        "--output", str(vcf_dir),
        "--cores", str(config["cores"]),
        "--threads", str(config.get("threads", config["cores"]))
    ]
    
    # Add Docker image option if specified
    if "vg_docker" in config:
        cmd.extend(["--vg-docker", config["vg_docker"]])
    
    run_command(cmd, "Calling variants from pangenome")
    
    return vcf_dir


def write_config_summary(config: Dict[str, Any], output_dir: Path) -> None:
    """Write configuration summary to output directory."""
    config_path = output_dir / "generation_config.yaml"
    
    with open(config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    
    logging.info(f"Configuration written to {config_path}")


def print_summary(config: Dict[str, Any], output_dir: Path) -> None:
    """Print summary of generated data."""
    print("\n" + "="*60)
    print("TEST DATA GENERATION COMPLETE")
    print("="*60)
    print(f"Dataset: {config.get('description', 'Custom dataset')}")
    print(f"Output directory: {output_dir}")
    print()
    print("Generated files:")
    print(f"  📄 Reference genome: reference.fa")
    print(f"  📁 Sample variants: samples/ ({config['num_samples']} samples)")
    print(f"  🧬 Pangenome graph: pangenome/pangenome.vg")
    print(f"  📊 VCF files: vcf_output/ ({config['num_samples']} individual + 1 combined)")
    print()
    print("Dataset statistics:")
    print(f"  Chromosomes: {config['num_chromosomes']}")
    print(f"  Base chromosome length: {config['base_length']:,} bp")
    print(f"  Samples: {config['num_samples']}")
    print(f"  Variants per haplotype: {config['variants_per_haplotype']}")
    print(f"  Homozygous fraction: {config['homozygous_fraction']:.1%}")
    print(f"  Variant types: {', '.join(config['variant_types'])}")
    print()
    print("Quick validation commands:")
    print(f"  # Check reference genome")
    print(f"  grep -c '^>' {output_dir}/reference.fa")
    print()
    print(f"  # Check sample count")
    print(f"  ls {output_dir}/samples/")
    print()
    print(f"  # Check VCF files")
    print(f"  ls {output_dir}/vcf_output/*.vcf.gz")
    print()
    print(f"  # View variant statistics")
    print(f"  cat {output_dir}/vcf_output/calling_stats.txt")
    print("="*60)


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Generate complete test datasets for pangenome analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Presets:
  tiny    - 2 chromosomes, 2 samples, 10 variants/haplotype (quick testing)
  small   - 3 chromosomes, 3 samples, 25 variants/haplotype (development)
  medium  - 5 chromosomes, 5 samples, 50 variants/haplotype (validation)
  large   - 10 chromosomes, 10 samples, 100 variants/haplotype (performance)

Examples:
  # Generate small test dataset
  python scripts/generate_test_data.py --preset small

  # Generate medium dataset in custom location
  python scripts/generate_test_data.py --preset medium --output-dir data/my_test

  # Custom dataset with specific parameters
  python scripts/generate_test_data.py --custom \\
    --num-samples 4 --variants-per-haplotype 75 \\
    --homozygous-fraction 0.4 --cores 8

  # Only generate up to pangenome (skip variant calling)
  python scripts/generate_test_data.py --preset small --stop-after pangenome
        """
    )
    
    # Preset or custom mode
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "--preset",
        choices=["tiny", "small", "medium", "large"],
        help="Use predefined configuration preset"
    )
    mode_group.add_argument(
        "--custom",
        action="store_true",
        help="Use custom configuration (specify parameters below)"
    )
    
    # Output options
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/test"),
        help="Output directory for generated data (default: data/test)"
    )
    
    parser.add_argument(
        "--stop-after",
        choices=["reference", "variants", "pangenome", "vcf"],
        help="Stop pipeline after specified stage"
    )
    
    # Custom configuration options
    custom_group = parser.add_argument_group("Custom configuration options")
    custom_group.add_argument("--num-chromosomes", type=int, default=3, help="Number of chromosomes")
    custom_group.add_argument("--base-length", type=int, default=100000, help="Base chromosome length")
    custom_group.add_argument("--num-samples", type=int, default=3, help="Number of samples")
    custom_group.add_argument("--variants-per-haplotype", type=int, default=50, help="Variants per haplotype")
    custom_group.add_argument("--homozygous-fraction", type=float, default=0.3, help="Fraction of homozygous variants")
    custom_group.add_argument("--variant-types", default="snv,insertion,deletion", help="Variant types")
    custom_group.add_argument("--gc-content", type=float, default=0.42, help="GC content")
    custom_group.add_argument("--cores", type=int, default=8, help="Number of CPU cores")
    
    # Docker options
    docker_group = parser.add_argument_group("Docker options")
    docker_group.add_argument("--cactus-docker", default="quay.io/comparative-genomics-toolkit/cactus:v2.8.4")
    docker_group.add_argument("--vg-docker", default="quay.io/vgteam/vg:v1.65.0")
    
    # General options
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without executing")
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    # Get configuration
    if args.preset:
        config = get_preset_config(args.preset)
        logging.info(f"Using preset: {args.preset}")
        logging.info(f"Description: {config['description']}")
    else:
        # Build custom config from arguments
        config = {
            "description": "Custom dataset",
            "num_chromosomes": args.num_chromosomes,
            "base_length": args.base_length,
            "num_samples": args.num_samples,
            "variants_per_haplotype": args.variants_per_haplotype,
            "homozygous_fraction": args.homozygous_fraction,
            "variant_types": [vt.strip() for vt in args.variant_types.split(",")],
            "gc_content": args.gc_content,
            "cores": args.cores
        }
        logging.info("Using custom configuration")
    
    # Add additional options to config
    config.update({
        "seed": args.seed,
        "cactus_docker": args.cactus_docker,
        "vg_docker": args.vg_docker
    })
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Write configuration
    write_config_summary(config, args.output_dir)
    
    if args.dry_run:
        print("DRY RUN - Would execute the following pipeline:")
        print(f"1. Generate reference genome ({config['num_chromosomes']} chromosomes)")
        print(f"2. Generate variant samples ({config['num_samples']} samples)")
        print(f"3. Build pangenome graph")
        print(f"4. Call variants to generate VCF files")
        print(f"\nOutput directory: {args.output_dir}")
        return
    
    try:
        # Pipeline execution
        logging.info("Starting test data generation pipeline")
        
        # Stage 1: Generate reference genome
        logging.info("Stage 1/4: Generating reference genome")
        reference_path = generate_reference_genome(config, args.output_dir)
        
        if args.stop_after == "reference":
            logging.info("Stopping after reference generation as requested")
            return
        
        # Stage 2: Generate variant samples
        logging.info("Stage 2/4: Generating variant samples")
        samples_dir = generate_variant_samples(config, reference_path, args.output_dir)
        
        if args.stop_after == "variants":
            logging.info("Stopping after variant generation as requested")
            return
        
        # Stage 3: Build pangenome
        logging.info("Stage 3/4: Building pangenome graph")
        pangenome_dir = build_pangenome(config, reference_path, samples_dir, args.output_dir)
        
        if args.stop_after == "pangenome":
            logging.info("Stopping after pangenome construction as requested")
            return
        
        # Stage 4: Call variants
        logging.info("Stage 4/4: Calling variants from pangenome")
        vcf_dir = call_variants(config, reference_path, pangenome_dir, args.output_dir)
        
        # Success!
        logging.info("Pipeline completed successfully")
        print_summary(config, args.output_dir)
        
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
