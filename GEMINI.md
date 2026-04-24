# Hapli Project Context

## Overview
**Hapli** is a genotype-centric variant analysis toolkit designed to move beyond isolated variant analysis. It constructs complete haplotype sequences from phased VCFs and reference genomes, then hierarchically aligns genomic features (Genes -> Transcripts -> Exons) to these haplotypes to assess functional preservation.

## Architecture
The project is a Python-based CLI and library with the following core components:

*   **CLI (`hapli/cli.py`)**: Built with `typer`. Entry point via `main.py`.
*   **Core (`hapli/core/`)**: IO handling (`GFFProcessor`, `SequenceExtractor`) and data models.
*   **Variation (`hapli/variation/`)**: Haplotype generation logic (`HaplotypeGenerator`).
*   **Alignment (`hapli/alignment/`)**: Hierarchical alignment logic (`HierarchicalAligner`) using `minimap2`.
*   **TUI (`hapli/tui.py`)**: Interactive results explorer built with `textual`.

## Build & Dependencies
Managed by `uv` (modern Python package manager).

*   **Python Version**: >=3.10.16
*   **Key Python Dependencies**:
    *   `biopython`: Sequence manipulation.
    *   `pysam`: VCF/BAM handling.
    *   `gffutils`: GFF3 parsing.
    *   `rich`, `textual`: CLI output and TUI.
    *   `snakemake`: Workflow automation.
    *   `typer`: CLI construction.
*   **External Dependencies**:
    *   `minimap2`: Required in `$PATH` for alignment.
    *   `bcftools`/`tabix`: Recommended for VCF operations.

### Installation
```bash
uv sync  # Install dependencies
# OR
pip install .
```

## Key Commands
Run commands via `uv run main.py <command>` or `hapli <command>` if installed.

1.  **Generate Haplotypes**:
    ```bash
    uv run main.py generate --vcf <vcf> --reference <ref> --region <chr:start-end> --sample <name> --output <out.fa>
    ```
2.  **Hierarchical Alignment**:
    ```bash
    uv run main.py align --haplotypes <hap.fa> --gff <features.gff> --reference <ref> --gene <gene_name> --output <res.json>
    ```
3.  **Explore Results (TUI)**:
    ```bash
    uv run main.py explore <results.json>
    ```

## Automated Workflows
Located in `workflows/`.
*   **Engine**: `snakemake`
*   **Execution**:
    ```bash
    uv run snakemake -s workflows/Snakefile --configfile workflows/config.yaml --cores <N>
    ```
    *Note: The Snakefile command definitions may need verification against the current `main.py` CLI interface.*

## Development
*   **Testing**: `uv run pytest`
*   **Conventions**:
    *   Adhere to Python type hinting.
    *   Use `pathlib` for file paths.
    *   `rich` is used for console output.
