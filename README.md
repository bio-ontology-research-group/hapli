# Hapli: Pangenome Test Data Generation and Variant Calling

## Overview

Hapli is a comprehensive toolkit for generating synthetic genomic test data and performing variant calling using pangenome graphs. The pipeline supports the complete workflow from reference genome generation through variant calling, producing realistic test datasets for pangenome analysis.

## Quick Start

Generate a complete test dataset with a single command:

`python scripts/generate_test_data.py --preset <size> --output-dir data/large_test`
where `size` can be `small`, `medium`, or `large`.

## Core Components

### GAM File Parser (`hapli/gam_parser.py`)

The `GAMParser` class provides functionality to parse GAM (Graph Alignment/Map) files and organize alignment data by sample and haplotype path names. This is useful for analyzing pangenome alignments and understanding how features align across different samples and haplotypes.

#### Features:
- Reads GAM alignments using vg's JSON output format
- Extracts sample and haplotype information from path names
- Groups alignments by sample/haplotype combinations
- Builds hierarchical feature structures from alignment metadata
- Returns structured data as nested dictionaries: `{sample: {haplotype: {feature_type: [alignment_data]}}}`
- Provides alignment statistics and quality metrics

#### Usage:

