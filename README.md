# Hapli: Pangenome Test Data Generation and Variant Calling

## Overview

Hapli is a comprehensive toolkit for generating synthetic genomic test data and performing variant calling using pangenome graphs. The pipeline supports the complete workflow from reference genome generation through variant calling, producing realistic test datasets for pangenome analysis.

## Quick Start

Generate a complete test dataset with a single command:

`python scripts/generate_test_data.py --preset <size> --output-dir data/large_test`
where `size` can be `small`, `medium`, or `large`.

## Core Analysis Workflow

Hapli provides a four-step analysis workflow for detecting variant impacts on genomic features using pangenome graphs:

### 1. GFF Alignment (`hapli/gff_alignment.py`)

Aligns genomic features from a GFF3 file to paths in a pangenome graph.

**Inputs:**
- GFF3 file with genomic features (genes, exons, etc.)
- Reference genome FASTA file
- Pangenome graph file (GFA format)

**Outputs:**
- GAM file containing feature alignments to graph paths

**Usage:**
