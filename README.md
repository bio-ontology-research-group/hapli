# Hapli: Pangenome Test Data Generation and Variant Calling

## Overview

Hapli is a comprehensive toolkit for generating synthetic genomic test data and performing variant calling using pangenome graphs. The pipeline supports the complete workflow from reference genome generation through variant calling, producing realistic test datasets for pangenome analysis.

## Quick Start

Generate a complete test dataset with a single command:

`python scripts/generate_test_data.py --preset <size> --output-dir data/large_test`
where `size` can be `small`, `medium`, or `large`.

