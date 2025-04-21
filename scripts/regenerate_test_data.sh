#!/bin/bash
set -e

echo "Regenerating test data..."
mkdir -p data
python scripts/generate_example_data.py --ref-length 1000 --num-features 5 --output-dir data --basename example
echo "Test data regenerated successfully."
