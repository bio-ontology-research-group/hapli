# Integration Test Dataset

This directory contains a comprehensive test dataset for the Haplotype Annotation Tool. It's designed to test various aspects of feature annotation on variation graphs.

## Files

- `reference.fasta`: Reference genome sequence (~30kb) that corresponds to the GFF3 annotations
- `annotations.gff3`: Gene annotations with diverse feature types and hierarchical relationships
- `variation_graph.gfa`: GFA2 file representing a variation graph with multiple haplotypes
- `variants.vcf`: VCF file containing variants corresponding to differences in the graph paths

## Dataset Structure

### Samples and Haplotypes

The dataset contains two samples, each with two haplotypes:
- Sample1: haplotype1, haplotype2
- Sample2: haplotype1, haplotype2

Path naming follows the pattern: `sample{sample_id}_hap{haplotype_id}`

### Feature Types

The GFF3 file includes a variety of feature types with hierarchical relationships:
- Genes (10 genes total)
- mRNAs (1-3 transcript isoforms per gene)
- Exons (2-8 per transcript)
- CDS (coding sequences)
- 5' and 3' UTRs

### Test Cases

This dataset is designed to test the following scenarios:

1. **Perfect Feature Alignment**
   - Gene1 and Gene2 align perfectly to all haplotypes
   - Tests basic alignment functionality

2. **SNP Variants**
   - Gene3 contains SNPs in Sample1-Hap2 and Sample2-Hap1
   - These correspond to entries in the VCF file
   - Tests variant detection for single nucleotide changes

3. **Insertion/Deletion Variants**
   - Gene4 contains small indels in Sample2 haplotypes
   - Gene5 contains a large deletion in Sample1-Hap1
   - Tests insertion/deletion detection and impact classification

4. **Complex Structural Variants**
   - Gene6 is duplicated in Sample2-Hap2
   - Gene7 is inverted in Sample1-Hap2
   - Tests handling of complex rearrangements

5. **Partially Present Features**
   - Gene8 is truncated in Sample2 haplotypes
   - Tests proper impact classification of truncated features

6. **Absent Features**
   - Gene9 is completely absent from Sample1-Hap1
   - Tests detection of missing features

7. **Alternative Splicing**
   - Gene10 has different exon usage across haplotypes
   - Tests alignment of features with alternative structures

8. **Feature Hierarchy Reconciliation**
   - Several genes have child features (exons, CDS) with boundary conflicts
   - Tests the feature reconciliation algorithms

## Usage in Testing

This dataset is ideal for integration testing of the entire annotation pipeline:

```bash
python -m src.main --gfa-file data/integration_test/variation_graph.gfa \
                  --gff3-file data/integration_test/annotations.gff3 \
                  --reference-fasta data/integration_test/reference.fasta
```

Validate the test data integrity with:

```bash
python -m tests.validate_test_data
```

## Regenerating Test Data

If you need to regenerate or modify the test data, you can use the included generator script:

```bash
python data/integration_test/generate_test_data.py
```

This script uses appropriate bioinformatics libraries (Biopython, GFApy) to ensure all files are valid and properly 
cross-referenced. It creates:

1. A reference FASTA sequence
2. GFF3 file with gene annotations
3. GFA2 variation graph with variant paths
4. VCF file with corresponding variants

The generator ensures that variants introduced in the graph are properly aligned with features and coordinates match between all files.
