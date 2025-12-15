import argparse
import logging
import sys
import json
from pathlib import Path

from hapli.core.io import GFFProcessor, SequenceExtractor
from hapli.alignment.hierarchical import HierarchicalAligner
from hapli.variation.haplotype import HaplotypeGenerator

def setup_logging(verbose: bool):
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

def main():
    parser = argparse.ArgumentParser(description="Hapli: Genotype-centric analysis tool.")
    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Command: generate-haplotypes
    gen_parser = subparsers.add_parser("generate-haplotypes", help="Generate haplotype sequences from VCF")
    gen_parser.add_argument("--vcf", required=True, type=Path, help="Phased VCF file")
    gen_parser.add_argument("--reference", required=True, type=Path, help="Reference FASTA")
    gen_parser.add_argument("--region", required=True, help="Region (chr:start-end)")
    gen_parser.add_argument("--sample", required=True, help="Sample name")
    gen_parser.add_argument("--output", required=True, type=Path, help="Output FASTA file")

    # Command: align-gene
    align_parser = subparsers.add_parser("align-gene", help="Align a gene hierarchically to haplotypes")
    align_parser.add_argument("--haplotypes", required=True, type=Path, help="Haplotype FASTA")
    align_parser.add_argument("--gff", required=True, type=Path, help="GFF3 file")
    align_parser.add_argument("--reference", required=True, type=Path, help="Reference FASTA (for feature extraction)")
    align_parser.add_argument("--gene", required=True, help="Gene Name/ID")
    align_parser.add_argument("--output", required=True, type=Path, help="Output JSON")
    align_parser.add_argument("--threads", type=int, default=4)

    args = parser.parse_args()
    setup_logging(True)

    if args.command == "generate-haplotypes":
        chrom, pos_str = args.region.split(':')
        start, end = map(int, pos_str.split('-'))
        
        generator = HaplotypeGenerator(args.reference, args.vcf)
        # Clear output file first
        if args.output.exists():
            args.output.unlink()
            
        generator.generate_haplotypes_for_region(chrom, start, end, args.sample, args.output)
        logging.info(f"Haplotypes generated at {args.output}")

    elif args.command == "align-gene":
        gff_proc = GFFProcessor(args.gff, target_gene=args.gene)
        seq_ext = SequenceExtractor(args.reference)
        aligner = HierarchicalAligner(gff_proc, seq_ext, threads=args.threads)
        
        gene_feature = None
        for feat in gff_proc.features_by_id.values():
             if feat.featuretype == 'gene': # This might be weak if multiple genes loaded
                 # But GFFProcessor with target_gene only loads one gene hierarchy essentially
                 if feat.id == args.gene or feat.attributes.get('Name', [''])[0] == args.gene:
                     gene_feature = feat
                     break
        
        # Fallback if ID didn't match exactly but we loaded it
        if not gene_feature:
             # Find ANY gene
             for feat in gff_proc.features_by_id.values():
                 if feat.featuretype == 'gene':
                     gene_feature = feat
                     break
        
        if not gene_feature:
            logging.error(f"Gene {args.gene} not found in loaded features.")
            sys.exit(1)
            
        results = aligner.align_gene(gene_feature, args.haplotypes)
        
        with open(args.output, 'w') as f:
            json.dump({k: v.to_dict() for k, v in results.items()}, f, indent=2)
        logging.info(f"Alignment results saved to {args.output}")

    else:
        parser.print_help()

if __name__ == "__main__":
    main()
