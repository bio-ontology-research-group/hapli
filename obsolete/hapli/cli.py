import typer
import logging
from pathlib import Path
from typing import Optional
from typing_extensions import Annotated
import sys

from hapli.core.io import GFFProcessor, SequenceExtractor
from hapli.alignment.hierarchical import HierarchicalAligner
from hapli.variation.haplotype import HaplotypeGenerator
from hapli.tui import HapliExplorer

app = typer.Typer(
    name="hapli",
    help="Genotype-centric pangenome variant analysis toolkit.",
    add_completion=False,
    no_args_is_help=True
)

def setup_logging(verbose: bool):
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

@app.command()
def generate(
    vcf: Annotated[Path, typer.Option(..., help="Phased VCF file")],
    reference: Annotated[Path, typer.Option(..., help="Reference FASTA")],
    region: Annotated[str, typer.Option(..., help="Region (chr:start-end)")],
    sample: Annotated[str, typer.Option(..., help="Sample name")],
    output: Annotated[Path, typer.Option(..., help="Output FASTA file")],
    verbose: bool = False
):
    """Generate haplotype sequences from VCF."""
    setup_logging(verbose)
    
    try:
        chrom, pos_str = region.split(':')
        start, end = map(int, pos_str.split('-'))
    except ValueError:
        typer.echo("Error: Region must be in format chr:start-end", err=True)
        raise typer.Exit(code=1)

    generator = HaplotypeGenerator(reference, vcf)
    if output.exists():
        output.unlink()
        
    generator.generate_haplotypes_for_region(chrom, start, end, sample, output)
    typer.echo(f"Haplotypes generated at {output}")

@app.command()
def align(
    haplotypes: Annotated[Path, typer.Option(..., help="Haplotype FASTA")],
    gff: Annotated[Path, typer.Option(..., help="GFF3 file")],
    reference: Annotated[Path, typer.Option(..., help="Reference FASTA")],
    gene: Annotated[str, typer.Option(..., help="Gene Name/ID")],
    output: Annotated[Path, typer.Option(..., help="Output JSON")],
    threads: int = 4,
    verbose: bool = False
):
    """Align a gene hierarchically to haplotypes."""
    setup_logging(verbose)
    
    import json
    
    gff_proc = GFFProcessor(gff, target_gene=gene)
    seq_ext = SequenceExtractor(reference)
    aligner = HierarchicalAligner(gff_proc, seq_ext, threads=threads)
    
    # Find the gene object
    gene_feature = None
    for feat in gff_proc.features_by_id.values():
         if feat.featuretype == 'gene':
             # Check ID, Name, gene_name
             attrs = feat.attributes
             names = [feat.id] + attrs.get('Name', []) + attrs.get('gene_name', [])
             if gene in names:
                 gene_feature = feat
                 break
    
    # Fallback: if only one gene loaded, use it
    if not gene_feature:
        genes = [f for f in gff_proc.features_by_id.values() if f.featuretype == 'gene']
        if len(genes) == 1:
            gene_feature = genes[0]
            logging.info(f"Exact match not found, using loaded gene: {gene_feature.id}")

    if not gene_feature:
        typer.echo(f"Error: Gene {gene} not found in loaded features.", err=True)
        raise typer.Exit(code=1)
        
    results = aligner.align_gene(gene_feature, haplotypes)
    
    with open(output, 'w') as f:
        json.dump({k: v.to_dict() for k, v in results.items()}, f, indent=2)
    typer.echo(f"Alignment results saved to {output}")

@app.command()
def explore(
    result_json: Annotated[Path, typer.Argument(help="JSON alignment result file to explore")]
):
    """Explore alignment results interactively."""
    if not result_json.exists():
        typer.echo(f"Error: File {result_json} does not exist.", err=True)
        raise typer.Exit(code=1)
        
    tui = HapliExplorer(result_json)
    tui.run()

def main():
    app()

if __name__ == "__main__":
    main()