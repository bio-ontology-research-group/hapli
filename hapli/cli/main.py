import typer
import logging
from pathlib import Path
from typing import Optional
from typing_extensions import Annotated
from hapli.workflow.pipeline import HapliPipeline

app = typer.Typer(help="Hapli: Genotype-centric variant analysis.", no_args_is_help=True)

def setup_logging(verbose: bool):
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

@app.command()
def analyze(
    gene: str = typer.Option(..., help="Target gene name/ID"),
    sample: str = typer.Option(..., help="Sample name in VCF"),
    vcf: Path = typer.Option(..., help="Path to phased VCF"),
    reference: Path = typer.Option(..., help="Path to reference FASTA"),
    gff: Path = typer.Option(..., help="Path to GFF3 annotation"),
    output_dir: Path = typer.Option(..., help="Directory for results"),
    with_esm: bool = typer.Option(
        False, "--with-esm",
        help="Compute the ESM2 epistasis residual (requires the `ml` extra).",
    ),
    esm_model: str = typer.Option(
        "esm2_t6_8M_UR50D",
        help="fair-esm checkpoint name (ignored unless --with-esm).",
    ),
    gnomad_constraint: Optional[Path] = typer.Option(
        None, "--gnomad-constraint",
        help="gnomAD constraint TSV (canonical per-gene pLI / MisZ / o/e).",
    ),
    clingen_dosage: Optional[Path] = typer.Option(
        None, "--clingen-dosage",
        help="ClinGen dosage sensitivity TSV.",
    ),
    alphamissense_table: Optional[Path] = typer.Option(
        None, "--alphamissense-table",
        help="bgzipped+tabix-indexed AlphaMissense TSV (per-variant lookup).",
    ),
    clinvar_vcf: Optional[Path] = typer.Option(
        None, "--clinvar-vcf",
        help="ClinVar VCF (bgzipped+tabix-indexed) — annotates consequences with clnsig.",
    ),
    verbose: bool = False,
):
    """
    Analyze a specific gene for a sample, generating haplotypes and checking feature preservation.
    """
    setup_logging(verbose)

    pipeline = HapliPipeline(
        reference, vcf, gff, output_dir,
        with_esm=with_esm, esm_model=esm_model,
        gnomad_constraint=gnomad_constraint,
        clingen_dosage=clingen_dosage,
        alphamissense_table=alphamissense_table,
        clinvar_vcf=clinvar_vcf,
    )
    try:
        pipeline.run_gene_analysis(gene, sample)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1) from e

@app.command()
def assess(
    gene: str = typer.Option(..., help="Target gene name/ID"),
    sample: str = typer.Option(..., help="Sample name used in output filenames"),
    hap1: Path = typer.Option(..., help="Pre-assembled haplotype 1 FASTA (indexed or will be indexed)"),
    hap2: Path = typer.Option(..., help="Pre-assembled haplotype 2 FASTA"),
    reference: Path = typer.Option(..., help="Reference FASTA matching the GFF"),
    gff: Path = typer.Option(..., help="Reference GFF3 annotation to lift onto each haplotype"),
    output_dir: Path = typer.Option(..., help="Directory for results"),
    with_esm: bool = typer.Option(
        False, "--with-esm",
        help="Compute the ESM2 epistasis residual (requires the `ml` extra).",
    ),
    esm_model: str = typer.Option(
        "esm2_t6_8M_UR50D",
        help="fair-esm checkpoint name (ignored unless --with-esm).",
    ),
    gnomad_constraint: Optional[Path] = typer.Option(
        None, "--gnomad-constraint",
        help="gnomAD constraint TSV (canonical per-gene pLI / MisZ / o/e).",
    ),
    clingen_dosage: Optional[Path] = typer.Option(
        None, "--clingen-dosage",
        help="ClinGen dosage sensitivity TSV.",
    ),
    alphamissense_table: Optional[Path] = typer.Option(
        None, "--alphamissense-table",
        help="bgzipped+tabix-indexed AlphaMissense TSV (per-variant lookup).",
    ),
    clinvar_vcf: Optional[Path] = typer.Option(
        None, "--clinvar-vcf",
        help="ClinVar VCF (bgzipped+tabix-indexed) — annotates consequences with clnsig.",
    ),
    verbose: bool = False,
):
    """
    Mode B: assess a gene on pre-assembled haplotype FASTAs (HPRC, hifiasm, Verkko).

    Runs Liftoff to transfer the reference GFF3 onto each haplotype, extracts
    and diffs the per-haplotype protein against the reference, optionally
    computes the ESM2 epistasis residual, and emits a schema-v2 analysis JSON
    plus the usual per-haplotype protein / annotation artifacts.

    No VCF is consumed. Variant-equivalent information is recovered from the
    protein diff and the lifted GFF coordinates.
    """
    setup_logging(verbose)
    pipeline = HapliPipeline(
        reference, None, gff, output_dir,
        with_esm=with_esm, esm_model=esm_model,
        gnomad_constraint=gnomad_constraint,
        clingen_dosage=clingen_dosage,
        alphamissense_table=alphamissense_table,
        clinvar_vcf=clinvar_vcf,
    )
    try:
        pipeline.run_gene_assess(gene, sample, hap1, hap2)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1) from e


@app.command()
def explore(
    result_json: Annotated[Path, typer.Argument(help="JSON alignment result file to explore")]
):
    """Explore alignment results interactively."""
    from hapli.cli.tui import HapliExplorer
    
    if not result_json.exists():
        typer.echo(f"Error: File {result_json} does not exist.", err=True)
        raise typer.Exit(code=1)
        
    tui = HapliExplorer(result_json)
    tui.run()

@app.command()
def interpret(
    result_json: Annotated[Path, typer.Argument(help="JSON alignment result file to interpret")],
    api_key: str = typer.Option(None, envvar="OPENROUTER_API_KEY", help="OpenRouter API Key"),
    model: str = "openai/gpt-3.5-turbo"
):
    """Generate LLM-based functional interpretation."""
    from hapli.interpretation.llm import LLMInterpreter
    import json
    
    if not result_json.exists():
        typer.echo(f"Error: File {result_json} does not exist.", err=True)
        raise typer.Exit(code=1)
        
    with open(result_json) as f:
        data = json.load(f)
        
    interpreter = LLMInterpreter(api_key, model)
    report = interpreter.interpret(data)
    
    typer.echo("--- Interpretation Report ---")
    typer.echo(report)
    typer.echo("-----------------------------")

@app.command()
def aggregate(
    glob_pattern: str = typer.Argument(..., help="Glob over hapli analysis JSONs, e.g. 'results/*_*_analysis.json'"),
    out_per_sample: Path = typer.Option(..., "--per-sample", help="Per-(sample, gene) TSV output"),
    out_per_gene: Path = typer.Option(..., "--per-gene", help="Per-gene population-aggregate TSV"),
    lof_threshold: float = typer.Option(
        0.5, help="Per-haplotype score at or below which the haplotype counts as LoF (default: 0.5)"
    ),
):
    """
    Aggregate many `analyze` / `assess` outputs into population-level TSVs.

    Two outputs:
      - per-sample TSV: one row per (sample, gene) with key signals.
      - per-gene TSV: one row per gene with population aggregates
        (LoF allele frequency, compound-het-LoF count, epistasis-flag rate).

    The natural HPRC-scale consumer of `analyze` outputs.
    """
    from hapli.workflow.aggregate import aggregate_glob, write_per_sample_tsv, write_per_gene_tsv

    rows, aggregates = aggregate_glob(glob_pattern, lof_threshold=lof_threshold)
    write_per_sample_tsv(rows, out_per_sample)
    write_per_gene_tsv(aggregates, out_per_gene)
    typer.echo(
        f"Aggregated {len(rows)} per-(sample,gene) records over {len(aggregates)} genes.\n"
        f"  per-sample → {out_per_sample}\n"
        f"  per-gene   → {out_per_gene}"
    )


@app.command()
def version():
    """Show version."""
    try:
        from importlib.metadata import version as _pkg_version
        typer.echo(f"hapli {_pkg_version('hapli')}")
    except Exception:
        typer.echo("hapli (development, uninstalled)")

if __name__ == "__main__":
    app()
