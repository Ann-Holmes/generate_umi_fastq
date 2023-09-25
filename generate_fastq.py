#!/usr/bin/env python3
import numpy as np
import click
from Bio


@click.group()
def cli():
    pass


@click.command(help="Generate DNA sequence by nt size")
@click.argument("nt", type=int)
@click.option("--debug", is_flag=True, help="Debug mode")
def generate_DNA(nt: int, debug: bool = False) -> str:
    """Generate DNA sequence by nt size"""
    dna = "".join(np.random.choice(list("ATCG"), nt, replace=True))
    if debug:
        click.echo(dna)
    return dna


@click.command(help="Generate reference genome by number chromosomes and each chromosome size")
@click.argument("n_chrs", type=int)
@click.argument("chrs_size", type=list)
@click.option("--debug", is_flag=True, help="Debug mode")
def generate_reference_genome(n_chrs: int, chrs_size: list, debug: bool = False) -> str:
    """Generate reference genome by number chromosomes and each chromosome size"""
    chrs = []
    for i in range(n_chrs):
        chrs.append(generate_DNA(chrs_size[i], debug))
    if debug:
        click.echo(chrs)
    return chrs



cli.add_command(generate_DNA)
cli.add_command(generate_reference_genome)
if __name__ == '__main__':
    # cli.add_command(generate_DNA)
    cli()
