#!/usr/bin/env python3
from typing import List
import numpy as np
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from generate_kmers_dna import generate_kmers, generate_DNA_cy
# Before import generate_kmers and generate_DNA_cy， 
# please run "python setup.py build_ext --inplace" for build cython code


@click.group()
def cli():
    pass


def generate_DNA(nt: int) -> str:
    """Generate DNA sequence by nt size"""
    dna = "".join(np.random.choice(list("ATCG"), nt, replace=True))
    return dna


@click.command(help="Generate reference genome by number chromosomes and each chromosome size")
@click.argument("fasta", type=str)
@click.option("--n_chrs", type=int, default=5, help="Number of chromosomes")
@click.option("--chrs_size", type=click.INT, multiple=True, help="Chromosomes size")
def generate_reference_genome(fasta: str, n_chrs: int, chrs_size: tuple) -> None:
    """Generate reference genome by number chromosomes and each chromosome size"""
    if not chrs_size:
        chrs_size = np.random.randint(10000, 100000, size=n_chrs).tolist()

    chrs = []
    for i, chr_size in zip(range(n_chrs), chrs_size):
        chromosome = SeqRecord(
            Seq(generate_DNA_cy(chr_size)),
            id=f"chr{i + 1}",
            description=""
        )
        chrs.append(chromosome)

    with open(fasta, "w") as output_handle:
        SeqIO.write(chrs, output_handle, "fasta")


def generate_fastq(n_reads: int, read_size: int = 144, n_umis: int = 1000,
                   umi_length: int = 6) -> List[SeqRecord]:
    """Generate fastq file content by number reads and read size"""
    umis = np.array([generate_DNA_cy(umi_length) for i in range(n_umis)])
    fastq = []
    umis = umis[np.random.randint(0, n_umis, size=(n_reads))]
    reads = np.array([generate_DNA_cy(read_size) for i in range(n_reads)])

    reads_seq = np.char.add(umis, reads)
    for i, read_seq in enumerate(reads_seq):
        read = SeqRecord(
            Seq(read_seq),
            id=f"FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{str(i)}",
            description="",
            letter_annotations={"phred_quality": [40] * (read_size + umi_length)}
        )
        fastq.append(read)

    return fastq


@click.command(help="Generate single-end fastq file content by number reads and read size")
@click.argument("fastq", type=str)
@click.option("--n_reads", type=int, default=10000, help="Number of reads")
@click.option("--read_size", type=int, default=144, help="Read size")
@click.option("--n_umis", type=int, default=1000, help="Number of UMIs")
@click.option("--umi_length", type=int, default=6, help="UMI length, set to zero will not add UMI")
def generate_fastq_SE(fastq: str, n_reads: int, read_size: int, n_umis: int,
                      umi_length: int) -> None:
    """Generate fastq file content by number reads and read size"""
    output_fastq = fastq
    fastq = generate_fastq(n_reads, read_size, n_umis, umi_length)

    for read in fastq:
        read.description = "1:N:0"

    with open(output_fastq, "w") as output_handle:
        SeqIO.write(fastq, output_handle, "fastq")


@click.command(help="Generate paired-end fastq file content by number reads and read size")
@click.argument("fastq_bn", type=str)
@click.option("--n_reads", type=int, default=10000, help="Number of reads")
@click.option("--read_size", type=int, default=144, help="Read size")
@click.option("--n_umis", type=int, default=1000, help="Number of UMIs")
@click.option("--umi_length", type=int, default=6, help="UMI length, set to zero will not add UMI")
def generate_fastq_PE(fastq_bn: str, n_reads: int, read_size: int, n_umis: int,
                      umi_length: int) -> None:
    """Generate fastq file content by number reads and read size"""
    output_fastq_R1 = fastq_bn + "_R1.fq"
    fastq_R1 = generate_fastq(n_reads, read_size, n_umis, umi_length)

    # for read in fastq_R1:
    #     read.description = "#0/1"

    with open(output_fastq_R1, "w") as output_handle:
        SeqIO.write(fastq_R1, output_handle, "fastq")

    output_fastq_R2 = fastq_bn + "_R2.fq"
    fastq_R2 = generate_fastq(n_reads, read_size, n_umis, umi_length)

    # for read in fastq_R2:
    #     read.description = "#0/2"

    with open(output_fastq_R2, "w") as output_handle:
        SeqIO.write(fastq_R2, output_handle, "fastq")


def generate_mapped_fastq(fasta: str, n_reads: int, read_size: int = 146, n_umis: int = 1000,
                          umi_length: int = 6) -> List[SeqRecord]:
    """Generate fastq file by reference genome, number reads and read size"""
    kmers = []
    for record in SeqIO.parse(fasta, "fasta"):
        record = str(record.seq)
        kmers += generate_kmers(record, read_size)

    umis = np.array([generate_DNA_cy(umi_length) for i in range(n_umis)])
    umis = umis[np.random.randint(0, n_umis, n_reads)]
    kmers = np.array(kmers)
    reads = kmers[np.random.randint(0, len(kmers), n_reads)]
    reads_seq = np.char.add(umis, reads)

    fastq = []
    for i, read_seq in enumerate(reads_seq):
        read = SeqRecord(
            Seq(read_seq),
            id=f"FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:{str(i)}",
            description="",
            letter_annotations={"phred_quality": [40] * (read_size + umi_length)}
        )
        fastq.append(read)

    return fastq


@click.command(help="Generate single-end fastq file content by reference genome, number reads and "
               "read size")
@click.argument("fasta", type=str)
@click.argument("fastq", type=str)
@click.option("--n_reads", type=int, default=10000, help="Number of reads")
@click.option("--read_size", type=int, default=144, help="Read size")
@click.option("--n_umis", type=int, default=1000, help="Number of UMIs")
@click.option("--umi_length", type=int, default=6, help="UMI length, set to zero will not add UMI")
def generate_mapped_fastq_SE(fasta: str, fastq: str, n_reads: int, read_size: int, n_umis: int,
                             umi_length: int) -> None:
    """Generate fastq file by reference genome, number reads and read size"""
    output_fastq = fastq
    fastq = generate_mapped_fastq(fasta, n_reads, read_size, n_umis, umi_length)

    for read in fastq:
        read.description = "1:N:0"

    with open(output_fastq, "w") as output_handle:
        SeqIO.write(fastq, output_handle, "fastq")


@click.command(help="Generate paired-end fastq file content by reference genome, number reads and "
               "read size")
@click.argument("fasta", type=str)
@click.argument("fastq_bn", type=str)
@click.option("--n_reads", type=int, default=10000, help="Number of reads")
@click.option("--read_size", type=int, default=144, help="Read size")
@click.option("--n_umis", type=int, default=1000, help="Number of UMIs")
@click.option("--umi_length", type=int, default=6, help="UMI length, set to zero will not add UMI")
def generate_mapped_fastq_PE(fasta: str, fastq_bn: str, n_reads: int, read_size: int, n_umis: int,
                             umi_length: int) -> None:
    """Generate fastq file by reference genome, number reads and read size"""
    output_fastq_R1 = fastq_bn + "_R1.fq"
    fastq_R1 = generate_mapped_fastq(fasta, n_reads, read_size, n_umis, umi_length)

    # for read in fastq_R1:
    #     read.description = "#0/1"

    with open(output_fastq_R1, "w") as output_handle:
        SeqIO.write(fastq_R1, output_handle, "fastq")

    output_fastq_R2 = fastq_bn + "_R2.fq"
    fastq_R2 = generate_mapped_fastq(fasta, n_reads, read_size, n_umis, umi_length)

    # for read in fastq_R2:
    #     read.description = "#0/2"

    with open(output_fastq_R2, "w") as output_handle:
        SeqIO.write(fastq_R2, output_handle, "fastq")


cli.add_command(generate_reference_genome, "generate_reference_genome")
cli.add_command(generate_fastq_SE, "generate_fastq_SE")
cli.add_command(generate_fastq_PE, "generate_fastq_PE")
cli.add_command(generate_mapped_fastq_SE, "generate_mapped_fastq_SE")
cli.add_command(generate_mapped_fastq_PE, "generate_mapped_fastq_PE")
if __name__ == '__main__':
    cli()
