"""Script to fetch the sequences given by an imodulon."""

import sys
import pandas as pd
import gzip
import typer

from Bio.SeqIO.FastaIO import SimpleFastaParser

COMPL = {"A": "T", "T": "A", "G": "C", "C": "G"}


def fetch_upstream(genome: str, position: int, upstream_len: int = 200) -> str:
    slice = genome[max(position - upstream_len, 0) : position]
    return slice


def read_genome(path: str):
    """Read a gzip genome (assuming only/first in fasta)."""
    opener = gzip.open if path.endswith("gz") else open
    flags = "rt" if path.endswith("gz") else "r"
    with opener(sys.argv[2], flags) as f:
        return next(seq for _, seq in SimpleFastaParser(f))


def run(input_csv: str, genome_fasta_gz: str, out_csv: str):
    """Transform the input form process_imodulons.py to sequences."""
    df = pd.read_csv(input_csv, index_col=0)
    genome = read_genome(genome_fasta_gz)
    df["seq"] = df.apply(lambda x: fetch_upstream(genome, x.position, 200), axis=1)
    df["geneid"] = df.index
    df.to_csv(out_csv, index=False)


if __name__ == "__main__":
    # the input is the output of scripts/process_imodulons.py
    typer.run(run)
