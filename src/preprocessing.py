
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path
from typing import List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _filter_non_acgt(genes: List[SeqRecord]) -> List[SeqRecord]:
    return [gene for gene in genes if set(str(gene.seq)).issubset("ACGT")]


def _bin_by_length(genes: List[SeqRecord], bin_width: int) -> List[List[SeqRecord]]:
    bins = defaultdict(list)
    for gene in genes:
        gene_length = len(gene.seq)
        bin_ = gene_length - (gene_length % bin_width)
        bins[bin_].append(gene)
    return list(bins.values())


def preprocess(
    genes: List[SeqRecord],
    bin_width: int,
    min_bin_size: int,
) -> List[List[SeqRecord]]:
    genes = _filter_non_acgt(genes)
    bins = _bin_by_length(genes, bin_width)
    bins = [bin_ for bin_ in bins if len(bin_) > min_bin_size]
    return bins


if __name__ == "__main__":
    arg_parser = ArgumentParser()
    arg_parser.add_argument(
        "--in_", type=str, help="Genes in FASTA format to preprocess"
    )
    arg_parser.add_argument("--out", type=str, help="Directory to save to")
    arg_parser.add_argument("--bin-width", type=int, help="Size of bins")
    arg_parser.add_argument(
        "--min-bin-size", type=int, help="Minimum number of sequences per bin"
    )
    args = arg_parser.parse_args()
    genes = SeqIO.parse(args.in_, format="fasta")
    bins = preprocess(list(genes), args.bin_width, args.min_bin_size)
    for bin_ in bins:
        SeqIO.write(
            bin_,
            Path(args.out) / f"{Path(args.in_).stem}_{len(bin_[0])}.fasta",
            format="fasta",
        )
