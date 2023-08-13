from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable, List

from Bio.SeqRecord import SeqRecord


def clip_lenghts(
    genes: List[SeqRecord], from_percentile: float, to_percentile: float
) -> List[SeqRecord]:
    sorted_by_len = sorted(genes, key=lambda seq_object: len(seq_object.seq))
    from_idx = int(len(genes) * from_percentile)
    to_idx = int(len(genes) * to_percentile)
    return sorted_by_len[from_idx:to_idx]


def filter_non_acgt(genes: List[SeqRecord]) -> List[SeqRecord]:
    return [gene for gene in genes if set(str(gene.seq)).issubset("ACGT")]


def preprocess(genes: List[SeqRecord]) -> List[SeqRecord]:
    genes = filter_non_acgt(genes)
    genes = clip_lenghts(genes, from_percentile=0.1, to_percentile=0.9)
    return genes
