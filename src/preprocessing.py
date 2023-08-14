from collections import defaultdict
from typing import List

from Bio.SeqRecord import SeqRecord


def _filter_non_acgt(genes: List[SeqRecord]) -> List[SeqRecord]:
    return [gene for gene in genes if set(str(gene.seq)).issubset("ACGT")]


def _clip_lenghts(
    genes: List[SeqRecord], from_percentile: float, to_percentile: float
) -> List[SeqRecord]:
    sorted_by_len = sorted(genes, key=lambda gene: len(gene.seq))
    from_idx = int(len(genes) * from_percentile)
    to_idx = int(len(genes) * to_percentile)
    return sorted_by_len[from_idx:to_idx]


def _bin_by_length(
    genes: List[SeqRecord], bin_width: int
) -> List[List[SeqRecord]]:
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
    from_percentile: float = 0.1,
    to_percentile: float = 0.9,
) -> List[List[SeqRecord]]:
    genes = _filter_non_acgt(genes)
    genes = _clip_lenghts(genes, from_percentile, to_percentile)
    bins = _bin_by_length(genes, bin_width)
    bins = [bin_ for bin_ in bins if len(bin_) > min_bin_size]
    return bins
