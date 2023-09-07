from __future__ import annotations
from collections import Counter
from itertools import pairwise
from pathlib import Path
from typing import List, Tuple


from Bio import SeqIO
import numpy as np
from numpy.typing import NDArray


def load_msa(path: Path) -> NDArray[np.generic]:
    alignment = np.array(list(SeqIO.parse(path, format="fasta")))
    return alignment


def get_dense_column_idxs(
    alignment: NDArray[np.generic], min_column_density: float
) -> List[int]:
    dense_column_idx = [
        i
        for i, column in enumerate(alignment.T)
        if Counter(column)["-"] == (1 - min_column_density)
    ]
    return dense_column_idx


def get_contiguous_intervals(indices: List[int]) -> List[Tuple[int, int]]:
    intervals = []
    indices = list(set(sorted(indices)))
    beginning = indices[0]
    for former, current in pairwise(indices[1:]):
        if former == current - 1:
            continue
        else:
            intervals.append((beginning, former))
            beginning = current
    return intervals


def get_blocks(msa_path: Path, min_block_len: int = 5, min_column_density: float = 0.5):
    alignment = load_msa(msa_path)
    dense_column_idxs = get_dense_column_idxs(alignment, min_column_density)
    if not dense_column_idxs:
        return []
    intervals = get_contiguous_intervals(dense_column_idxs)
    intervals = [
        (start, end) for start, end in intervals if start <= end - min_block_len
    ]
    blocks = [alignment[:, start:end] for start, end in intervals]
    blocks = [block for block in blocks if np.all(block != "-")]
    return blocks
