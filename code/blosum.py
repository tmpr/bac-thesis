from typing import Sequence

import numpy as np
from networkx import Graph, connected_components
from sklearn.preprocessing import OneHotEncoder

block = ["AGT", "AGG", "ACC", "AGA"]
x = 0.62

ONE_HOT = {
    "A": np.array([1, 0, 0, 0]),
    "C": np.array([0, 1, 0, 0]),
    "G": np.array([0, 0, 1, 0]),
    "T": np.array([0, 0, 0, 1]),
}


def compute_blosum_matrix(blocks: Sequence[str], x: float = 62):
    enc = OneHotEncoder()
    pairwise_counts = []
    for block in blocks:
        if len(set(len(seq) for seq in block)) != 1:
            raise ValueError(
                f"Received block with uneven sequence length. {block=}"
            )
        block = np.array([[char for char in seq] for seq in block])
        sim_matrix = np.array(
            [
                [sum(seq_a == seq_b) / len(seq_a) for seq_a in block]
                for seq_b in block
            ]
        )
        to_merge = [(i, j) for i, j in zip(*np.where(sim_matrix > x)) if i >= j]
        one_hot_block = [
            np.array([ONE_HOT[char] for char in seq]) for seq in block
        ]
        components = list(connected_components(Graph(to_merge)))
        if len(components) == 1:
            raise ValueError(f"Block contains only similar sequences: {block=}")

        clustered_block = [
            np.mean([one_hot_block[i] for i in component], axis=0)
            for component in components
        ]

        # Counting_table
        indices = list(range(len(clustered_block)))
        pairwise_counts.append(
            sum(
                clustered_block[i].T @ clustered_block[j]
                for i in indices
                for j in indices
                if i != j
            )
        )

    pairwise_count = sum(pairwise_counts)
    count_sum = sum(
        pairwise_count[i, j] for i in range(4) for j in range(i + 1)
    )

    def q(i, j):
        return pairwise_count[i, j] / count_sum

    def p(i):
        return q(i, i) + sum(
            q(i, j) / 2 for i in range(4) for j in range(4) if i < j
        )

    log_odds = np.log2(
        np.array(
            [[q(i, j) / (p(i) * p(j)) for i in range(4)] for j in range(4)]
        )
    )

    scoring_matrix = np.round(log_odds) * 2

    return scoring_matrix
