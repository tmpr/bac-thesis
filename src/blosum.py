import pdb
from functools import cache
from typing import Sequence

import numpy as np
from networkx import Graph, connected_components
from numpy.typing import NDArray
from tqdm import tqdm

x = 0.62

# For the very simplicity, we handcode the one-hot encoding.
ONE_HOT = {
    "-": np.array([0, 0, 0, 0]),
    "A": np.array([1, 0, 0, 0]),
    "C": np.array([0, 1, 0, 0]),
    "G": np.array([0, 0, 1, 0]),
    "T": np.array([0, 0, 0, 1]),
}


def compute_blosum_matrix(
    blocks: Sequence[NDArray[np.character]], x: float = 0.62
):
    pairwise_counts = []
    for block in tqdm(
        blocks, total=len(blocks), desc="Computing BLOSUM matrix"
    ):
        # (n sequences, n sequences)
        similarity_matrix = np.array(
            [
                [sum(seq_a == seq_b) / len(seq_a) for seq_a in block]
                for seq_b in block
            ]
        )

        # If A and B are similar and B and C are similar, A and C get
        # clustered. This means that we can view the pairs of similar
        # sequences as edges in a graph, of which the resulting components
        # are to be clustered.
        similar_sequences_edges = list(zip(*np.where(similarity_matrix > x)))

        # We include self edges to make sure all nodes are in the graph.
        self_edges = [(i, i) for i in range(len(block))]

        similarity_graph = Graph(similar_sequences_edges + self_edges)
        components = list(connected_components(similarity_graph))

        if len(components) == 1:
            print(f"Block contains only similar sequences. Skipping.")
            continue

        # (n sequences, n chars, 4 (n nucleotides))
        one_hot_block = [
            np.array([ONE_HOT[char] for char in seq]) for seq in block
        ]
        clustered_block = [
            np.mean([one_hot_block[i] for i in component], axis=0)
            for component in components
        ]

        indices = list(range(len(clustered_block)))

        pairwise_counts.append(
            sum(
                # The outer product of the vectors gives us a 4x4 matrix
                # with the counts of the nucleotide pairs.
                clustered_block[i].T @ clustered_block[j]
                for i in indices
                for j in indices
                if i != j
            )
        )

    pairwise_count = (
        np.sum(pairwise_counts, axis=0)
        if len(pairwise_counts) > 1
        else pairwise_counts[0]
    )

    # The paper hides away that i <= j, i.e, the matrix and all steps are triangular.
    Q = pairwise_count / np.sum(np.tril(pairwise_count))
    P = np.array(
        [
            Q[i, i] + sum(Q[i, j] / 2 for j in range(4) if i != j)
            for i in range(4)
        ],
        dtype=np.float32,
    )

    anti_identity_matrix = -(np.eye(4) - 1)
    P_P = np.outer(P, P)
    E = P_P * np.eye(4) + 2 * P_P * anti_identity_matrix
    log_odds = np.log2(Q / E)

    scoring_matrix = np.round(log_odds) * 2

    return scoring_matrix
