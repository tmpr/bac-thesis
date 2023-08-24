from argparse import ArgumentParser

from pathlib import Path
from typing import Sequence

import numpy as np
import swalign
from networkx import Graph, connected_components
from numpy.typing import NDArray
from tqdm import tqdm

from blocks import get_blocks

# For simplicity, we handcode the one-hot encoding.
ONE_HOT = {
    "-": np.array([0, 0, 0, 0]),
    "A": np.array([1, 0, 0, 0]),
    "C": np.array([0, 1, 0, 0]),
    "G": np.array([0, 0, 1, 0]),
    "T": np.array([0, 0, 0, 1]),
}


def compute_blosum_matrix(blocks: Sequence[NDArray[np.character]], x: float = 0.62):
    """Compute BLOSUMx matrix.

    Args:
        blocks (Sequence[NDArray[np.character]]):
            Sequence of character matrices with the letters ACGT-.
        x (float, optional):
            Float between 0 and 1 - sequences with more or equal
            than x are being clustered. Defaults to 0.62.

    Returns:
        NDArray[int]: BLOSUMx scoring matrix for nucleotides.
    """
    counting_tables = []

    clustered_blocks = [
        _cluster_block(block, x)
        for block in tqdm(blocks, total=len(blocks), desc="Clustering blocks")
    ]
    counting_tables = [
        _compute_counting_table(block)
        for block in tqdm(
            clustered_blocks,
            total=len(clustered_blocks),
            desc="Computing counting_tables",
        )
    ]

    counting_table = (
        np.sum(counting_tables, axis=0)
        if len(counting_tables) > 1
        else counting_tables[0]
    )
    Q = _compute_Q(counting_table)
    P = _compute_p(Q)
    log_odds = _compute_log_odds(Q, P)
    # print(f"{counting_table=}\n, {Q=}, {P=}, {log_odds=}")
    return np.round(log_odds)


def _cluster_block(block: NDArray[np.character], x: float) -> NDArray:
    # (n sequences, n sequences)
    similarity_matrix = np.array(
        [[sum(seq_a == seq_b) / len(seq_a) for seq_a in block] for seq_b in block]
    )

    # If A and B are similar and B and C are similar, A and C get
    # clustered. This means that we can view the pairs of similar
    # sequences as edges in a graph, of which the resulting components
    # are to be clustered.
    similar_sequences_edges = list(zip(*np.where(similarity_matrix >= x)))

    # We include self edges to make sure all nodes are in the graph.
    self_edges = [(i, i) for i in range(len(block))]

    # TODO: Implement this yourself, a fullblown graph is probably overkill.
    similarity_graph = Graph(similar_sequences_edges + self_edges)
    components = list(connected_components(similarity_graph))

    if len(components) == 1:
        print("Block contains only similar sequences. Skipping.")
        # Return two sequences with only zeros to maintain typical shape.
        empty_clustered_block = np.zeros((2, 1, 4))
        return empty_clustered_block

    # (n sequences, n chars, 4 (n nucleotides))
    one_hot_block = [np.array([ONE_HOT[char] for char in seq]) for seq in block]
    clustered_block = np.array(
        [
            np.mean([one_hot_block[i] for i in component], axis=0)
            for component in components
        ]
    )
    return clustered_block


def _compute_counting_table(
    clustered_block: NDArray[np.character],
) -> NDArray[np.int32]:
    counting_table_assymmetric = np.sum(  # type: ignore
        # The outer prpoduct of the vectors computes
        # for each position the percentage of overlap # TODO
        [
            clustered_block[i].T @ clustered_block[j]
            for i in list(range(len(clustered_block)))
            for j in range(i)
        ],
        axis=0,
    )
    counting_table = counting_table_assymmetric + counting_table_assymmetric.T
    counting_table[np.diag_indices_from(counting_table)] /= 2
    counting_table = np.tril(counting_table)
    return counting_table


def _compute_Q(counting_table: NDArray):
    # Note that `counting_table` is a lower triangular matrix.
    Q = counting_table / np.sum(counting_table)
    for i in range(4):
        for j in range(i + 1, 4):
            Q[i, j] = Q[j, i]
    return Q


def _compute_p(Q: NDArray[np.float32]):
    return np.array(
        [
            Q[i, i] + sum(Q[i, j] / 2 for j in range(len(Q)) if i != j)
            for i in range(len(Q))
        ],
        dtype=np.float32,
    )


def _compute_log_odds(Q, p):
    return 2 * np.log2(
        np.array(
            [
                [
                    Q[i, j] / (p[i] ** 2) if i == j else Q[i, j] / (2 * p[i] * p[j])
                    for i in range(len(p))
                ]
                for j in range(len(p))
            ]
        )
    )


def _matrix_to_ssw_format(m: NDArray[np.int16]):
    return "\n".join(
        [
            "  A C G T",
            f"A {m[0, 0]} {m[0, 1]} {m[0, 2]} {m[0, 3]}",
            f"C {m[1, 0]} {m[1, 1]} {m[1, 2]} {m[1, 3]}",
            f"G {m[2, 0]} {m[2, 1]} {m[2, 2]} {m[2, 3]}",
            f"T {m[3, 0]} {m[3, 1]} {m[3, 2]} {m[3, 3]}",
        ]
    )


if __name__ == "__main__":
    arg_parser = ArgumentParser()
    arg_parser.add_argument("protein_family_code", type=str)
    arg_parser.add_argument(
        "x",
        type=float,
        help="Float between 0 and 1 - sequences"
        " with more or equal than x are being clustered",
        default=0.62,
    )
    arg_parser.add_argument("origin", type=str)
    args = arg_parser.parse_args()
    paths = list(Path("data/aligned/").glob(f"{args.protein_family_code}*"))
    blocks = []
    for path in paths:
        blocks.extend(get_blocks(path, min_column_density=1, min_block_len=20))

    Path(f"blocks/{args.protein_family_code}").write_text(
        "\n\n".join("\n".join("".join(row) for row in block) for block in blocks)
    )

    matrix = compute_blosum_matrix(blocks, args.x)
    text_matrix = _matrix_to_ssw_format(matrix)
    matrix_name = (
        f"BLOSUM{int(args.x*100)}_{args.protein_family_code}_{args.origin}.matrix"
    )
    matrix_path = f"matrices/{matrix_name}"
    Path(matrix_path).write_text(text_matrix)
