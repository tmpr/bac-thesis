import numpy as np
import pytest
from numpy.typing import NDArray

from blosum import (
    _cluster_block,
    _compute_counting_table,
    _compute_log_odds,
    _compute_p,
    _compute_Q,
    compute_blosum_matrix,
)


@pytest.fixture()
def hochreiter_example() -> NDArray[np.character]:
    # G <=> V, A <=> I, C <=> L
    return np.array(
        [
            ["G", "A", "A", "C"],
            ["G", "A", "A", "A"],
            ["C", "A", "G", "G"],
            ["C", "C", "G", "A"],
            ["G", "C", "C", "A"],
            ["G", "C", "C", "C"],
        ],
        dtype="<U1",
    )


def test_cluster_blocks(hochreiter_example):
    clustered = _cluster_block(hochreiter_example, x=0.75)
    expected = np.array(
        [
            [
                # First block has as last letter 0.5 A 0.5 C
                [0.0, 0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0, 0.0],
                [0.5, 0.5, 0.0, 0.0],
            ],
            [
                [0.0, 1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
            ],
            [
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0, 0.0],
            ],
            [
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.5, 0.5, 0.0, 0.0],
            ],
        ]
    )
    assert np.all(clustered == expected)


def test_compute_counting_table(hochreiter_example):
    clustered_block = _cluster_block(hochreiter_example, 0.75)
    counting_table = _compute_counting_table(clustered_block)
    expected = np.array(
        [
            [2.25, 0.0, 0.0, 0.0],
            [6.5, 2.25, 0.0, 0.0],
            [4.0, 7.0, 2.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
    )
    assert np.all(counting_table == expected)


def test_compute_Q(hochreiter_example):
    counting_table = _compute_counting_table(
        _cluster_block(hochreiter_example, 0.75)
    )
    Q = _compute_Q(counting_table)
    expected = np.array(
        [
            [0.094, 0.271, 0.167, 0.0],
            [0.271, 0.094, 0.292, 0.0],
            [0.167, 0.292, 0.083, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
    )
    assert np.all(np.round(Q, 3) == expected)


def test_compute_p(hochreiter_example):
    Q = _compute_Q(
        _compute_counting_table(_cluster_block(hochreiter_example, 0.75))
    )
    p = _compute_p(Q)
    expected = np.array([0.3125, 0.375, 0.3125, 0.0])
    assert np.all(p == expected)


def test_compute_log_odds(hochreiter_example):
    Q = _compute_Q(
        _compute_counting_table(_cluster_block(hochreiter_example, 0.75))
    )
    p = _compute_p(Q)
    # For stability reasons, we cut away the T entries.
    Q, p = Q[:3, :3], p[:3]
    log_odds = _compute_log_odds(Q, p)
    expected = np.array(
        [
            [-0.12, 0.42, -0.46],
            [0.42, -1.17, 0.63],
            [-0.46, 0.63, -0.46],
        ]
    )
    assert np.all(np.round(log_odds, 2) == expected)


def test_compute_blosum(hochreiter_example):
    blosum_75 = compute_blosum_matrix([hochreiter_example], 0.75)
    expected = np.array(
        [
            [0, 0, 0],
            [0, -1, 1],
            [0, 1, 0],
        ]
    )
    assert np.all(blosum_75[:3, :3] == expected)
