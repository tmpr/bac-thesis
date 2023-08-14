import numpy as np
import pytest

from blosum import compute_blosum_matrix


def test_blosum():
    block = [
        np.array([np.array(list(s)) for s in ("AGT", "AGG", "ACC", "AGA")])
    ]
    compute_blosum_matrix(block)
