import pytest
from longbow.module.readqv_cutoff import cutoff_qv


def test_cutoff_qv():
    # Test case 1: Standard case with a known cutoff value
    readqv = {i: 0 for i in range(0, 94)}
    readqv[4] = 10
    assert cutoff_qv(readqv) == 4


def test_edge_case():
    # Test case 2: Case where no QVs are zero
    readqv = {i: 10 for i in range(0, 94)}
    assert cutoff_qv(readqv) == 0


