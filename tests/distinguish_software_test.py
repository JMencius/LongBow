import pytest
from longbow.module.distinguish_software import guppy_or_dorado


def dorado_test_case():
    # Test case 1: Dorado predicted with Q scores from 51 onward as 0
    base_qv = [0] + [10] * 50 + [0] * 43
    assert guppy_or_dorado(base_qv) == "dorado"


def guppy_test_case():
    # Test case 2: Guppy predicted with non-zero Q scores from 51 onward
    read_qv = [0] + [10] * 93
    assert guppy_or_dorado(read_qv) == "guppy"



def ancient_guppy_version_case():
    # Test case 3: Guppy predicted with Q score ceil at 20
    read_qv = [10] * 20 + [0] * 74
    assert guppy_or_dorado(read_qv) == "guppy"



def edge_cases():
    # Test case 4: Edge case with list length not equal to 94
    with pytest.raises(AssertionError):
        guppy_or_dorado([0] * 93)

    # Test case 5: Edge case with negative Q scores
    with pytest.raises(AssertionError):
        guppy_or_dorado([0] * 50 + [-1] + [0] * 43)
