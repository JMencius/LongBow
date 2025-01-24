import pytest
from longbow.module.prediction_decode import decode


def test_qv_aspect():
    # Test for guppy software and 'qv' aspect
    assert decode(1, "guppy", "qv") == ("R9", "3or4", "FAST")
    assert decode(8, "guppy", "qv") == ("R10", "5or6", "SUP")

    # Test for dorado software and 'qv' aspect
    assert decode(3, "dorado", "qv") == ("R10", "FAST")


def test_mode_aspect():
    # Test for 'mode' aspect decoding
    assert decode(0, "guppy", "mode") == "HAC"
    assert decode(2, "dorado", "mode") == "FAST"


def test_edge_cases():
    # Test invalid aspect
    with pytest.raises(KeyError):
        decode(0, "guppy", "invalid_aspect")

    # Test invalid software
    with pytest.raises(KeyError):
        decode(1, "unknown_software", "qv")

