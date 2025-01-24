import pytest
import math
import numpy as np
from longbow.module.bhattacharyya import bhattacharyya


def test_bhattacharyya_normal_case():
    """
    Test the Bhattacharyya distance calculation for identical probability distributions.
    """
    # Example probability distributions with keys 0-93, where values sum to 1
    a = {i: 0 for i in range(94)}
    a[1] = 1

    b = a

    result = bhattacharyya(a, b)

    # Expect distance to be 0 for identical distributions
    assert result == 0, f"Expected 0, but got {result}"




def test_bhattacharyya_empty_distributions():
    """
    Test the Bhattacharyya distance when both distributions are empty.
    """
    a = {}
    b = {}
    
    with pytest.raises(ValueError):
        result = bhattacharyya(a, b)



