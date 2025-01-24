import math
import numpy as np


def bhattacharyya(a: dict, b: dict) -> float:
    """
    Calculate the Bhattacharyya distance between two probability distributions.

    The function computes the Bhattacharyya distance between two input distributions
    represented by dictionaries and then converts it to a distance metric. A smaller
    Bhattacharyya distance indicates that the distributions are more similar.

    Parameters:
    a (dict): The first probability distribution (key-value pairs, where keys are the items
              and values are the probabilities).
    b (dict): The second probability distribution (key-value pairs, where keys are the items
              and values are the probabilities).

    Returns:
    float: The Bhattacharyya distance between the two input distributions.
    """
    bhattacharyya_coefficient = 0
    for k in a:
        if k in b:
            bhattacharyya_coefficient += math.sqrt(a[k] * b[k])
    
    if bhattacharyya_coefficient > 0:
        bhattacharyya_distance = -math.log(bhattacharyya_coefficient)
    else:
        raise ValueError("Bhattacharyya coefficient < 0, please check the distribution")
    
    if np.isinf(bhattacharyya_distance):
        return 0


    return bhattacharyya_distance
