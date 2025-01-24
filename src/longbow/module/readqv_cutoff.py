def cutoff_qv(readqv : dict) -> int:
    """
    Calculates the read QV cutoff from a given read QV distribution.

    Parameters:
    readqv (dict): A Python dictionary containing the read QV distribution for a sample.

    Returns:
    int: The read QV cutoff.

    Example:
    >>> readqv = {0: 0, 1: 0, 2: 0, 3: 0, 4: 10, 5: 20, ..., 93: 0}
    >>> cutoff_qv(readqv)
    4
    """

    cutoff_qv = 0
    for i in range(0, 94):
        if readqv[i] == 0:
            cutoff_qv = i + 1
        else:
            break

    return cutoff_qv


