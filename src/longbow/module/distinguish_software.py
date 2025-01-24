def guppy_or_dorado(in_list : list) -> str:
    """
    Determines the basecalling software ('guppy' or 'dorado') based on Q score distribution.

    If all Q scores from Q51 onward are 0, and Q50 is not 0,
    the software is predicted as "dorado"; otherwise, it is predicted as "guppy".

    Parameters:
    in_list (list): A list of 94 integers representing the Q scores.

    Returns:
    str: The predicted basecalling software ("guppy" or "dorado").

    Raises:
    AssertionError: If the input list does not have length 94 or contains negative values.
    """

    assert len(in_list) == 94, "Q score list length does not equal 90"
    assert min(in_list) >= 0, "Q score list has negative numbers, check the input fastq file"

    if sum(in_list[51: ]) == 0 and in_list[50] != 0:
        software = "dorado"
    else:
        software = "guppy"

    return software
