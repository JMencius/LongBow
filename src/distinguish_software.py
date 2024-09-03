# quality control and basecalling software distinguish
def guppy_or_dorado(in_list : list) -> str:
    assert len(in_list) == 94, "Q score list length does not equal 90"
    assert min(in_list) >= 0, "Q score list has negative numbers, check the input fastq file"

    if sum(in_list[51: ]) == 0 and in_list[50] != 0:
        software_predict = "dorado"
    else:
        software_predict = "guppy"

    return software_predict
