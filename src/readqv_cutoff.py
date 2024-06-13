def cutoff_qv(readqv : dict) -> int:
    cutoff_qv = 0
    for i in range(1, 91):
        if readqv[i - 1] == 0:
            cutoff_qv = i
        else:
            break
    return cutoff_qv
