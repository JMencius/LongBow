def readqv_predict(readqv : dict) -> int:
    pred = None
    if sum([readqv[i] for i in range(1, 11)]) == 0:
        pred = 1
    elif sum([readqv[i] for i in range(1, 10)]) == 0:
        pred = 0

    return pred
