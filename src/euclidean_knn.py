import math
import numpy as np



def normalize(in_list : list) -> list:
    s = sum(in_list)
    assert s != 0, "Abnormal FASTQ format, Q score list sum is 0."
    normalized_list = [i/s for i in in_list]
    return normalized_list



def cal_euclidean_distance(autocorr1 : list, autocorr2 : list) -> float:
    assert len(autocorr1) == len(autocorr2)
    return math.sqrt(sum([(autocorr1[i] - autocorr2[i])**2 for i in range(len(autocorr1))])) 



def predict_hac_sup(subject : list, train_x, train_y, trim_lag : int, k = 3) -> int:
    # print(len(train_x))
    # print(len(train_y))

    edistance_list = list()
    for i in range(len(train_x)):
        clean_train_x = [float(j) for j in train_x[i]][: trim_lag]
        subject = subject[: trim_lag]
        dist = cal_euclidean_distance(subject, clean_train_x)
        edistance_list.append((list(train_y)[i], dist))

    edistance_list.sort(key = lambda s : s[1])
    top_k = edistance_list[ : k]
    label_k = [i[0] for i in top_k]
    # print(edistance_list)

    if len(set(label_k)) == len(label_k):
        return label_k[0]
    else:
        max_count_label = None
        max_count = 0
        for i in set(label_k):
            if label_k.count(i) > max_count:
                max_count = label_k.count(i)
                max_count_label = i
        return max_count_label

