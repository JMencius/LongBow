import math
import numpy as np
from dictances import bhattacharyya


def cal_bhattacharyya_sim(distri1 : dict, distri2 : dict) -> float:

    bhattacharyya_similarity = bhattacharyya(distri1, distri2)
    return bhattacharyya_similarity

def normalize(in_list : list) -> list:
    s = sum(in_list)
    assert s != 0, "Q score list sum is 0."
    normalized_list = [i/s for i in in_list]
    return normalized_list


def cal_euclidean_distance(autocorr1 : list, autocorr2 : list) -> float:
    assert len(autocorr1) == len(autocorr2)
    return math.sqrt(sum([(autocorr1[i] - autocorr2[i])**2 for i in range(len(autocorr1))])) 


def predict_hac_sup(subject : list, train_x, train_y, trim_lag : int, k = 3) -> int:
    #print(subject)
    edistance_list = list()
    for i in range(len(train_x)):
        clean_train_x = [float(j) for j in train_x[i]][: trim_lag]
        subject = subject[: trim_lag]
        dist = cal_euclidean_distance(subject, clean_train_x)
        edistance_list.append((list(train_y)[i], dist))

    edistance_list.sort(key = lambda s : s[1])
    #print("euclid:")
    #print(edistance_list)
    top_k = edistance_list[ : k]
    label_k = [i[0] for i in top_k]

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

def predict_amp_hac_sup(readqv : list, train_x, train_y, k = 3) -> int:
    # preprocess change into dict
    
    normalized_readqv = normalize(readqv)
    readqv_dict = {(i + 10) : normalized_readqv[i] for i in range(0, 81)}

    sim_list = list()
    
    for i in range(len(train_x)):
        train = list(train_x[i])
        train_readqv = {(j + 10) : train[j] for j in range(0, 81)}
        sim = cal_bhattacharyya_sim(readqv_dict, train_readqv)
        sim_list.append((list(train_y)[i], sim))

    sim_list.sort(key = lambda s : s[1])
    top_k = sim_list[ : k]    
    label_k = [i[0] for i in top_k]
    # print(sim_list)
    
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
