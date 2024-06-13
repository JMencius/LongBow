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


def predict_knn(baseqv : list, train_x, train_y, k = 3) -> int:
    # preprocess change into dict
    
    normalized_baseqv = normalize(baseqv)
    baseqv_dict = {(i + 1) : normalized_baseqv[i] for i in range(90)}

    sim_list = list()
    
    for i in range(len(train_x)):
        train = list(train_x[i])
        # print(train)
        train_baseqv = {(j+1) : train[j] for j in range(90)}
        sim = cal_bhattacharyya_sim(baseqv_dict, train_baseqv)
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
    
