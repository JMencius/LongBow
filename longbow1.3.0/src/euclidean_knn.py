import math

def cal_euclidean_distance(autocorr1 : list, autocorr2 : list) -> float:
    # print(len(autocorr1), len(autocorr2))
    assert len(autocorr1) == len(autocorr2)
    # print(autocorr1)
    # print(autocorr2)
    return math.sqrt(sum([(autocorr1[i] - autocorr2[i]) ** 2 for i in range(len(autocorr1))])) 


def predict_hac_sup(subject : list, train_x, train_y, k = 3) -> int:
    # print(subject)
    edistance_list = list()
    for i in range(len(train_x)):
        clean_train_x = [float(j) for j in train_x[i]]
        dist = cal_euclidean_distance(subject, clean_train_x)
        edistance_list.append((list(train_y)[i], dist))

    edistance_list.sort(key = lambda s : s[1])
    # print(edistance_list)
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
