import math


def normalize(in_list : list) -> list:
    """
    Normalizes the input list by dividing each element by the sum of all elements.

    Args:
        in_list (list): A list of numerical values.

    Returns:
        list: A normalized list where the sum of elements equals 1.
    """

    s = sum(in_list)
    assert s != 0, "Abnormal FASTQ format, Q score list sum is 0."
    normalized_list = [i/s for i in in_list]
    return normalized_list



def cal_euclidean_distance(autocorr1 : list, autocorr2 : list) -> float:
    """
    Calculates the Euclidean distance between two lists.

    Args:
        autocorr1 (list): First list of autocorrelation values.
        autocorr2 (list): Second list of autocorrelation values.

    Returns:
        float: The Euclidean distance between the two lists.
    """

    assert len(autocorr1) == len(autocorr2)
    return math.sqrt(sum([(autocorr1[i] - autocorr2[i])**2 for i in range(len(autocorr1))])) 



def predict_mode(subject : list, train_x : list, train_y : list, trim_lag : int, k : int = 3) -> int:
    """
    Predicts the most commom class label (mode) from the k nearest neighbors 
    based on Euclidean distance, using the provided training data.
    
    Args:
        subject (list): The query list of values to classify.
        train_x (list): The list of training data points (features).
        train_y (list): The list of class labels for each training data point.
        trim_lag (int): The number of initial features to consider from each data point.
        k (int, optional): The number of nearest neighbors to consider. Default is 3.

    Returns:
        tuple: A tuple containing the predicted label and the weighted confidence score.
    """


    edistance_list = list()
    for i in range(len(train_x)):
        clean_train_x = [float(j) for j in train_x[i]][: trim_lag]
        subject = subject[: trim_lag]
        dist = cal_euclidean_distance(subject, clean_train_x)
        edistance_list.append((list(train_y)[i], dist))

    edistance_list.sort(key = lambda s : s[1])
    top_k = edistance_list[ : k]
    label_k = [i[0] for i in top_k]
    #print(edistance_list)

    if len(set(label_k)) == len(label_k):
        pred_mode = label_k[0]
    else:
        max_count_label = None
        max_count = 0
        for i in set(label_k):
            if label_k.count(i) > max_count:
                max_count = label_k.count(i)
                max_count_label = i
        pred_mode = max_count_label
    
    # calculate weighted confidence
    total_weight = sum([1 / e[1] for e in edistance_list])
    predict_weight = 0
    for i in edistance_list:
        if i[0] == pred_mode:
            predict_weight += 1 / i[1]

    weighted_confidence = predict_weight / total_weight

    return (pred_mode, weighted_confidence)

