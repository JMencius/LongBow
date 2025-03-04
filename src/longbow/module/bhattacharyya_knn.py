import numpy as np
import math
try:
    from longbow.module.bhattacharyya import bhattacharyya
except ImportError:
    from module.bhattacharyya import bhattacharyya



def cal_bhattacharyya_sim(distri1 : dict, distri2 : dict) -> float:
    """
    Calculate the Bhattacharyya distance to represent the similarity between two probability distributions.

    Parameters:
    distri1 (dict): The first probability distribution represented as a dictionary.
    distri2 (dict): The second probability distribution represented as a dictionary.

    Returns:
    float: The Bhattacharyya distance value.
    """

    bhattacharyya_distance = bhattacharyya(distri1, distri2)
    return bhattacharyya_distance


def normalize(in_list : list) -> list:
    """
    Normalize a list of values by dividing each element by the sum of the list.
    This function ensures that the sum of the elements in the list becomes 1.

    Parameters:
    in_list (list): A list of numerical values to be normalized.

    Returns:
    list: A list of normalized values.
    
    Raises:
    AssertionError: If the sum of the input list is 0.
    """

    s = sum(in_list)
    assert s != 0, "Abnormal FASTQ format, Q score list sum is 0."
    normalized_list = [i/s for i in in_list]
    return normalized_list



def predict_knn(baseqv : list, train_x : list, train_y : list, software : str, k : int = 3) -> tuple:
    """
    Predict the label of a sample using k-nearest neighbors (k-NN) and compute a weighted confidence score.

    The function calculates the Bhattacharyya similarity between the input sample and each training sample,
    performs a majority vote among the top k nearest neighbors, and computes a confidence score based on 
    the inverse similarity (closer neighbors have more influence).

    Parameters:
    baseqv (list): Base QV count of the input sample (length 94).
    train_x: Training input samples (list of lists, each of length 94).
    train_y: Labels corresponding to the training samples.
    software: basecalling software (guppy or dorado)
    k (int, optional): Number of nearest neighbors to consider (default is 3).

    Returns:
    tuple: (predicted_label, weighted_confidence)
        predicted_label (int): Predicted label based on k-NN.
        weighted_confidence (float): Confidence score based on the weighted average of the top k neighbors' Bhattacharyya distance.
    """

    # preprocess change into dict
    normalized_baseqv = normalize(baseqv)
    baseqv_dict = {(i + 1) : normalized_baseqv[i] for i in range(94)}

    sim_list = list()
    
    for i in range(len(train_x)):
        train = list(train_x[i])
        train_baseqv = {(j+1) : train[j] for j in range(94)}
        sim = cal_bhattacharyya_sim(baseqv_dict, train_baseqv)
        sim_list.append((list(train_y)[i], sim))

    sim_list.sort(key = lambda s : s[1])
    top_k = sim_list[ : k]    
    label_k = [i[0] for i in top_k]
    
    
    if len(set(label_k)) == len(label_k):
        predicted_label = label_k[0]
    else:
        max_count_label = None
        max_count = 0
        for i in set(label_k):
            if label_k.count(i) > max_count:
                max_count = label_k.count(i)
                max_count_label = i
        predicted_label = max_count_label
    
    # calculate weighted confidence
    if software == "guppy":
        same_group_labels = ((0,), (1, 2), (3, 4, 5), (6, 7, 8))
    elif software == "dorado":
        same_group_labels = ((0, 1, 2), (3, 4, 5))
    
    for i in same_group_labels:
        if predicted_label in i:
            group = i

    total_weight = sum([1 / sim[1] for sim in sim_list])
    predict_weight = 0
    for i in sim_list:
        if i[0] in group:
            predict_weight += 1 / i[1]

    weighted_confidence = predict_weight / total_weight
    return (predicted_label, weighted_confidence)
    


