import numpy as np
from collections import Counter
import os
import sys
import argparse
import json
import pickle
from multiprocessing import Pool
import math
import time


def get_quality_score(qstring : str) -> float:
    assert type(qstring) == str
    score = 0
    count = 0
    for i in qstring:
        score += ord(i) - 33
        count += 1
    
    return score / count


def cal_quantile(in_list : list, Q : float) -> int:
    total = sum([i[1] for i in in_list])
    cut = total * Q
    count = 0
    subject = in_list[0][0]
    for i in in_list:
        if count > cut:
            return subject
        else:
            count += i[1]
            subject = i[0]
    return in_list[-1][0]


class KNN:
    def __init__(self, k=3):
        self.k = k

    def fit(self, X, y):
        self.X_train = X
        self.y_train = y

    def euclidean_distance(self, x1, x2):
        return np.sqrt(np.sum((x1 - x2) ** 2))

    def predict(self, X):
        y_pred = [self._predict(x) for x in X]
        return np.array(y_pred)

    def _predict(self, x):
        # 计算新数据点与所有参考数据点的距离
        distances = [self.euclidean_distance(x, x_train) for x_train in self.X_train]
        
        # 找到距离最近的k个数据点的索引
        k_indices = np.argsort(distances)[:self.k]
        
        # 找到这些数据点对应的标签
        k_nearest_labels = [self.y_train[i] for i in k_indices]
        
        # 选择最常见的标签作为预测结果
        most_common = Counter(k_nearest_labels).most_common(1)
        return most_common[0][0]


def is_sequence_line(line : str) -> bool:
    nucleo = {'A', 'T', 'G', 'C', 'N'}
    line = line.strip()
    for i in line:
        if i not in nucleo:
            return False
    return True



def process_chunck(filename : str, coreindex : int) -> None:
    base_level_result = dict()
    avg_read_level_result = list()
    linecount = 1
    with open(filename, 'r') as file:
        for line in file:
            if (linecount % 4) == 0:
                if (linecount // 4) % threads == coreindex:
                    # 去掉空格和换行符
                    m = line.strip()
                    read_score = 0
                    count = 0
                    for asci in m:
                        score = ord(asci) - 33
                        base_level_result[score] = base_level_result.get(score, 0) + 1
                        read_score += score
                        count += 1
                    if count != 0:
                        avg_read_level_result.append(read_score / count)
            linecount += 1

    """
    # Threading模块才用的全局变量来共享数据
    # 把运算结果返回给一个全局变量
    global results
    results[coreindex] = temp_result
    """
    # 使用multiprocessing时使用retun直接返回结果
    return (base_level_result, avg_read_level_result)
        


if __name__ == "__main__":
    start_time = time.time()

    version = "0.1"
    # 解析参数
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help = "Path to the input fastq file, including the fastq file name",required = True, type = str)
    parser.add_argument("-t", "--threads", help = "Number of threads", default = 12, type = int)
    parser.add_argument("-v", "--version", help = "Print software version info", action = "store_true")
    
    # print version info
    if "--version" in sys.argv[1 : ] or "-v" in sys.argv[1 : ]:
        print(f"Longbow version {version} based on Python 3.7+")
        sys.exit(0)
    
    args = parser.parse_args()

    threads = args.threads
    fastqfile = str(args.input)
    
    # print(type(threads), end = "\t")
    # print(threads)
    # print(type(fastqfile), end = "\t")
    print(fastqfile)

    if os.path.exists("KNN.pkl"):
        with open("KNN.pkl", "rb") as file:
            knn = pickle.load(file)
    
    else:

        train_data_dict = dict()
        with open("train_data.json", 'r') as file:
            data = json.load(file)
    
        convert = {"R9G2" : 0, "R9G3or4" : 1, "R9G5or6" : 2, "R10" : 3}
        x_data = list()
        y_data = list()

        for group in data:
            for species in data[group]:
                x_data.append(data[group][species])
                y_data.append(convert[group])

        X_train = np.array(x_data)
        y_train = np.array(y_data)    

        # 创建KNN分类器并拟合参考数据
        knn = KNN(k = 3)
        knn.fit(X_train, y_train)
        with open("KNN.pkl", "wb") as file:
            pickle.dump(knn, file)
    
    


    # 多线程计算fastq文件里面的质量分数
    # 用来保存结果的results词典
    
    
    # print(corerange)
    # print(pos)

    with Pool(threads) as p:
        output = p.starmap(process_chunck, [(fastqfile, coreindex) for coreindex in range(threads)])


    total_dict = dict()
    total_list = list()
    for d in output:
        for k in d[0].keys():
            if k not in total_dict:
                total_dict[k] = d[0][k]
            else:
                total_dict[k] += d[0][k]
    total_list = list(total_dict.items())
    total_list.sort(key = lambda K : K[0])

    max_avg_read_Q = 0
    for d in output:
        max_avg_read_Q = max(max_avg_read_Q, max(d[1]))

    print(total_list)
    print(f"max_avg_read_Q : {max_avg_read_Q}")


    quantile = [cal_quantile(total_list, 0), cal_quantile(total_list, 0.05), cal_quantile(total_list, 0.1), cal_quantile(total_list, 0.2),
                cal_quantile(total_list, 0.3), cal_quantile(total_list, 0.4), cal_quantile(total_list, 0.5), cal_quantile(total_list, 0.6),
                cal_quantile(total_list, 0.7), cal_quantile(total_list, 0.8), cal_quantile(total_list, 0.9), cal_quantile(total_list, 0.95), cal_quantile(total_list, 1)]
    print(quantile)


    # 预测新数据点的类别
    prediction = knn.predict([quantile])
    print(prediction)

    if prediction == 0:
        print("The predicted version is R9G2.")
    if prediction == 1:
        if max_avg_read_Q < 50:
            print("The predicted version is R9G3/4.")
        else:
            print("The predicted version is R9G5/6.")
    if prediction == 2:
        if max_avg_read_Q > 50:
            print("The predicted version is R9G5/6.")
        else:
            print("The predicted version is R9G3/4.") 
    if prediction == 3:
        print("The predicted version is R10G6.")
    end_time = time.time()

    print(f"Total time cost : {end_time - start_time} seconds")
    print("\n\n")


    

