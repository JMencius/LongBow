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



def process_chunck(filename : str, coreindex : int) -> tuple:
    base_level_result = dict()
    avg_base_level_result = list()
    read_level_result = list()


    linecount = 1
    with open(filename, 'r') as file:
        for line in file:
            if (linecount % 4) == 0:
                if (linecount // 4) % threads == coreindex:
                    # 去掉空格和换行符
                    m = line.strip()
                    score_sum = 0
                    read_score_sum = 0
                    count = 0
                    for asci in m:
                        score = ord(asci) - 33
                        base_level_result[score] = base_level_result.get(score, 0) + 1
                        score_sum += score
                        read_score_sum += 10 ** ((ord(asci) - 33)/-10)
                        count += 1
                    if count != 0:
                        avg_base_level_result.append(score_sum / count)
                        read_level_result.append(-10 * math.log(read_score_sum / count, 10))
            linecount += 1

    # 使用multiprocessing时使用return直接返回结果
    return (base_level_result, avg_base_level_result, read_level_result)
        


def train_KNN_model(trainfile : str):
    with open(trainfile, 'r') as file:
            data = json.load(file)

    convert = {"R9G2" : 0, "R9G3or4" : 1, "R9G5or6" : 2, "R10" : 3}
    x_data = list()
    y_data = list()

    for group in data:
        for species in data[group]:
            x_data.append(data[group][species])
            y_data.append(convert[group])

    X_train = np.array(x_data)
    Y_train = np.array(y_data)

    # 创建KNN分类器并拟合参考数据
    knn = KNN(k = 3)
    knn.fit(X_train, Y_train)
    with open(args.json + ".KNN.pkl", "wb") as file:
        pickle.dump(knn, file)

    return knn

def process_fastq(fastqfile : str, threads : int) -> tuple:
    # 多线程计算fastq文件里面的质量分数
    # 用来保存结果的results词典
    with Pool(threads) as p:
        output = p.starmap(process_chunck, [(fastqfile, coreindex) for coreindex in range(threads)])
    
    # split and assemble result
    per_base_total_dict = dict()
    for d in output:
        for k in d[0].keys():
            if k not in per_base_total_dict:
                per_base_total_dict[k] = d[0][k]
            else:
                per_base_total_dict[k] += d[0][k]
    per_base_total_list = list(per_base_total_dict.items())
    per_base_total_list.sort(key = lambda K : K[0])

    max_avg_per_base_Q = 0
    for d in output:
        max_avg_per_base_Q = max(max_avg_per_base_Q, max(d[1]))
    
    min_per_read_Q = None
    for d in output:
        if min_per_read_Q == None:
            min_per_read_Q = min(d[2])
        else:
            min_per_read_Q = min(min_per_read_Q, min(d[2]))

    return (per_base_total_list, max_avg_per_base_Q, min_per_read_Q)



def judgement(model, quantile, baseQ_list, max_avg_per_base_Q, min_per_read_Q) -> tuple:
    software, flowcell, guppy_version, basecalling_mode = None, None, None, None
    """    
    return tuple format: 
    (0, 1, 2, 3)
    0 : Basecalling software, Albacore or Guppy or Dorado ;
    1 : Nanopore flowcell version, R9, R10, NA(due to lack of Dorado training data);
    2 : Guppy main version, Guppy2, Guppy3/4, Guppy5/6, NA;
    3 : Basecalling mode, FAST, HAC, SUP, MIXED, NA;
    """
    if 40 <= max([i[0] for i in baseQ_list]) <= 50:
        software = "Dorado"
        return (software, flowcell, guppy_version, basecalling_mode)
    elif max([i[0] for i in baseQ_list]) < 40:
        software = "Albacore"
        flowcell = "R9"
        return (software, flowcell, guppy_version, basecalling_mode)
    else:
        software = "Guppy"
        
    prediction = model.predict([quantile])
    if prediction == 0:
        flowcell = "R9"
        guppy_version = "Guppy2"
    elif prediction == 3:
        flowcell = "R10"
        guppy_version = r"Guppy5/6"
    else:
        flowcell = "R9"
        if max_avg_per_base_Q > 50:
            guppy_version = r"Guppy5/6"
        else:
            guppy_version = r"Guppy3/4"
    print(min_per_read_Q)
    if guppy_version == r"Guppy5/6":
        if min_per_read_Q < 5:
            basecalling_mode = "MIXED"
        if 7.5 < min_per_read_Q < 8.5:
            basecalling_mode = "FAST"
        elif 8.5 < min_per_read_Q < 9.5:
            basecalling_mode = "HAC"
        elif 9.5 < min_per_read_Q:
            basecalling_mode = "SUP"


    return (software, flowcell, guppy_version, basecalling_mode)



if __name__ == "__main__":
    start_time = time.time()

    version = (0, 2, 1)
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # 解析参数
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help = "Path to the input fastq file, including the fastq file name", required = True, type = str)
    parser.add_argument("-o", "--output", help = "Output directory or file name, default output to current directory", default = ".KNN.txt", type = str)
    parser.add_argument("-t", "--threads", help = "Number of threads", default = 12, type = int)
    parser.add_argument("-v", "--version", help = "Print software version info", action = "store_true")
    parser.add_argument("-l", "--label", help = "Data label", default = "NA", type = str)
    parser.add_argument("-j", "--json", help = "Path to the training data json input; default : 0 means used pretrained model", default = '0', type = str)


    # print version info
    if "--version" in sys.argv[1 : ] or "-v" in sys.argv[1 : ]:
        print(f"Longbow version {'.'.join(version)} based on Python 3.7+")
        sys.exit(0)
    
    args = parser.parse_args()

    threads = args.threads
    fastqfile = str(args.input)
    trainfile = str(args.json)
    output_name = os.path.basename(fastqfile) + str(args.output)
    label = str(args.label)
    

    # print(fastqfile)
    
    if trainfile == '0':
        with open(os.path.join(script_dir, "bin", "KNN.pkl"), "rb") as pklfile:
            model = pickle.load(pklfile)
    
    else:
        model = train_KNN_model(trainfile)

    # process fastq file with multiprocessing and assemble the result into a dict 
    temp_result = process_fastq(fastqfile, threads)
    baseQ_list, max_avg_per_base_Q, min_per_read_Q = temp_result
    

    # calculate quantile
    quantile = [cal_quantile(baseQ_list, 0), cal_quantile(baseQ_list, 0.05), cal_quantile(baseQ_list, 0.1), cal_quantile(baseQ_list, 0.2),
                cal_quantile(baseQ_list, 0.3), cal_quantile(baseQ_list, 0.4), cal_quantile(baseQ_list, 0.5), cal_quantile(baseQ_list, 0.6),
                cal_quantile(baseQ_list, 0.7), cal_quantile(baseQ_list, 0.8), cal_quantile(baseQ_list, 0.9), cal_quantile(baseQ_list, 0.95), cal_quantile(baseQ_list, 1)]
    # print(quantile)
    final_result = judgement(model, quantile, baseQ_list, max_avg_per_base_Q, min_per_read_Q)
    

    end_time = time.time()
    # write to output txt
    with open(output_name, 'w') as f:
        f.write(f"Processed file : {fastqfile}" + '\n')
        f.write(f"File size is : {os.path.getsize(fastqfile) / 1024**3} G" + '\n')
        f.write(f"Label is : {label}" + '\n')
        f.write(f"Q score distribution list is : {baseQ_list}" + '\n')
        f.write(f"Q score quantile is : {quantile}" + '\n')
        f.write(f"Predicted software is : {final_result[0]}" + '\n')
        f.write(f"Predicted flowcell version is : {final_result[1]}" + '\n')
        f.write(f"Predicted subversion is : {final_result[2]}" + '\n')
        f.write(f"Predicted basecalling mode is : {final_result[3]}" + '\n')
        f.write(f"Total time cost : {end_time - start_time} seconds")
    
    print(f"{fastqfile} process completed")
