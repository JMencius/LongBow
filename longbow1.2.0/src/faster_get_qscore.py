import os
import sys
from multiprocessing import Pool
import math


def process_chunck(filename : str, coreindex : int, threads : int, qscore_cutoff : float) -> dict:
    base_level_result = dict()
    linecount = 1
    with open(filename, 'r') as file:
        for line in file:
            if (linecount % 4) == 0:
                if (linecount // 4) % threads == coreindex:
                    # 去掉空格和换行符
                    m = line.strip()
                    read_score_sum = 0
                    count = 0
                    temp_base = dict()
                    for asci in m:
                        count += 1
                        score = ord(asci) - 33
                        temp_base[score] = temp_base.get(score, 0) + 1
                        converted = 10 ** (score / -10)
                        read_score_sum += converted
                        read_score = -10 * math.log(read_score_sum / count, 10)
                    if read_score > qscore_cutoff:    
                        for t in temp_base:
                            base_level_result[t] = base_level_result.get(t, 0) + temp_base[t]
                                
            linecount += 1
    return base_level_result


def dict2sortlist(in_dict : dict) -> list:
    a = list(in_dict.items())
    a.sort(key = lambda K : K[0])
    return([i[1] for i in a])


def get_qscore(fastqfile : str, threads : int, qscore_cutoff : float) -> list:
    # split and multiprocessing fastq file
    with Pool(threads) as p:
        output = p.starmap(process_chunck, [(fastqfile, coreindex, threads, qscore_cutoff) for coreindex in range(threads)])
    
    # combine the result
    final_base_level = {i : 0 for i in range(1, 91)}
    for i in output:
        for k in i:
            final_base_level[k] += i[k]

    sorted_base_list = dict2sortlist(final_base_level)
    return sorted_base_list
    




