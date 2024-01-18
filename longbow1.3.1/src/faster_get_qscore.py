import os
import sys
from multiprocessing import Pool
import math
import pyfastx 
import numpy as np
from statsmodels.tsa.stattools import acf



def process_chunck(filename : str, coreindex : int, threads : int, qscore_cutoff : float, autocorr : bool) -> tuple:
    base_level_result = dict()
    itemcount = 1
    lags = range(1, 4)
    autocorr_summary = {i : [0, 0] for i in lags}

    for name, seq, qual in pyfastx.Fastq(filename, build_index = False):
        if itemcount % threads == coreindex:
            m = qual.strip()
            read_score_sum = 0
            count = 0
            temp_base = dict()
            base_qv = list()
            for asci in m:
                count += 1
                score = ord(asci) - 33
                base_qv.append(score)
                temp_base[score] = temp_base.get(score, 0) + 1
                converted = 10**(score / -10)
                read_score_sum += converted
                    
            read_score = -10 * math.log(read_score_sum / count, 10)    
            if read_score > qscore_cutoff:    
                for t in temp_base:
                    base_level_result[t] = base_level_result.get(t, 0) + temp_base[t]
                if autocorr:
                    max_lag = 3
                    if count > max_lag:
                        base_qv = np.array(base_qv)
                        autocor_value = acf(base_qv, nlags = 3)[1 : ]
                        for i in range(max_lag):
                            autocorr_summary[i + 1][0] += autocor_value[i] * (count - i - 1)
                            autocorr_summary[i + 1][1] += (count - i - 1)
                                
        itemcount += 1

    if autocorr:
        return (base_level_result, autocorr_summary)
    else:
        return (base_level_result)


def dict2sortlist(in_dict : dict) -> list:
    a = list(in_dict.items())
    a.sort(key = lambda K : K[0])
    return([i[1] for i in a])


def process_autocorr_dict(in_dict : dict) -> list:
    a = list(in_dict.items())
    a.sort(key = lambda K : K[0])
    return([i[1][0] / i[1][1] for i in a])


def get_qscore(fastqfile : str, threads : int, qscore_cutoff : float, autocorr : bool) -> tuple:
    # split and multiprocessing fastq file
    with Pool(threads) as p:
        output = p.starmap(process_chunck, [(fastqfile, coreindex, threads, qscore_cutoff, autocorr) for coreindex in range(threads)])
    
    # combine the result
    final_base_level = {i : 0 for i in range(1, 91)}
    final_corr = {i : [0, 0] for i in range(1, 4)}

    for i in output:
        for k in i[0]:
            if 1 <= k <= 90:
                final_base_level[k] += i[0][k]
        if autocorr:
            for j in i[1]:
                if 1 <= j <= 3:
                    final_corr[j][0] += i[1][j][0]
                    final_corr[j][1] += i[1][j][1]
    
    # print(final_base_level)
    # print(final_corr)


    sorted_base_list = dict2sortlist(final_base_level)
    sorted_autocorr_list = process_autocorr_dict(final_corr)
    if autocorr:
        return (sorted_base_list, sorted_autocorr_list)
    else:
        return (sorted_base_list)
    




