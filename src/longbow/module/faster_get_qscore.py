import os
import sys
from multiprocessing import Pool
import math
import pyfastx 
import numpy as np
from statsmodels.tsa.stattools import acf


def process_chunck(filename : str, coreindex : int, threads : int, qscore_cutoff : int, autocorr : bool) -> tuple:
    """
    Process a chunk of a given FASTQ file, calculating base quality scores, read quality scores,
    and optionally performing autocorrelation analysis.

    Parameters:
    filename (str): Path to the FASTQ file.
    coreindex (int): Index of the current processing thread for task distribution.
    threads (int): Total number of threads.
    qscore_cutoff (int): Quality score cutoff. Reads with read QV below this score are ignored.
    autocorr (bool): Whether to calculate autocorrelation. If True, autocorrelation data is returned.

    Returns:
    tuple: A tuple containing the following elements:
        - base_level_result (dict): Frequency of each base quality score.
        - autocorr_summary (dict, optional): Autocorrelation results, returned only if autocorr is True.
        - read_level_result (dict): Frequency of each read quality score.
        - pass_reads (int): Number of reads passing the quality score cutoff.
        - qv_outliner (int): Number of bases with outlier quality values.
    """
    base_level_result = dict()
    read_level_result = dict()
    pass_reads = 0
    itemcount = 1
    lags = range(1, 101)
    autocorr_summary = {i : [0, 0] for i in lags}
    prophet = {chr(i + 33) : i for i in range(0, 94)}
    prophet_score = {chr(i + 33) : 10**(i / -10) for i in range(0, 94)}
    qv_outliner = 0
    
    try:
        fastq_file = pyfastx.Fastq(filename, build_index = False)
    except:
        print("FAIL to open input FASTQ file, check input FASTQ")
        sys.exit(1)
    

    for name, seq, qual in fastq_file:
        if itemcount % threads == coreindex:
            m = qual.strip()
            read_score_sum = 0
            count = 0
            temp_base = dict()
            base_qv = list()
            for asci in m:
                if asci in prophet:
                    count += 1
                    score = prophet[asci]
                    base_qv.append(score)
                    temp_base[score] = temp_base.get(score, 0) + 1
                    converted = prophet_score[asci]
                    read_score_sum += converted
                else:
                    qv_outliner += 1
            # skip if no base
            if count == 0:
                continue
            read_score = -10 * math.log(read_score_sum / count, 10)    
            read_qv = round(read_score)
            read_level_result[read_qv] = read_level_result.get(read_qv, 0) + 1
            
            if read_score > qscore_cutoff:
                for t in temp_base:
                    base_level_result[t] = base_level_result.get(t, 0) + temp_base[t]
                if autocorr:
                    max_lag = 100
                    if len(base_qv) > max_lag:
                        base_qv = np.array(base_qv)
                        try:
                            autocor_value = acf(base_qv, nlags = max_lag)[1 : ]
                            for i in range(max_lag):
                                autocorr_summary[i + 1][0] += autocor_value[i] * (count - i - 1)
                                autocorr_summary[i + 1][1] += (count - i - 1)
                            pass_reads += 1
                        except:
                            pass
                else:
                    pass_reads += 1
                            
        itemcount += 1

    if autocorr:
        return (base_level_result, autocorr_summary, read_level_result, pass_reads, qv_outliner)
    else:
        return (base_level_result, read_level_result, pass_reads, qv_outliner)


def dict2sortlist(in_dict : dict) -> list:
    """
    Convert a dictionary to a sorted list based on the keys.

    Parameters:
    in_dict (dict): The dictionary to be sorted.

    Returns:
    list: A list of values sorted by their corresponding keys.
    """

    a = list(in_dict.items())
    a.sort(key = lambda K : K[0])
    return([i[1] for i in a])


def process_autocorr_dict(in_dict : dict) -> list:
    """
    Process an autocorrelation dictionary to compute the autocorrelation values.

    Parameters:
    in_dict (dict): The dictionary containing autocorrelation values.

    Returns:
    list: A list of processed autocorrelation values.
    """

    a = list(in_dict.items())
    a.sort(key = lambda K : K[0])
    autocorr_result = []
    for i in a:
        if i[1][1] == 0:
            autocorr_result.append(0)
        else:
            autocorr_result.append(i[1][0] / i[1][1])
    return autocorr_result


def get_qscore(fastqfile : str, threads : int, qscore_cutoff : float, autocorr : bool) -> tuple:
    """
    Process the FASTQ file in parallel, calculate base and read quality scores, and optionally
    calculate autocorrelation for the reads.

    Parameters:
    fastqfile (str): Path to the FASTQ file.
    threads (int): Number of threads for parallel processing.
    qscore_cutoff (float): Quality score cutoff. Reads below this score are ignored.
    autocorr (bool): Whether to calculate autocorrelation. If True, autocorrelation data is returned.

    Returns:
    tuple: A tuple containing:
        - sorted_base_list (list): Sorted list of base quality scores.
        - sorted_autocorr_list (list, optional): Sorted list of autocorrelation values, returned only if autocorr is True.
        - sorted_read_list (list): Sorted list of read quality scores.
        - total_outliner (int): Total number of bases with outlier quality values.
    """

    # split and multiprocessing fastq file
    with Pool(threads) as p:
        output = p.starmap(process_chunck, [(fastqfile, coreindex, threads, qscore_cutoff, autocorr) for coreindex in range(threads)])
    
    # combine the result
    final_base_level = {i : 0 for i in range(0, 94)}
    final_read_level = {i : 0 for i in range(1, 94)}
    final_corr = {i : [0, 0] for i in range(1, 101)}
    total_pass = 0
    total_outliner = 0

    for i in output:
        for k in i[0]:
            if 0 <= k <= 93:
                final_base_level[k] += i[0][k]
        if autocorr:
            for j in i[1]:
                if 1 <= j <= 100:
                    if not(math.isnan(i[1][j][0])):
                        final_corr[j][0] += i[1][j][0]
                        final_corr[j][1] += i[1][j][1]
            for l in i[2]:
                if 1 <= l <= 93:
                    final_read_level[l] += i[2][l]
            total_pass += i[3]
            total_outliner += i[4]
        else:
            for l in i[1]:
                if 1 <= l <= 93:
                    final_read_level[l] += i[1][l]
            total_pass += i[2]
            total_outliner += i[3]
    
    if total_pass == 0:
        print("Abnormal FASTQ format, please check input FASTQ")
        sys.exit(1)


    sorted_base_list = dict2sortlist(final_base_level)
    sorted_read_list = dict2sortlist(final_read_level)
    if autocorr:
        sorted_autocorr_list = process_autocorr_dict(final_corr)
        return (sorted_base_list, sorted_autocorr_list, sorted_read_list, total_outliner)
    else:
        return (sorted_base_list, sorted_read_list, total_outliner)
    

