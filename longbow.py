#!/usr/bin/env python3

import time
import os
import sys
import argparse
import warnings
import json
from src.faster_get_qscore import get_qscore
from src.distinguish_software import guppy_or_dorado
from src.read_train import read_qv_train_file
from src.read_train import read_autocorr_train_file
from src.bhattacharyya_knn import predict_knn
from src.euclidean_knn import predict_hac_sup
from src.readqv_classificator import readqv_predict
from src.prediction_decode import decode


def main():
    start_time = time.time()
    warnings.simplefilter(action = "ignore", category = FutureWarning)
    version = ('1', '4', '1')
    script_dir = os.path.dirname(os.path.realpath(__file__))
    current_dir = os.path.dirname(os.path.realpath(__file__))

    # parse parameter
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help = "Path to the input fastq/fastq.gz file, including the fastq file name", required = True, type = str)
    parser.add_argument("-o", "--output", help = "Output directory or file name", default = "", type = str)
    parser.add_argument("-t", "--threads", help = "Number of parallel threads", default = 12, type = int)
    parser.add_argument("-q", "--qscore", help = "Read-level qscore filter", default = 10, type = int)
    parser.add_argument("-m", "--model", help = "Path to the training model csv data", default = os.path.join(current_dir, "model"), type = str)
    parser.add_argument("-c", "--corr", help = r"Do read-qv based correction or autocorrelation for hac/sup config", action = "store_true", default = True)
    parser.add_argument("-b", "--buf", help = "Output intermediate results of QV and autocorrelation", action = "store_true")
    parser.add_argument("-V", "--verbose", help = "Verbose mode, print the result", action = "store_true")
    parser.add_argument("-v", "--version", help = "Print software version info", action = "store_true")

    
    # print version info
    if "--version" in sys.argv[1 : ] or "-v" in sys.argv[1 : ]:
        print(f"Longbow version {'.'.join(version)} based on Python 3.7+")
        sys.exit(0)
    
    args = parser.parse_args()
    threads = args.threads
    fastqfile = args.input
    output_name = args.output
    model_path = args.model
    qscore_cutoff = args.qscore
    autocorr = args.corr
    buf = args.buf
    verbose = args.verbose
    
    # input not empty
    if os.path.getsize(fastqfile) == 0:
        print("Input file is empty!")
        sys.exit(1)

    # Get property from fastq file
    if autocorr:
        baseqv, autocorr, readqv = get_qscore(fastqfile, threads, qscore_cutoff, args.corr)
    else:
        baseqv, readqv = get_qscore(fastqfile, threads, qscore_cutoff, args.corr)



    # QV score distribution prediction
    pred_dict = {"Sample" : fastqfile, "Flowcell" : None, "Software" : None, "Version" : None, "Mode" : None}

    if guppy_or_dorado(baseqv) == "guppy":
        pred_dict["Software"] = "guppy"
        train_x, train_y = read_qv_train_file(f"{model_path}/train_guppy.csv")
        predict = predict_knn(baseqv, train_x, train_y)
        pred_dict["Flowcell"], pred_dict["Version"], pred_dict["Mode"] = decode(predict, "guppy", "qv")
    else:
        pred_dict["Software"] = "dorado"
        pred_dict["Version"] = '0'
        train_x, train_y = read_qv_train_file(f"{model_path}/train_dorado.csv")
        predict = predict_knn(baseqv, train_x, train_y)
        pred_dict["Flowcell"], pred_dict["Mode"] = decode(predict, "dorado", "qv")

    # autocorrelation
    if autocorr:
        autocorr_flag = 0
        if pred_dict["Software"] == "dorado" and pred_dict["Mode"] != "FAST":
            autocorr_flag = 1
            if pred_dict["Flowcell"] == "R9":
                model_csv_file = "dorado_autocorr_R9.csv"
                preset_lag = 10
            else:
                model_csv_file = "dorado_autocorr_R10.csv"
                preset_lag = 2

        elif pred_dict["Software"] == "guppy" and pred_dict["Version"] == '5or6' and pred_dict["Mode"] != "FAST":
            autocorr_flag = 1
            if pred_dict["Flowcell"] == "R9":
                model_csv_file = "guppy_autocorr_R9.csv"
                preset_lag = 10
            else:
                model_csv_file = "guppy_autocorr_R10.csv"
                preset_lag = 100
        
        if autocorr_flag:
            readqv_hac_sup = readqv_predict(readqv)
            if readqv_hac_sup != None:
                pred_dict["Mode"] = decode(readqv_hac_sup, "dorado", "hac_sup")
            else:
                train_x, train_y = read_autocorr_train_file(os.path.join(model_path, model_csv_file))
                hac_sup_pred = predict_hac_sup(autocorr, train_x, train_y, preset_lag)
                pred_dict["Mode"] = decode(hac_sup_pred, "dorado", "hac_sup")
    
    end_time = time.time()

    # if verbose mode or output not set, print the prediction result interactively
    if verbose or (not output_name):
        print(pred_dict)
    
    # output to json
    if output_name:
        with open(output_name, 'w') as json_file:
            json.dump(pred_dict, json_file, indent = 4, separators=(',', ': '))
            json_file.write("\n\n")

            # output intermediate results
            if buf:
                json.dump({"Run info" : {"LongBow version" : f"{'.'.join(version)}",
                                         "input" : fastqfile,
                                         "output" : output_name,
                                         "model" : model_path,
                                         "threads" : threads,
                                         "run time" : f"{end_time - start_time} s",
                                         "qscore cutoff" : qscore_cutoff,
                                         "autcorrelation" : bool(args.corr),
                                         "detail output" : bool(args.buf),
                                         "verbose mode" : bool(verbose)}}, json_file, indent = 4, separators=(',', ': '))
                json_file.write("\n\n")
                json.dump({"baseqv" : {i + 1 : baseqv[i] for i in range(len(baseqv))}}, json_file)
                json_file.write("\n\n")
                json.dump({"readqv" : {i + 1 : readqv[i] for i in range(len(readqv))}}, json_file)
                json_file.write("\n\n")
                if autocorr:
                    json.dump({"autocorrelation" : {i + 1 : autocorr[i] for i in range(len(autocorr))}}, json_file)
                    json_file.write("\n")

if __name__ == "__main__":
    main()



