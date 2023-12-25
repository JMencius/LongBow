#!/usr/bin/env python3

import time
import os
import sys
import argparse
import warnings
import json
from src.faster_get_qscore import get_qscore
from src.distinguish_software import guppy_or_dorado
from src.read_train import read_train_file
from src.knn import predict_knn
from src.prediction_decode import decode

def main():
    warnings.simplefilter(action = "ignore", category = FutureWarning)
    version = ('1', '2', '1')
    script_dir = os.path.dirname(os.path.realpath(__file__))
    current_dir = os.path.dirname(os.path.realpath(__file__))

    # parse parameter
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help = "Path to the input fastq file, including the fastq file name", required = True, type = str)
    parser.add_argument("-o", "--output", help = "Output directory or file name", required = True, type = str)
    parser.add_argument("-t", "--threads", help = "Number of threads", default = 12, type = int)
    parser.add_argument("-q", "--qscore", help = "Read-level qscore filter", default = 10, type = int)
    parser.add_argument("-m", "--model", help = "Path to the output model", default = f"{current_dir}/model", type = str)
    parser.add_argument("-v", "--version", help = "Print software version info", action = "store_true")
    
    # print version info
    if "--version" in sys.argv[1 : ] or "-v" in sys.argv[1 : ]:
        print(f"Longbow version {'.'.join(version)} based on Python 3.7+")
        sys.exit(0)
    
    args = parser.parse_args()
    threads = args.threads
    fastqfile = str(args.input)
    output_name = str(args.output)
    assert output_name.lower().endswith('.json'), "-o output must end with .json file"

    model_path = str(args.model)
    qscore_cutoff = int(args.qscore)
    

    # print(fastqfile, threads, qscore_cutoff)
    subject = get_qscore(fastqfile, threads, qscore_cutoff)
    pred_dict = {"Sample" : fastqfile, "Flowcell" : None, "Software" : None, "Version" : None, "Mode" : None}
    if guppy_or_dorado(subject) == "guppy":
        pred_dict["Software"] = "guppy"
        train_x, train_y = read_train_file(f"{model_path}/train_guppy.csv")
        predict = predict_knn(subject, train_x, train_y)
        pred_dict["Flowcell"], pred_dict["Version"], pred_dict["Mode"] = decode(predict, "guppy")
    else:
        pred_dict["Software"] = "dorado"
        pred_dict["Version"] = '0'
        train_x, train_y = read_train_file(f"{model_path}/train_dorado.csv")
        predict = predict_knn(subject, train_x, train_y)
        pred_dict["Flowcell"], pred_dict["Mode"] = decode(predict, "dorado")

    # output to json
    if os.path.isfile(output_name):
        with open(output_name, 'w') as json_file:
            json.dump(pred_dict, json_file, indent = 4, separators=(',', ': '))
            json_file.write('\n')



if __name__ == "__main__":
    main()



