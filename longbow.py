import time
import os
import sys
import argparse
import warnings
import json
import math
from src.faster_get_qscore import get_qscore
from src.distinguish_software import guppy_or_dorado
from src.read_train import read_qv_train_file
from src.read_train import read_autocorr_train_file
from src.bhattacharyya_knn import predict_knn
from src.euclidean_knn import predict_hac_sup
from src.readqv_cutoff import cutoff_qv
from src.prediction_decode import decode


def main():
    start_time = time.time()
    warnings.simplefilter(action = "ignore", category = FutureWarning)
    warnings.simplefilter(action = "ignore", category = RuntimeWarning)
    version = ('2', '1', '1')
    script_dir = os.path.dirname(os.path.realpath(__file__))
    current_dir = os.path.dirname(os.path.realpath(__file__))

    # parse parameter
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help = "Input fastq/fastq.gz file", required = True, type = str)
    parser.add_argument("-o", "--output", help = "Output file [default : NONE, print to Stdout]", default = "", type = str)
    parser.add_argument("-t", "--threads", help = "Number of parallel threads [default : 12]", default = 12, type = int)
    parser.add_argument("-q", "--qscore", help = "Read-level qscore filter [default : 0]", default = 0, type = int)
    parser.add_argument("-m", "--model", help = "Path to the training model [default = ./model]", default = os.path.join(current_dir, "model"), type = str)
    parser.add_argument("-a", "--ar", help = r"Autocorrelation mode, fhs: use autocorrelation to distinguish FAST/HAC/SUP [default : fhs] [default : fhs]", default = "fhs", type = str)
    parser.add_argument("-b", "--buf", help = "Output intermediate results of QV and autocorrelation", action = "store_true")
    parser.add_argument("-V", "--verbose", help = "Verbose mode, print the result to Stdout", action = "store_true")
    parser.add_argument("-v", "--version", help = "Print software version info", action = "store_true")

    
    # print version info
    if "--version" in sys.argv[1 : ] or "-v" in sys.argv[1 : ]:
        print(f"Longbow version {'.'.join(version)} based on Python 3.7+")
        sys.exit(0)
    
    args = parser.parse_args()
    threads = args.threads
    fastqfile = os.path.abspath(args.input)

    if args.output:
        output_name = os.path.abspath(args.output)
    else:
        output_name = False
    
    model_path = args.model
    qscore_cutoff = args.qscore

    autocorr = args.ar
    if autocorr not in ("off", "hs", "fhs"):
        raise ValueError(r"-a or --ar input error, must be off, hs, or fhs")
    if autocorr == "off":
        autocorr = False

    buf = args.buf
    verbose = args.verbose
    
    # input not empty
    if os.path.getsize(fastqfile) == 0:
        print("Input file is empty!")
        sys.exit(1)


    # Get property from fastq file
    if autocorr:
        baseqv, corr, readqv, outliner = get_qscore(fastqfile, threads, qscore_cutoff, args.ar)
    else:
        baseqv, readqv, outliner = get_qscore(fastqfile, threads, qscore_cutoff, args.ar)
    
    # if fail to calculate autocorrelation
    for i in corr:
        if math.isnan(i):
            print("Abnormal QV, fail to calculate QV autocorrelation, check QV in input FASTQ file.")
            sys.exit(0)
       

    # calculate readqv cutoff
    readqv_cutoff = cutoff_qv(readqv)
    
    # QV score distribution prediction
    pred_dict = {"Sample" : fastqfile, "Flowcell" : None, "Software" : None, "Version" : None, "Mode" : None}

    if guppy_or_dorado(baseqv) == "guppy":
        pred_dict["Software"] = "guppy"
        train_x, train_y = read_qv_train_file("guppy", readqv_cutoff, model_path)
        predict = predict_knn(baseqv, train_x, train_y)
        pred_dict["Flowcell"], pred_dict["Version"], pred_dict["Mode"] = decode(predict, "guppy", "qv")
    else:
        pred_dict["Software"] = "dorado"
        pred_dict["Version"] = '0'
        train_x, train_y = read_qv_train_file(f"dorado", readqv_cutoff, model_path)
        predict = predict_knn(baseqv, train_x, train_y)
        pred_dict["Flowcell"], pred_dict["Mode"] = decode(predict, "dorado", "qv")

    # autocorrelation
    autocorr_flag = 0
    if autocorr == "hs":
        if pred_dict["Software"] == "dorado" and pred_dict["Mode"] != "FAST":
            autocorr_flag = 1
            if pred_dict["Flowcell"] == "R9":
                model_file = "R9D0"
                preset_lag = 10
            else:
                model_file = "R10D0"
                preset_lag = 21

        elif pred_dict["Software"] == "guppy" and pred_dict["Version"] == '5or6' and pred_dict["Mode"] != "FAST":
            autocorr_flag = 1
            if pred_dict["Flowcell"] == "R9":
                model_file = "R9G6"
                preset_lag = 9
            else:
                model_file = "R10G6"
                preset_lag = 100

    elif autocorr == "fhs":
        if pred_dict["Software"] == "dorado":
            autocorr_flag = 1
            if pred_dict["Flowcell"] == "R9":
                model_file = "R9D0"
                preset_lag = 2
            else:
                model_file = "R10D0"
                preset_lag = 3

        elif pred_dict["Software"] == "guppy" and pred_dict["Version"] == '5or6':
            autocorr_flag = 1
            if pred_dict["Flowcell"] == "R9":
                model_file = "R9G6"
                preset_lag = 10
            else:
                model_file = "R10G6"
                preset_lag = 100

        elif pred_dict["Software"] == "guppy" and pred_dict["Flowcell"] == "R9" and pred_dict["Version"] == '3or4':
            autocorr_flag = 1
            model_file = "R9G4"
            preset_lag = 10


    # model detail classifcation
    if autocorr_flag:
        readqv_mode = None
        if pred_dict["Version"] in ('5or6', '0'):
            if readqv_cutoff == 7:
                readqv_mode = "FAST"
            if readqv_cutoff == 8:
                readqv_mode = "HAC"
            if readqv_cutoff == 9:
                readqv_mode = "SUP"

        if readqv_mode != None:
            pred_dict["Mode"] = readqv_mode
        else:
            train_x, train_y = read_autocorr_train_file(model_file, readqv_cutoff, model_path, autocorr)
            hac_sup_pred = predict_hac_sup(corr, train_x, train_y, preset_lag)
            pred_dict["Mode"] = decode(hac_sup_pred, pred_dict["Software"], "hac_sup")


    end_time = time.time()

    # if verbose mode or output not set, print the prediction result interactively
    if verbose or (not output_name):
        print(pred_dict)
    
    if outliner != 0:
        print(f"In total {outliner} base out of range of 1-90")

    # output to json
    if output_name:
        with open(output_name, 'w') as json_file:
            json.dump(pred_dict, json_file, indent = 4, separators=(',', ': '))
            json_file.write("\n\n")

            # output intermediate results
            if buf:
                json.dump({"Run info" : {"LongBow version" : f"{'.'.join(version)}",
                                         "Input" : os.path.basename(fastqfile),
                                         "Output" : os.path.basename(output_name),
                                         "Model" : model_path,
                                         "Threads" : threads,
                                         "Run time" : f"{end_time - start_time} s",
                                         "Read QV cutoff" : f"Q{readqv_cutoff + 1}",
                                         "Base QV outliner count" : f"{outliner}",
                                         "Autcorrelation" : args.ar,
                                         "Detail output" : bool(args.buf),
                                         "Verbose mode" : bool(verbose)}}, json_file, indent = 4, separators=(',', ': '))
                json_file.write("\n\n")
                json.dump({"baseqv" : {i : baseqv[i] for i in range(len(baseqv))}}, json_file)
                json_file.write("\n\n")
                json.dump({"readqv" : {i + 1 : readqv[i] for i in range(len(readqv))}}, json_file)
                json_file.write("\n\n")
                if autocorr:
                    json.dump({"autocorrelation" : {i + 1 : corr[i] for i in range(len(corr))}}, json_file)
                    json_file.write("\n")

if __name__ == "__main__":
    main()



