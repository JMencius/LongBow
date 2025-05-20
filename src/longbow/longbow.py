import time
import os
import sys
import argparse
import warnings
import json
import math

try:
    from longbow.module.faster_get_qscore import get_qscore
    module_path = 'longbow.module'
except ImportError:
    # Fallback if longbow is not available
    module_path = 'module'

# Now use module_path dynamically for subsequent imports
parse_args = __import__(f"{module_path}.cli", fromlist=["parse_args"]).parse_args
get_qscore = __import__(f"{module_path}.faster_get_qscore", fromlist=["get_qscore"]).get_qscore
guppy_or_dorado = __import__(f"{module_path}.distinguish_software", fromlist=["guppy_or_dorado"]).guppy_or_dorado
read_qv_train_file = __import__(f"{module_path}.read_train", fromlist=["read_qv_train_file"]).read_qv_train_file
read_autocorr_train_file = __import__(f"{module_path}.read_train", fromlist=["read_autocorr_train_file"]).read_autocorr_train_file
predict_knn = __import__(f"{module_path}.bhattacharyya_knn", fromlist=["predict_knn"]).predict_knn
predict_mode = __import__(f"{module_path}.euclidean_knn", fromlist=["predict_mode"]).predict_mode
cutoff_qv = __import__(f"{module_path}.readqv_cutoff", fromlist=["cutoff_qv"]).cutoff_qv
decode = __import__(f"{module_path}.prediction_decode", fromlist=["decode"]).decode


def main():
    start_time = time.time()
    warnings.simplefilter(action = "ignore", category = FutureWarning)
    warnings.simplefilter(action = "ignore", category = RuntimeWarning)
    version = ('2', '3', '1')
    script_dir = os.path.dirname(os.path.realpath(__file__))
    current_dir = os.path.dirname(os.path.realpath(__file__))

    args = parse_args(version, script_dir)    
    
    # transfer parameters to variables
    fastqfile = args.input
    output_name = args.output
    threads = args.threads
    qscore_cutoff = args.qscore
    model_path = args.model
    autocorr = args.ar
    buf = args.buf
    readqvcorrect = args.rc
    stdout = args.stdout
    
    # Check if input is empty
    if os.path.getsize(fastqfile) == 0:
        print("Input file is empty!")
        sys.exit(1)


    # Get property from fastq file
    if autocorr:
        baseqv, corr, readqv, outliner = get_qscore(fastqfile, threads, qscore_cutoff, autocorr)
    else:
        baseqv, readqv, outliner = get_qscore(fastqfile, threads, qscore_cutoff, autocorr)
    
    # Check if autocorrelation calculation is normal
    if autocorr:
        for i in corr:
            if math.isnan(i):
                print("Abnormal QV, fail to calculate QV autocorrelation, check QV in input FASTQ file.")
                sys.exit(0)
       

    # calculate readqv cutoff
    readqv_cutoff = cutoff_qv(readqv)
    
    # QV score distribution prediction
    pred_dict = {"Sample" : os.path.basename(fastqfile), "Flowcell" : None, "Software" : None, "Version" : None, "Mode" : None, "Confidence level": None}
    
    l1_confidence, l2_confidence = None, None
    if guppy_or_dorado(baseqv) == "guppy":
        pred_dict["Software"] = "guppy"
        train_x, train_y = read_qv_train_file("guppy", readqv_cutoff, model_path)
        predict, l1_confidence  = predict_knn(baseqv, train_x, train_y, "guppy")
        pred_dict["Flowcell"], pred_dict["Version"], pred_dict["Mode"] = decode(predict, "guppy", "qv")
    else:
        pred_dict["Software"] = "dorado"
        pred_dict["Version"] = '0'
        train_x, train_y = read_qv_train_file(f"dorado", readqv_cutoff, model_path)
        predict, l1_confidence = predict_knn(baseqv, train_x, train_y, "dorado")
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
                preset_lag = 10
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
        if readqvcorrect:
            if pred_dict["Version"] in ('5or6'):
                if readqv_cutoff == 7:
                    readqv_mode = "FAST"
                    l2_confidence = 1
                if readqv_cutoff == 8:
                    readqv_mode = "HAC"
                    l2_confidence = 1
                if readqv_cutoff == 9:
                    readqv_mode = "SUP"
                    l2_confidence = 1
        
        if readqv_mode != None:
            pred_dict["Mode"] = readqv_mode
        else:
            train_x, train_y = read_autocorr_train_file(model_file, readqv_cutoff, model_path, autocorr)
            pred_mode, l2_confidence = predict_mode(corr, train_x, train_y, preset_lag)
            pred_dict["Mode"] = decode(pred_mode, pred_dict["Software"], "mode")

    # Confidence level 
    overall_confidence = 0
    if l2_confidence:
        overall_confidence = (l1_confidence + l2_confidence) / 2
    else:
        overall_confidence = l1_confidence

    if overall_confidence < 0.2:
        pred_dict["Confidence level"] = "very low"
    elif overall_confidence < 0.4:
        pred_dict["Confidence level"] = "low"
    elif overall_confidence < 0.6:
        pred_dict["Confidence level"] = "medium"
    elif overall_confidence < 0.8:
        pred_dict["Confidence level"] = "high"
    else:
        pred_dict["Confidence level"] = "very high"

    end_time = time.time()

    # if verbose mode or output not set, print the prediction result interactively
    if stdout or (not output_name):
        print(pred_dict)
    
    if outliner != 0:
        print(f"WARNING in total {outliner} base out of range of 1-90")

    # output to json
    if output_name:
        with open(output_name, 'w') as json_file:
            if buf:
                conf_dict = {"Overall confidence score": overall_confidence}
                pred_dict.update(conf_dict)
                buf_dict = {"Run info" : {"LongBow version" : f"{'.'.join(version)}",
                                         "Input" : os.path.basename(fastqfile),
                                         "Output" : os.path.basename(output_name),
                                         "Model" : model_path,
                                         "Threads" : threads,
                                         "Run time" : f"{end_time - start_time} s",
                                         "Read QV cutoff" : f"Q{readqv_cutoff + 1}",
                                         "Read QV for mode correction" : bool(readqvcorrect),
                                         "Base QV outliner count" : f"{outliner}",
                                         "Autcorrelation" : args.ar,
                                         "Detail output" : bool(args.buf),
                                         "Stdout" : bool(stdout)}}
                pred_dict.update(buf_dict)
                baseqv_dict = {"baseqv" : {i : baseqv[i] for i in range(len(baseqv))}}
                pred_dict.update(baseqv_dict)
                if autocorr:
                    autodict = {"autocorrelation" : {i + 1 : corr[i] for i in range(len(corr))}}
                pred_dict.update(autodict)

            json.dump(pred_dict, json_file, indent = 4, separators=(',', ': '))
        


if __name__ == "__main__":
    main()



