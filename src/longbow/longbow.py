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
faster_get_qscore = __import__(f"{module_path}.faster_get_qscore", fromlist=["get_qscore"]).get_qscore
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
    version = ('2', '3', '0')
    script_dir = os.path.dirname(os.path.realpath(__file__))
    current_dir = os.path.dirname(os.path.realpath(__file__))

    # parse parameter
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help = "Input fastq/fastq.gz file", required = True, type = str)
    parser.add_argument("-o", "--output", help = "Output file [default : None]", default = None, type = str)
    parser.add_argument("-t", "--threads", help = "Number of parallel threads [default : 12]", default = 12, type = int)
    parser.add_argument("-q", "--qscore", help = "Read-level qscore filter [default : 0]", default = 0, type = int)
    parser.add_argument("-m", "--model", help = "Path to the training model [default = ./model]", default = os.path.join(current_dir, "model"), type = str)
    parser.add_argument("-a", "--ar", help = r"Do read-qv based correction or autocorrelation for hac/sup config [default : fhs]", default = "fhs", type = str)
    parser.add_argument("-b", "--buf", help = "Output intermediate results of QV and autocorrelation", action = "store_true")
    parser.add_argument("-c", "--rc", help = "Use Read QV cutoff to conduct mode correction for R9G5/6 or not [default = on]", default = "on", type = str)
    parser.add_argument("--stdout", help = "Verbose mode, print the result to StdOut", action = "store_true")
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
    stdout = args.stdout
    readqvcorrect = args.rc
    if readqvcorrect not in ("off", "on"):
        raise ValueError(r"-c or --rc input error, must be either off or on")
    if readqvcorrect == "on":
        mc = True
    else:
        mc = False
    
    # input not empty
    if os.path.getsize(fastqfile) == 0:
        print("Input file is empty!")
        sys.exit(1)


    # Get property from fastq file
    if autocorr:
        baseqv, corr, readqv, outliner = get_qscore(fastqfile, threads, qscore_cutoff, autocorr)
    else:
        baseqv, readqv, outliner = get_qscore(fastqfile, threads, qscore_cutoff, autocorr)
    
    # if fail to calculate autocorrelation
    ## print(corr)
    if autocorr:
        for i in corr:
            if math.isnan(i):
                print("Abnormal QV, fail to calculate QV autocorrelation, check QV in input FASTQ file.")
                sys.exit(0)
       

    # calculate readqv cutoff
    readqv_cutoff = cutoff_qv(readqv)
    
    # QV score distribution prediction
    pred_dict = {"Sample" : os.path.basename(fastqfile), "Flowcell" : None, "Software" : None, "Version" : None, "Mode" : None}
    
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
        if mc:
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


    end_time = time.time()

    # if verbose mode or output not set, print the prediction result interactively
    if stdout or (not output_name):
        print(pred_dict)
    
    if outliner != 0:
        print(f"In total {outliner} base out of range of 1-90")

    # output to json
    if output_name:
        with open(output_name, 'w') as json_file:
            confidence_dict = {"Flowcell and basecaller confidence": l1_confidence, "Basecalling mode confidence": l2_confidence}
            pred_dict.update(confidence_dict)
            if buf:
                buf_dict = {"Run info" : {"LongBow version" : f"{'.'.join(version)}",
                                         "Input" : os.path.basename(fastqfile),
                                         "Output" : os.path.basename(output_name),
                                         "Model" : model_path,
                                         "Threads" : threads,
                                         "Run time" : f"{end_time - start_time} s",
                                         "Read QV cutoff" : f"Q{readqv_cutoff + 1}",
                                         "Read QV for mode correction" : bool(mc),
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



