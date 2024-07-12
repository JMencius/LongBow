
import os

guppy_label = {"R9G2" : 0, 
               "R9G4FAST" : 1,
               "R9G4HAC" : 2,
               "R9G6FAST" : 3,
               "R9G6HAC" : 4,
               "R9G6SUP" : 5,
               "R10G6FAST" : 6,
               "R10G6HAC" : 7,
               "R10G6SUP" : 8}

dorado_label = {"R9D0FAST" : 0,
                "R9D0HAC" : 1,
                "R9D0SUP" : 2,
                "R10D0FAST" : 3,
                "R10D0HAC" : 4,
                "R10D0SUP" : 5}


def split_filename(filename : str) -> str:
    # e.g.  human_R10G6SUP.csv -> R10G6SUP
    temp = filename.split('.')[0]
    extract = temp.split('_')[-1]

    return extract



def csv2qvlist(filename : str, readqv_cutoff : int) -> list:
    qv = [0 for i in range(90)]
    linecount = 0
    with open(filename, 'r') as f:
        for line in f:
            if linecount == 0:
                linecount += 1
                continue
            m = line.split(',')
            if len(m) != 192:
                continue
            if int(m[0]) > readqv_cutoff:
                if int(m[1]) != 0:
                    for idx in range(90):
                        qv[idx] += int(m[idx + 2])

    return qv
                

def csv2autocorr(filename : str, readqv_cutoff : int) -> list:
    autocorr = [[0, 0] for i in range(100)]
    linecount = 0
    with open(filename, 'r') as f:
        for line in f:
            if linecount == 0:
                linecount += 1
                continue
            m = line.split(',')
            if len(m) != 192:
                continue
            if int(m[0]) > readqv_cutoff:
                if int(m[1]) != 0:
                    for i in range(100):
                        temp = (m[i + 92]).split('|')
                        if temp[0] != "nan":
                            autocorr[i][0] += float(temp[0])
                            autocorr[i][1] += int(temp[1])

    autocorr_output = [i[0] / i[1] for i in autocorr]

    return autocorr_output



def read_qv_train_file(software : str, readqv_cutoff : int, model_path : str) -> tuple:
    if software == "guppy":
        target = guppy_label
    else:
        target = dorado_label
    
    train_X = list()
    train_Y = list()
    for csvfile in os.listdir(model_path):
        csvtag = split_filename(csvfile)
        if csvtag in target:
            readX = csv2qvlist(os.path.join(model_path, csvfile), readqv_cutoff)
            if sum(readX) != 0:
                normalized_readX = [i / sum(readX) for i in readX]
                train_X.append(normalized_readX)
                train_Y.append(target[csvtag])
    return (train_X, train_Y)



def read_autocorr_train_file(target_tag : str, readqv_cutoff : int, model_path : str, autocorr : str) -> tuple:
    train_X = list()
    train_Y = list()
    for csvfile in os.listdir(model_path):
        if target_tag in csvfile:
            if autocorr == "fhs":
                train_X.append(csv2autocorr(os.path.join(model_path, csvfile), readqv_cutoff))
                if "HAC" in csvfile:
                    train_Y.append(0)
                if "SUP" in csvfile:
                    train_Y.append(1)
                if "FAST" in csvfile:
                    train_Y.append(2)
            if autocorr == "hs":
                if "HAC" in csvfile:
                    train_X.append(csv2autocorr(os.path.join(model_path, csvfile), readqv_cutoff))
                    train_Y.append(0)
                if "SUP" in csvfile:
                    train_X.append(csv2autocorr(os.path.join(model_path, csvfile), readqv_cutoff))
                    train_Y.append(1)

    return (train_X, train_Y)
    


























