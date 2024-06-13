# Decode the prediction numbers to the corresponding human understandable config
def decode(code : int, software : str, aspect : str) -> tuple:
    if aspect == "qv":
        if software == "guppy":
            guppy_map_code = {0 : ("R9",  '2', "NONE"),
                              1 : ("R9",  "3or4", "FAST"),
                              2 : ("R9",  "3or4", "HAC"),
                              3 : ("R9",  "5or6", "FAST"),
                              4 : ("R9",  "5or6", "HAC"),
                              5 : ("R9",  "5or6", "SUP"),
                              6 : ("R10",  "5or6", "FAST"),
                              7 : ("R10",  "5or6", "HAC"),
                              8 : ("R10",  "5or6", "SUP")
                             }
            return guppy_map_code[code]
        
        else:
            dorado_map_code = {0 : ("R9", "FAST"),
                               1 : ("R9", "HAC"),
                               2 : ("R9", "SUP"),
                               3 : ("R10", "FAST"),
                               4 : ("R10", "HAC"),
                               5 : ("R10", "SUP")
                              }
            return dorado_map_code[code]

    elif aspect == "hac_sup":
        hac_sup_map_code = {0 : "HAC",
                            1 : "SUP",
                            2 : "FAST"
                            }
        return hac_sup_map_code[code]

    else:
        raise KeyError("No specific decoding aspect")
