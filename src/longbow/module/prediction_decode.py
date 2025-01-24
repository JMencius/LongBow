def decode(code : int, software : str, aspect : str):
    """
    Decodes a numerical code into a human-readable configuration based on the software (guppy or dorado)
    and the aspect to be decoded (either 'qv' or 'mode').

    Parameters:
    code (int): The numerical code to decode.
    software (str): The basecalling software ('guppy' or 'dorado').
    aspect (str): The aspect to decode ('qv' for read quality values or 'hac_sup' for configurations).

    Returns:
    tuple or str: A tuple containing the decoded values for 'qv' aspect (e.g., ("R9", "3or4", "HAC"))
                  or a string for 'mode' (e.g., "HAC").

    Raises:
    KeyError: If the aspect or basecalling software is not recognized or the code does not exist in the mapping.
    """

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
        
        elif software == "dorado":
            dorado_map_code = {0 : ("R9", "FAST"),
                               1 : ("R9", "HAC"),
                               2 : ("R9", "SUP"),
                               3 : ("R10", "FAST"),
                               4 : ("R10", "HAC"),
                               5 : ("R10", "SUP")
                              }
            return dorado_map_code[code]
        else:
            raise KeyError("No specific basecalling software")

    elif aspect == "mode":
        hac_sup_map_code = {0 : "HAC",
                            1 : "SUP",
                            2 : "FAST"
                            }
        return hac_sup_map_code[code]

    else:
        raise KeyError("No specific decoding aspect")
