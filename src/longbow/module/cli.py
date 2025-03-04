import argparse
import os
import sys


def parse_args(version: tuple, current_dir: str):
    """
    Parse command-line arguments for the Longbow software.

    This function sets up and parses the command-line interface (CLI) arguments for Longbow.
    It handles options such as input/output file paths, thread count, quality score filtering,
    model paths, and more. Additionally, if --version is set, print the software version and exit.

    Args:
        version (tuple): A tuple representing the software version.
        current_dir (str): The current directory path, used to set default model paths.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.
    """
    
    # print version info
    if "--version" in sys.argv[1 : ] or "-v" in sys.argv[1 : ]:
        print(f"Longbow version {'.'.join(version)} based on Python 3.7+")
        sys.exit(0)

    # parse parameter
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help = "Path to the input fastq/fastq.gz file (required)", required = True, type = str)
    parser.add_argument("-o", "--output", help = "Path to the output json file [default : None]", default = None, type = str)
    parser.add_argument("-t", "--threads", help = "Number of parallel threads to use [default: 12]", default = 12, type = int)
    parser.add_argument("-q", "--qscore", help = "Minimum read QV to filter reads [default: 0]", default = 0, type = int)
    parser.add_argument("-m", "--model", help = "Path to the trained model [default: {longbow_code_base}/model]", default = os.path.join(current_dir, "model"), type = str)
    parser.add_argument("-a", "--ar", help = r"Enable autocorrelation for basecalling mode prediction HAC/SUP(hs) or FAST/HAC/SUP (fhs) (Options: hs, fhs, off) [default: fhs]", default = "fhs", type = str)
    parser.add_argument("-b", "--buf", help = "Output intermediate QV, autocorrelation results, and detailed run info to output json file", action = "store_true")
    parser.add_argument("-c", "--rc", help = "Enable read QV cutoff for mode correction in Guppy5/6 [default: on]", default = "on", type = str)
    parser.add_argument("--stdout", help = "Print results to standard output", action = "store_true")
    parser.add_argument("-v", "--version", help = "Print software version info and exit", action = "store_true")
    
    
    args = parser.parse_args()

    # sequential check and clean paramters
    args.input = os.path.abspath(args.input)

    if args.output:
        args.output = os.path.abspath(args.output)
    

    if args.ar not in ("off", "hs", "fhs"):
        raise ValueError(r"-a or --ar input error, must be off, hs, or fhs")
    if args.ar == "off":
        args.ar = False

    if args.buf:
        if not args.output:
            raise ValueError("-b or --buf must be set with -o and --output")
    
    if args.rc not in ("off", "on"):
        raise ValueError(r"-c or --rc input error, must be either off or on")
    if args.rc == "on":
        args.rc = True
    else:
        args.rc = False


    return args
