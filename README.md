# Longbow (L<br><u>_</u>ucid d\_o\_rado a\_n\_d \_g\_uppy \_b\_asecalling c\_o\_nfig predictor)
##### Author : Mencius Jun @ Fudan University
##### Date : Nov/3/2023

## Briefing
Longbow is a python-based predictor for quality control of basecalling output of oxford nanopore sequencing.

It accept fastq file for input and predict :
1. the basecalling software was used (Dorado or Guppy);
2. the Nanopore flowcell version (R9 / R10);
3. the major guppy basecaller version(Guppy2, Guppy3/4, Guppy5/6);
4. the basecalling mode (FAST, HAC, SUP)


## Installation

Longbow can operate in any environment with Python 3.7+. The majority of the dependent package is offical built-in package.

Official built-in package used in Longbow:
collections, sys, multiprocessing, argparse, json
Third-party package used in Longbow:
Numpy

Environment installation is simple.

$ conda create -n longbow python=3.7; <br>
$ conda activate nagp; <br>
$ conda install numpy; <br>

Then you are ready to go.

## Usage

Usage: python longbow.py -i {fastq input path} -r {reference genome input path} [options]

Options: <br>
Mandatory Options: <br>
-i, --input    Path to the input combined fastq file <br>

Other Options: <br>
-t, --threads    How many CPU threads to use [DEFAULT 12] <br>
-h, --help       Display this help and exit <br>
-v, --version    Print longbow version info <br>


## Features
1. Alignment-free / reference free.
2. In theory compatabible with ONT basecalled data for almost all species, and all kinds of DNA (cDNA, cfDNA, mtDNA).
3. Fast predictor thanks to CPU parallism.
4. User-friendly graphical html output.

## Contributing author
PENDING
## License
PENDING
