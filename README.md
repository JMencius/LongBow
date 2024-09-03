# Longbow - Lucid dOrado aNd Guppy Basecalling cOnfig predictor

## Introduction
Longbow is a python-based predictor for quality control of basecalling output of oxford nanopore sequencing.

It accept only the `FASTQ` file for input and predict:
1. the basecalling software was used (Dorado or Guppy);
2. the Nanopore flowcell version (R9 / R10);
3. the major guppy basecaller version(Dorado0, Guppy2, Guppy3/4, Guppy5/6);
4. the basecalling mode (FAST, HAC, SUP, NONE)

## Installation
Longbow can operate in most modern operation system with Python 3.7+ environment. 
### Option 1. Build an anaconda virtural environment
The `ont-longbow.yaml` file is also included in the release, which you can recreate the author's python environment using the following command.
```bash
$ conda env create -f ont-longbow.yaml&&bash post_install.sh;
```

### Option 2. Build the environment manually
For users using Microsoft Windows operating system or having troubles with the following installion, try to install it manually.
```
conda create -n ont-longbow python=3.7;
conda activate ont-longbow;
conda install numpy statsmodels;
pip install pyfastx;
pip install dictances;
```


## Usage
Only one parameter is mandatory 
- `INPUT` which is the input `fastq`/`fastq.gz` file

Other usage of `longbow` is listed in below. 
```
usage: longbow.py [-h] -i INPUT [-o OUTPUT] [-t THREADS] [-q QSCORE]
                  [-m MODEL] [-a AR] [-b] [-V] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input fastq/fastq.gz file
  -o OUTPUT, --output OUTPUT
                        Output file [default : NONE
  -t THREADS, --threads THREADS
                        Number of parallel threads [default : 12]
  -q QSCORE, --qscore QSCORE
                        Read-level qscore filter [default : 0]
  -m MODEL, --model MODEL
                        Path to the training model [default = ./model]
  -a AR, --ar AR        Do read-qv based correction or autocorrelation for
                        hac/sup config [default : fhs]
  -b, --buf             Output intermediate results of QV and autocorrelation
  -V, --verbose         Verbose mode, print the result
  -v, --version         Print software version info
```

## Features
- Fast speed and low RAM usage due to CPU parallism, alignment-free.
- Good compatability : in theory compatabible with ONT basecalled data for almost all species, and all kinds of DNA (cDNA, cfDNA, mtDNA, cpDNA, metaDNA, etc.).
- User-friendly graphical website (LongBowDB) for quick search.
