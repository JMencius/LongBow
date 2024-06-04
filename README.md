# Longbow - Lucid dOrado aNd Guppy Basecalling cOnfig predictor

## Introduction
Longbow is a python-based predictor for quality control of basecalling output of oxford nanopore sequencing.

It accept only the `FASTQ` file for input and predict:
1. the basecalling software was used (Dorado or Guppy);
2. the Nanopore flowcell version (R9 / R10);
3. the major guppy basecaller version(Dorado0, Guppy2, Guppy3/4, Guppy5/6);
4. the basecalling mode (FAST, HAC, SUP, NONE)


## Features
- Fast speed and low RAM usage due to CPU parallism, alignment-free.
- Good compatability : in theory compatabible with ONT basecalled data for almost all species, and all kinds of DNA (cDNA, cfDNA, mtDNA, cpDNA, metaDNA, etc.).
- User-friendly graphical website (LongBowDB) for quick search.


## Installation
Longbow can operate in most modern operation system with Python 3.7+ environment. 
### Option 1. Bioconda
__PENDING__


### Option 2. Build an anaconda virtural environment
The `ont-longbow.yaml` file is also included in the release, which you can simply recreate the author's python environment.
```bash
$ conda env create -f ont-longbow.yaml;
```

### Option 3. Build the environment manually
For users using Microsoft Windows operating system, the installation ways provided above mihgt result in chaos, try to install it manually.
```
conda install numpy pandas statsmodels;
pip install pyfastx;
pip install dictances;
```


## Usage
Only two parameter is mandatory 
- `INPUT` which is the input `FASTQ` file
- `OUTPUT` which is the outout `JSON` file is


Other usage of `longbow` is listed in below. 
```
usage: longbow.py [-h] -i INPUT -o OUTPUT [-t THREADS] [-q QSCORE] [-m MODEL] [-c] [-V] [-v]

arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input fastq file, including the fastq filename
  -o OUTPUT, --output OUTPUT
                        Output directory or file name
  -t THREADS, --threads THREADS
                        Number of parallel threads
  -q QSCORE, --qscore QSCORE
                        Read-level qscore filter
  -m MODEL, --model MODEL
                        Path to the training model csv data
  -c, --corr            Do autocorrelation of hac/sup config or not
  -V, --verbose         Verbose mode, print the result to stdout
  -v, --version         Print software version info
```
You can add more `THREADS` to tackle with large `FASTQ` file

