# Longbow
Lucid dOrado aNd Guppy Basecalling cOnfig predictor


## Brief
Longbow is a python-based predictor for quality control of basecalling output of oxford nanopore sequencing.

It accept fastq file for input and predict :
1. the basecalling software was used (Dorado or Guppy);
2. the Nanopore flowcell version (R9 / R10);
3. the major guppy basecaller version(Dorado0, Guppy2, Guppy3/4, Guppy5/6);
4. the basecalling mode (FAST, HAC, SUP, NONE)


## Installation

Longbow can operate in any operation system with Python 3.7+ environment. The majority of the dependent package is offical built-in package.
Environment installation is simple.</br>

```sh
$ conda env create -f ont-longbow.yaml
```
Then you are ready to go.



## Usage
```sh
usage: longbow.py [-h] -i INPUT -o OUTPUT [-t THREADS] [-q QSCORE] [-m MODEL] [-c] [-V] [-v]

arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input fastq file, including the fastq file
                        name
  -o OUTPUT, --output OUTPUT
                        Output directory or file name
  -t THREADS, --threads THREADS
                        Number of parallel threads
  -q QSCORE, --qscore QSCORE
                        Read-level qscore filter
  -m MODEL, --model MODEL
                        Path to the training model csv data
  -c, --corr            Do autocorrelation of hac/sup config or not
  -V, --verbose         Verbose mode, print the result
  -v, --version         Print software version info
```



## Features
1. Alignment-free / reference free.
2. In theory compatabible with ONT basecalled data for almost all species, and all kinds of DNA (cDNA, cfDNA, mtDNA, cpDNA).
3. Fast predictor due to CPU parallism and low RAM usage.
4. User-friendly graphical website for quick search.


## License
PENDING
