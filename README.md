# Longbow (Lucid dOrado aNd Guppy Basecalling cOnfig predictor)
##### Author : Mencius Jun @ Fudan University
##### Date : Jan. 18. 2023

## Briefing
Longbow is a python-based predictor for quality control of basecalling output of oxford nanopore sequencing.

It accept fastq file for input and predict :
1. the basecalling software was used (Dorado or Guppy);
2. the Nanopore flowcell version (R9 / R10);
3. the major guppy basecaller version(Dorado0, Guppy2, Guppy3/4, Guppy5/6);
4. the basecalling mode (FAST, HAC, SUP, NONE)


## Installation

Longbow can operate in any operation system with Python 3.7+ environment. The majority of the dependent package is offical built-in package.
Environment installation is simple.</br>

$ conda env create -f ont-longbow.yml

Then you are ready to go.

## Usage
usage: longbow.py [-h] -i INPUT -o OUTPUT [-t THREADS] [-q QSCORE] [-m MODEL] [-v] </br>

optional arguments: </br>
usage: longbow.py [-h] -i INPUT -o OUTPUT [-t THREADS] [-q QSCORE] [-m MODEL] [-c] [-V] [-v] </br>

optional arguments:
  -h, --help                        show this help message and exit  </br>
  -i INPUT, --input INPUT           Path to the input fastq file, including the fastq filename  </br>
  -o OUTPUT, --output OUTPUT        Output directory or file name  </br>
  -t THREADS, --threads THREADS     Number of parallel threads  </br>
  -q QSCORE, --qscore QSCORE        Read-level qscore filter  </br>
  -m MODEL, --model MODEL           Path to the training model csv data  </br>
  -c, --corr                        Do autocorrelation of hac/sup config or not  </br>
  -V, --verbose                     Verbose mode, print the result  </br>
  -v, --version                     Print software version info  </br>



## Features
1. Alignment-free / reference free.
2. In theory compatabible with ONT basecalled data for almost all species, and all kinds of DNA (cDNA, cfDNA, mtDNA).
3. Fast predictor thanks to CPU parallism and low RAM usage.
4. User-friendly graphical html output.

## Contributing author
PENDING

## Pulication

## License
PENDING
