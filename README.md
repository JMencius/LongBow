# Longbow (Lucid dOrado aNd Guppy Basecalling cOnfig predictor)
##### Author : Mencius Jun @ Fudan University
##### Date : Nov. 17. 2023

## Briefing
Longbow is a python-based predictor for quality control of basecalling output of oxford nanopore sequencing.

It accept fastq file for input and predict :
1. the basecalling software was used (Dorado or Guppy);
2. the Nanopore flowcell version (R9 / R10);
3. the major guppy basecaller version(Guppy2, Guppy3/4, Guppy5/6);
4. the basecalling mode (FAST, HAC, SUP) under Guppy version


## Installation

Longbow can operate in any operation system with Python 3.7+ environment. The majority of the dependent package is offical built-in package.

Official built-in package used in Longbow:
collections, sys, multiprocessing, argparse, json, pickles, os, time

Third-party package used in Longbow:
numpy

Environment installation is simple.

$ conda create -n longbow python=3.7; <br>
$ conda activate longbow; <br>
$ conda install numpy; <br>

Then you are ready to go.

## Usage
usage: python longbow.py [-h] -i INPUT [-o OUTPUT] [-t THREADS] [-v] [-l LABEL] [-j JSON]

optional arguments:</br>
  -h, --help                     show this help message and exit </br>
  -i INPUT, --input INPUT        Path to the input fastq file, including the fastq file name </br>
  -o OUTPUT, --output OUTPUT     Output directory or file name, default output to current directory </br>
  -t THREADS, --threads THREADS  Number of threads </br>
  -v, --version                  Print software version info </br>
  -l LABEL, --label LABEL        Data label </br>
  -j JSON, --json JSON  Path to the training data json input; default : 0 means used pretrained model </br>



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
