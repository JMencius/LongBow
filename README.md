# Longbow - Lucid dOrado aNd Guppy Basecalling cOnfig predictor

## Introduction
Longbow is a Python-based tool designed for quality control of basecalling output from Oxford Nanopore sequencing.

It accepts a `FASTQ` file as input and predicts:
1. The basecalling software used (Dorado or Guppy);
2. The Nanopore flowcell version (R9 / R10);
3. The major basecaller version (Guppy2, Guppy3/4, Guppy5/6, Dorado0);
4. The basecalling mode (FAST, HAC, SUP).

## Installation
Longbow is compatible with most Linux operating systems and requires a Python 3.7+ environment.

### Option 1. Install LongBow via Bioconda
```bash
conda create -n longbow python=3.7;
conda install -c bioconda longbow;
```

### Option 2. Install LongBow through pip
**Due to name conflicts, online pip installation is not supported.**
First, download and unzip LongBow release, then navigate to the source code root directory containing setup.py.
```bash
conda create -n longbow python=3.7;

## Download and unzip LongBow, enter the source code root directory containing setup.py

pip install .;
```


## Usage
Only one parameter is mandatory:
- `-i or --input` which is the input `fastq`/`fastq.gz` file


Full parameters of `longbow` is listed in below. 
```
usage: longbow [-h] -i INPUT [-o OUTPUT] [-t THREADS] [-q QSCORE] [-m MODEL]
               [-a AR] [-b] [-c RC] [--stdout] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input fastq/fastq.gz file (required)
  -o OUTPUT, --output OUTPUT
                        Path to the output json file [default: None]
  -t THREADS, --threads THREADS
                        Number of parallel threads to use [default: 12]
  -q QSCORE, --qscore QSCORE
                        Minimum read QV to filter reads [default: 0]
  -m MODEL, --model MODEL
                        Path to the trained model [default:
                        {longbow_code_base}/model]
  -a AR, --ar AR        Enable autocorrelation for basecalling mode prediction
                        HAC/SUP(hs) or FAST/HAC/SUP (fhs) (Options: hs, fhs,
                        off) [default: fhs]
  -b, --buf             Output intermediate QV, autocorrelation results,
                        confidence score (experimental) and detailed run info
                        to output json file
  -c RC, --rc RC        Enable read QV cutoff for mode correction in Guppy5/6
                        [default: on]
  --stdout              Print results to standard output
  -v, --version         Print software version info and exit
```


## Examples
1. (Standard) Predict basecalling configuration of `reads.fastq.gz` and save the results to `pred.json`.
```
longbow -i reads.fastq.gz -o pred.json;
```

2. (Minimal input) Only input file is specified, and the prediction results are printed to standard output.
```
longbow -i reads.fastq;
```

3. (Dual output) Save results to `pred.json` and simultaneously print them to standard output.
```
longbow -i reads.fastq -o pred.json --stdout; 
```

4. (Detailed output) Save intermediate QV, autocorrelation results, confidence score (experimental), along with detailed parameters, to `pred.json`.
```
longbow -i reads.fastq -o pred.json -b;
```



## Resource consumption
Longbow can process 10,000 reads of ONT sequencing within seconds using 32 threads on modern Desktop CPU or Server CPU. 

In our tests with a large dataset (10<sup>7</sup> reads, approximately 100 GB in uncompressed format), LongBow completed processing within one hour using 32 threads.

The actual performance may vary depending on factors such as I/O speed, memory speed, and CPU capabilities.
