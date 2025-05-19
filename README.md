# Longbow - Lucid dOrado aNd Guppy Basecalling cOnfig predictor

## Introduction
Longbow is a Python-based tool designed for quality control of basecalling output from Oxford Nanopore sequencing.

It accepts a `FASTQ` file as input and predicts:
1. The basecalling software used (Dorado, Guppy);
2. The Nanopore flowcell version (R9, R10);
3. The major basecaller version (Guppy2, Guppy3/4, Guppy5/6, Dorado0);
4. The basecalling mode (FAST, HAC, SUP).
5. (In development) The basecalling model version of Dorado (pre-V4.2.0, V4.2.0/V4.3.0, V5.0.0).

## Installation
Longbow is compatible with most Linux operating systems and requires a Python 3.7+ environment.

### Option 1. Install LongBow through pip [recommended]
**Due to name conflict on PyPi, we have to use epg-longbow for pip installation.** It is just a name change
```bash
conda create -n longbow python=3.7;

pip install epg-longbow;
```


### Option 2. Install LongBow via Bioconda
```bash
conda create -n longbow python=3.7;
conda install -c bioconda longbow;
```


### Option 3. Install LongBow through Docker
The Docker Hub link for LongBow is in <https://hub.docker.com/r/jmencius/longbow/tags>.
```bash
docker pull jmencius/longbow:latest;
```

### Option 4. Local installation
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

## Test data
A small FASTQ file is provided in [here](./tests/data). You can use it as input for testing


## Resource consumption
Longbow can process 10,000 reads of ONT sequencing within seconds using 32 threads on modern Desktop CPU or Server CPU. 

In our tests with a large dataset (10<sup>7</sup> reads, approximately 100 GB in uncompressed format), LongBow completed processing within one hour using 32 threads.

The actual performance may vary depending on factors such as I/O speed, memory speed, and CPU capabilities.

## Known limitation
1. Basecalling results from early versions of Guppy 3 (e.g., Guppy 3.0.3) may be classified as Guppy 2 by LongBow. However, the impact of this misclassification is minor.


## Citation
Mencius, J., Chen, W., Zheng, Y. et al. Restoring flowcell type and basecaller configuration from FASTQ files of nanopore sequencing data. Nat Commun 16, 4102 (2025). 

<https://doi.org/10.1038/s41467-025-59378-x>
```
@article{mencius_restoring_2025,
	title = {Restoring flowcell type and basecaller configuration from {FASTQ} files of nanopore sequencing data},
	volume = {16},
	issn = {2041-1723},
	url = {https://doi.org/10.1038/s41467-025-59378-x},
	doi = {10.1038/s41467-025-59378-x},
	abstract = {As nanopore sequencing has been widely adopted, data accumulation has surged, resulting in over 700,000 public datasets. While these data hold immense potential for advancing genomic research, their utility is compromised by the absence of flowcell type and basecaller configuration in about 85\% of the data and associated publications. These parameters are essential for many analysis algorithms, and their misapplication can lead to significant drops in performance. To address this issue, we present LongBow, designed to infer flowcell type and basecaller configuration directly from the base quality value patterns of FASTQ files. LongBow has been tested on 66 in-house basecalled FAST5/POD5 datasets and 1989 public FASTQ datasets, achieving accuracies of 95.33\% and 91.45\%, respectively. We demonstrate its utility by reanalyzing nanopore sequencing data from the COVID-19 Genomics UK (COG-UK) project. The results show that LongBow is essential for reproducing reported genomic variants and, through a LongBow-based analysis pipeline, we discovered substantially more functionally important variants while improving accuracy in lineage assignment. Overall, LongBow is poised to play a critical role in maximizing the utility of public nanopore sequencing data, while significantly enhancing the reproducibility of related research.},
	number = {1},
	journal = {Nature Communications},
	author = {Mencius, Jun and Chen, Wenjun and Zheng, Youqi and An, Tingyi and Yu, Yongguo and Sun, Kun and Feng, Huijuan and Feng, Zhixing},
	month = may,
	year = {2025},
	pages = {4102},
}
```

