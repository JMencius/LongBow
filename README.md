# NAGP
##### Author : Mencius Jun @ Fudan University
##### Date : Sep/12/2023

## Briefing

NAGP stands for Nanopore And Guppy version Predictor, does not stands for Not Available Good Price :).

It accept fastq file for input and predict the Nanopore version (R9 / R10) and the guppy basecaller version(Guppy2, Guppy3/4, Guppy5/6) use to generate the fastq files.

The output result in shorten in format such as R10G6 (Nanopore version R9, basecalled with Guppy6).


## Pipeline
1. [optional] Random downsampling to 1000 or use-defined sequences.
2. Minimap2 alignment to reference genome.
3. Samtools filter.
4. Python script calculates error rate in samtools-filtered sam files.
5. Python script calculates max Phred score in samtools-filtered sam files.
6. Predict based on error rate and Phred score.


## Installation

NAGP utilizes minimap2 for sequence alignment to the reference genome and samtools for filtering the alignment files. 
Subsequent python scripts utilize multiprocessing to accelerate sam files processing.

Environment installation is simple.

$ conda create -n nagp python=3.7; <br>
$ conda activate nagp; <br>
$ conda install -c bioconda minimap2==2.17 samtools==1.6;   # Other minimap or samtools version should work as well. <br>

Then you are ready to go.

## Usage

Usage: ./nagp -i {fastq input path} -r {reference genome input path} [options]

Options: <br>
Mandatory Options: <br>
-r, --ref      Path to reference genome in fasta format <br>
-i, --input    Path to the input combined fastq file <br>

Other Options: <br>
-n, --numbers    Numbers of sequences to be sampled [DEFAULT 1000] <br>
-a, --all        Use all the sequences in the input files <br>
-t, --threads    How many CPU threads to use [DEFAULT 12] <br>
-h, --help       Display this help and exit <br>
-v, --version    Print nagp version info <br>


## Features

1. Fast prediction, due to harnessment of multicores. For an 11G Escherichia coli fastq file, nagp gives predictions in approximately 1 minute using typical 12 CPU cores.
2. Support random downsampling fastq sequences. The random downsampling is based on the Fisher-Yates shuffle algorithm using the random.shuffle() in python. Or for small fastq files, user can use -a/--all to use all sequences for prediction.



## Contributing
PENDING
## License
PENDING
