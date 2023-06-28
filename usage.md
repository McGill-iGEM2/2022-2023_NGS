# Usage
## src directory
To implement the NGS pipeline, several bash scripts can be found in the `src` directory. Each script is numbered to indicate the order in which they should be run. The scripts are as follows:
```
.
├── 1-align
│   └── alignment.sh
├── 2-clean
│   └── clean.sh
├── 3-variant
│   └── variant.sh
└── 4-frequency
    └── empty.txt
```
This MD file documents the usage of each script. Details on the pipeline can be found in the [README.md](README.md) file.

## 1-align
The `alignment.sh` script is used to align the fastq files to the reference genome. The script takes in the following arguments:

```
Usage: ./alignment.sh [options] -r <reference genome> -p -n <normal sample 1> <normal sample 2> -t <tumor sample 1> <tumor sample 2>
Options:
    -r  Reference genome
    -p  Whether paired-end reads concatenated in one file
    -n  Normal sample(s)
    -t  Tumor sample(s)
```

The `-p` flag is used to indicate whether the paired-end reads are concatenated in one file. If the flag is not used, the script will assume that the paired-end reads are in two separate files. If only one file is passed and the `-p` flag is not used, the script will assume that the reads are single-end.

The `-n` and `-t` flags are used to indicate the normal and tumor samples, respectively. The script can take in multiple samples. The samples should be passed in the order in which they appear in the fastq file. For example, if the fastq file contains the following samples:

```
sample1_normal
sample1_tumor
sample2_normal
sample2_tumor
```

The script should be run as follows:

```
./alignment.sh -r reference.[fa|fasta|*.gz] -p -n sample1_normal.[fastq|*.gz] sample2_normal.[fastq|*.gz] -t sample1_tumor.[fastq|*.gz] sample2_tumor.[fastq|*.gz]
```

## 2-clean
The `clean.sh` script is used to clean the SAM files. The script takes in the following arguments:

```
Usage: ./clean.sh [options] -n <aligned normal sample> -t <aligned tumour sample>
Options:
    -n  Aligned normal sample
    -t  Aligned tumour sample
```

The `-n` and `-t` flags are used to indicate the normal and tumor samples, respectively. 
The script should be run as follows:

```
./clean.sh -n aligned_normal.sam -t aligned_tumor.sam
```

## 3-variant
The `variant.sh` script is used to call variants. This can be done either in tumour-only mode or in matched (tumour-normal) mode. The script takes in the following arguments:

```
Usage: ./variant.sh [options] -r <reference genome> -n <cleaned normal sample> -t <cleaned tumour sample> -o <output directory>
Options:
    -r  Reference genome
    -n  Cleaned normal sample
    -t  Cleaned tumour sample
    -o  Output directory
```

The `-n` and `-t` flags are used to indicate the normal and tumor samples, respectively. They should be called with the same reference genome that was used during alignment. The `-o` flag is used to indicate the output directory. By default, this script calls `Mutect2` in matched tumour-normal mode.

The script should be run as follows:

```
./variant.sh -r reference.[fa|fasta|*.gz] -n cleaned_normal.bam -t cleaned_tumor.bam -o output_directory
```