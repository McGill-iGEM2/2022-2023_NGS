# McGill iGEM NGS Bioinformatics Pipeline
A next-generation sequencing bioinformatics pipeline developed for the discovery and detection of somatic SNPs and indels in cancer patients. This repository consists of individual `bash` files that can be run locally as well as a series of `SLURM` scripts for execution on the cloud.

This repository outlines the exploratory process undertaken by 2022-2023 McGill iGEM's Dry Lab team. Here, we document the development of an NGS bioinformatics pipeline from start to finish, and explore its utility as an adjacent tool for the <i>PROTEUS</i> system presented at the 2023 iGEM Jamboree in Paris, France.

## General Overview

<i>PROTEUS</i> is a CRISPR-mediated gene-editing tool designed as a oncological precision medicine tool. The heterogeneous nature of cancer highlights the need of patient-specific guide RNAs (gRNAs) in targeted therapies to increase efficacy and minimize off-target effects. 

This pipeline is designed to take in whole-exome sequencing data of a patient's healthy and diseased tissue to produce genomic fragments (reads). These raw sequences will undergo several processes including (1) read mapping, (2) post-processing quality control, (3) variant calling, and (4) variant annotation and filtering. These fragments are aligned against a human reference genome for comparison. Deviations present in the tumour sample and absent in the normal sample are flagged as candidate variants, which can then be cross-referenced with existing variation databases to gauge the degree of pathogenicity.

## Repository Structure
```
.
├── README.md
├── data
│   ├── normal
│   │   ├── bqsr
│   │   └── raw
│   ├── reference
│   ├── tumor_sample_1
│   │   ├── bqsr
│   │   └── raw
│   ├── tumor_sample_2
│   │   ├── bqsr
│   │   └── raw
│   └── variant
│   │   ├── bqsr
│   │   └── pcgr
├── scripts
│   ├── alignment.sh
│   ├── annotate.sh
│   ├── clean.sh
│   ├── variant.sh
└── src
    ├── 1_align.sh
    ├── 2_clean.sh
    ├── 3_variant.sh
    ├── 4_annotate.sh
    └── batch.sh
  ```

## Local Installation
Many of the software used are available through `bioconda`, a distribution of bioinformatics software as part of the Conda package manager. Conda simplifies the steps required to run a program, and set-up is often as simple as running one line of code.

- miniconda https://docs.conda.io/projects/miniconda/en/latest/
- sratoolkit https://anaconda.org/bioconda/sra-tools
- bwa https://anaconda.org/bioconda/bwa
- samtools https://anaconda.org/bioconda/samtools
- gatk https://anaconda.org/bioconda/gatk4
- pcgr https://anaconda.org/pcgr/pcgr 

## Cloud Cluster

It is recommended to run the pipeline on the cloud, where memory and compute resources are more extensive than local machines. The scripts located in the `src` directory are executed through the `SLURM` workload manager, where the configurations can be adjusted in the header lines of each file. We have also included a `batch.sh` script that can be used to run the pipeline on multiple samples at once.

Each script will load the necessary modules required to run. For annotation, we used PCGR, which is available through the Alliance's cvmfs distribution. This can be quickly setup by running the `setup.sh` script in the main directory, which will add the necessary lines to your `~/.bashrc` file. This will allow access to the MUGQIC software repository, which contains a suite of bioinformatics tools including PCGR. Running this will also install the data bundle in the `PCGR/`directory that is referenced by PCGR during annotation. By default, it will install the bundle for the GRCh38 build, but this can be changed by editing the `setup.sh` file.

## Usage
- using the batch script
- more detailed explanation on each of the src files here (link)