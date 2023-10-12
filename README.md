# McGill iGEM NGS Bioinformatics Pipeline
A next-generation sequencing bioinformatics pipeline developed for the discovery and detection of somatic SNPs and indels in cancer patients. This repository consists of individual `bash` files that can be run locally as well as a series of `SLURM` scripts for execution on the cloud.

This repository outlines the exploratory process undertaken by 2022-2023 McGill iGEM's Dry Lab team. Here, we document the development of an NGS bioinformatics pipeline from start to finish, and explore its utility as an adjacent tool for the <i>PROTEUS</i> system presented at the 2023 iGEM Jamboree in Paris, France.

## Description

<i>PROTEUS</i> is a CRISPR-mediated gene-editing tool designed as a oncological precision medicine tool. The heterogeneous nature of cancer highlights the need of patient-specific guide RNAs (gRNAs) in targeted therapies to increase efficacy and minimize off-target effects. 

This pipeline is designed to take in whole-exome sequencing data of a patient's healthy and diseased tissue to produce genomic fragments (reads). These raw sequences will undergo several processes including (1) read mapping, (2) post-processing quality control, (3) variant calling, and (4) variant annotation and filtering. These fragments are aligned against a human reference genome for comparison. Deviations present in the tumour sample and absent in the normal sample are flagged as candidate variants, which can then be cross-referenced with existing variation databases to gauge the degree of pathogenicity.

Included in this repository are all of the scripts we developed throughout the year. The `data` directory contains the raw sequencing data, as well as the reference genome used for alignment. The `scripts` directory contains command line calls to each software we used, and is referenced by the `src` directory to be run on the cloud. For obvious storage reasons, the data we used is not included, but we've included links to the data we used [here](4.-References.md).

## Repository Structure
```
.
├── README.md
├── data
│   ├── normal
│   ├── reference
│   ├── tumor_sample_1
│   │   ├── tumor
│   │   └── variant
│   │       └── pcgr
│   └── tumor_sample_2
│       ├── tumor
│       └── variant
│           └── pcgr
├── PCGR
├── scripts
│   ├── alignment.sh
│   ├── annotate.sh
│   ├── clean.sh
│   └── variant.sh
├── src
│   ├── 1_align.sh
│   ├── 2_clean.sh
│   ├── 3_variant.sh
│   ├── 4_annotate.sh
│   └── batch.sh
├── README.md
└── setup.sh
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
To start using the pipeline, clone the repository to your local machine. The following is a brief overview of each directory:

- `data` - contains the raw sequencing data, as well as the reference genome used for alignment
- `PCGR` - contains the data bundle used by PCGR for annotation
- `scripts` - contains command line calls to each software we used, and is referenced by the `src` directory to be run on the cloud
- `src` - contains the `SLURM` scripts that are executed on the cloud

This pipeline can be run using the scripts in the `scripts` directory. The following is a brief overview of each script:

- `alignment.sh` - aligns the reads to the reference genome
- `clean.sh` - post-processing quality control
- `variant.sh` - variant calling
- `annotate.sh` - variant annotation and filtering

If running the pipeline locally, the scripts should be run in this order. A more detailed description of each script can be found in the [Usage](2.-Usage.md) page.

If running the pipeline on the cloud, the scripts can be called using the `src` directory, which are numbered in order of execution.

- `1_align.sh`
- `2_clean.sh`
- `3_variant.sh`
- `4_annotate.sh`

Inside each script is a header section that informs SLURM of the resources required to run the script. This includes the number of CPUs, memory, and time required to run the script. Some default values are provided, but these can be adjusted to suit your needs. The scripts can be run by calling `sbatch` on the script, e.g. `sbatch 1_align.sh` followed by the necessary arguments.

Further down are is the command line call to the corresponding script in the `scripts` directory. These can be adjusted as needed (see [Usage](2.-Usage.md) for more details on optional flags).

To run multiple samples at once, the `batch.sh` script can be used. This script will call the numbered scripts in the `src` directory, and can be modified to run on any number of samples. The content in this script can be modified as needed, but we've provided an example on how to get started once you've downloaded all the necessary data.

## Next Steps

There is still much to be done in terms of improving the pipeline. The following are some of the next steps we hope to take in the future:

- Deploy a workflow manager to automate the pipeline
- Integrate the Cancer Predisposition Sequencing Reporter (CPSR), which annotates for germline variants
  - It is compatible with PCGR, as its output can be used to inform PCGR of the germline variants, providing another layer of assurance that somatic variants are being called
- Bioinformatics can only go so far without experimental results. As we continue to develop <i>PROTEUS</i> in the lab, we can better finetune the pipeline to suit our needs