# McGill iGEM NGS Bioinformatics Pipeline
A next-generation sequencing bioinformatics pipeline developed for the discovery and detection of somatic SNPs and indels in cancer patients. This repository consists of individual `bash` files that can be run locally as well as a series of `SLURM` scripts for execution on the cloud.

## Repository Structure
```
.
├── README.md
├── data
│   ├── normal
│   ├── reference
│   ├── tumour
│   └── variant
├── scripts
│   ├── alignment.sh
│   ├── clean.sh
│   ├── variant.sh
└── src
    ├── 1_align.sh
    ├── 2_clean.sh
    ├── 3_variant.sh
    └── 4_annotate.sh
  ```

## Setup
### Local Installation
Many of the software used are available through `bioconda`, a distribution of bioinformatics software as part of the Conda package manager.

- miniconda https://docs.conda.io/projects/miniconda/en/latest/
- bwa https://anaconda.org/bioconda/bwa
- samtools https://anaconda.org/bioconda/samtools
- gatk https://anaconda.org/bioconda/gatk4
