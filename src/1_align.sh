#!/bin/bash
#SBATCH --mem-per-cpu=12G
#SBATCH --time=60:0:0
#SBATCH --mail-type=ALL

module load bwa/0.7.17
module load gatk/4.2.5.0
module load samtools/1.17

echo "Starting run at: `date`"

cd ../scripts

./alignment.sh \
    -r # path_to_ref \
    -n # path_1 path_2 \
    -s # sample_name

echo "Job finished with exit code $? at: `date`"
