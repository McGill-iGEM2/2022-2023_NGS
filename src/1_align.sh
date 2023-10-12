#!/bin/bash
#SBATCH --mem-per-cpu=12G
#SBATCH --time=60:0:0

module load bwa/0.7.17
module load gatk/4.2.5.0
module load samtools/1.17

echo "Starting run at: `date`"

cd ../scripts

# Static
ref=#path_to_ref

# Reading from command line
left=$1
right=$2        # path_to_input
sample=$3       # sample_name

./alignment.sh \
    -r $ref \
    -n $left $right \
    -s $sample

echo "Job finished with exit code $? at: `date`"
