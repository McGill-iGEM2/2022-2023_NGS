#!/bin/bash
#SBATCH --mem-per-cpu=12G
#SBATCH --time=120:0:0

module load gatk/4.2.5.0
module load samtools/1.17
module load r/4.3.1

echo "Starting run at: `date`"

cd ../scripts

# Static
ref=#path_to_ref

# Reading from command line
in=$1         # path_to_input
sample=$2     # sample_name

./clean.sh \
    -r $ref \
    -i $in \
    -s $sample \
    -b 

echo -e "\nDone!\n"

echo "Job finished with exit code $? at: `date`"
