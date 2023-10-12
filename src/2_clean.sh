#!/bin/bash
#SBATCH --mem-per-cpu=12G
#SBATCH --time=120:0:0
#SBATCH --mail-type=ALL

module load gatk/4.2.5.0
module load samtools/1.17
module load r/4.3.1

echo "Starting run at: `date`"

cd ../scripts

./clean.sh \
    -r # path_to_ref \
    -i # path_to_sam \
    -s # sample_name \
    -o # path_to_out_dir \
    -b 

echo -e "\nDone!\n"

echo "Job finished with exit code $? at: `date`"
