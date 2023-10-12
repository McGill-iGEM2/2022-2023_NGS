#!/bin/bash
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0:30:0
#SBATCH --mail-type=ALL

module add mugqic/pcgr/1.2.0

echo "Starting run at: `date`"

cd ../scripts

./annotate.sh \
    -i # path_to_vcf \
    -s # sample_name \
    -p # path_to_pcgr \
    -o # path_to_out_dir

echo "Job finished with exit code $? at: `date`"