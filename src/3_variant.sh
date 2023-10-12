#!/bin/bash
#SBATCH --mem-per-cpu=32G
#SBATCH --time=80:0:0
#SBATCH --mail-type=ALL

module load gatk/4.2.5.0
module load tabix/0.2.6

echo "Starting run at: `date`"

cd ../scripts

./variant.sh \
	-r # path_to_ref \
	-t # path_to_tumor \
	-n # path_to_normal \
	-s # sample_name \
	-o # path_to_out_dir \
	-g # path_to_germline_resource \
	-p # path_to_PoN

echo "Job finished with exit code $? at: `date`"
