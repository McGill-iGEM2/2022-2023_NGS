#!/bin/bash
#SBATCH --mem-per-cpu=12G
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:0:0
#SBATCH --mail-user=jennifer.tramsu@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=annotate.out

echo "Starting run at: `date`"

cd ../data/variant/bqsr
ls -l

singularity exec ~/vep.sif vep \
	--dir $HOME/vep_data \
	--dir_cache $HOME/vep_data \
	-i /home/jts/jts/NGS/data/variant/raw/variant_call.vcf \
	--cache \
	--offline \
	--force_overwrite \
	--everything \
	--tab \
	-o $HOME/jts/NGS/data/variant/raw/variant_output.txt
