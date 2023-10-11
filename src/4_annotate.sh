#!/bin/bash
#SBATCH --mem-per-cpu=12G
#SBATCH --ntasks-per-node=1
#SBATCH --time=0:30:0
#SBATCH --mail-user=jennifer.tramsu@mail.mcgill.ca
#SBATCH --mail-type=ALL

module add mugqic/pcgr/1.2.0

echo "Starting run at: `date`"

cd ../scripts
folder=$1
name=$2
./annotate.sh $folder $name

echo "Job finished with exit code $? at: `date`"