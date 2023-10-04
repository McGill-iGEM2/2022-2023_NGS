#!/bin/bash
#SBATCH --mem-per-cpu=12G
#SBATCH --ntasks-per-node=1
#SBATCH --time=60:0:0
#SBATCH --mail-user=jennifer.tramsu@mail.mcgill.ca
#SBATCH --mail-type=ALL

module purge
module load bwa/0.7.17
module load gatk/4.2.5.0
module load samtools/1.17

echo "Starting run at: `date`"

cd ../scripts

./alignment.sh \
    -r ../data/reference/Homo_sapiens_assembly38.fasta \
    -t $1 $2

echo "Job finished with exit code $? at: `date`"
