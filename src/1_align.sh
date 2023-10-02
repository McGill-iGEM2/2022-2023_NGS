#!/bin/bash
#SBATCH --mem-per-cpu=12G
#SBATCH --ntasks-per-node=1
#SBATCH --time=60:0:0
#SBATCH --mail-user=jennifer.tramsu@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=align.out

module purge
module load bwa/0.7.17
module load gatk/4.2.5.0
module load samtools/1.17

echo "Starting run at: `date`"

cd ../scripts

./alignment.sh \
    -r ../data/reference/Homo_sapiens_assembly38.fasta \
    -n ../data/normal/SRR12331371_1.fastq.gz ../data/normal/SRR12331371_2.fastq.gz \
    -t ../data/tumour/SRR12331408_1.fastq.gz  ../data/tumour/SRR12331408_2.fastq.gz

echo "Job finished with exit code $? at: `date`"
