#!/bin/bash
#SBATCH --mem-per-cpu=12G
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:0:0
#SBATCH --mail-user=jennifer.tramsu@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=clean.out

module purge
module load gatk/4.2.5.0
module load samtools/1.17
module load r/4.3.1

echo "Starting run at: `date`"

cd ../scripts

n_path=../data/normal/
t_path=../data/tumour/
ref=../data/reference/Homo_sapiens_assembly38.fasta
tmp=~/scratch

echo -e "\nNo bqsr...\n"

dir=raw

./clean.sh \
    -r $ref \
    -d $tmp \
    -f $dir \
    -n $n_path/n_aligned.sam \
    -t $t_path/t_aligned.sam

echo -e "\nDone!\n"

echo -e "\nBqsr...\n"

dir=bqsr

./clean.sh \
    -r $ref \
    -d $tmp \
    -f $dir \
    -c \
    -t $t_path/t_aligned.sam

echo -e "\nDone!\n"

echo "Job finished with exit code $? at: `date`"
