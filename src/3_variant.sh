#!/bin/bash
#SBATCH --mem=8G
#SBATCH --time=80:0:0
#SBATCH --mail-user=jennifer.tramsu@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=variant.out

module purge
module load gatk/4.2.5.0
module load tabix/0.2.6

echo "Starting run at: `date`"

cd ../scripts

echo -e "\nNo bqsr...\n"

./variant.sh \
	-r ../data/reference/Homo_sapiens_assembly38.fasta \
	-n ../data/normal/raw/n_csorted.bam \
	-t ../data/tumour/raw/t_csorted.bam \
	-g ../data/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
	-p ../data/somatic-hg38_1000g_pon.hg38.vcf.gz \
	-o ../data/variant/test/raw/variant_call \
	-d ~/scratch

echo -e "\nBqsr...\n"

./variant.sh \
	-r ../data/reference/Homo_sapiens_assembly38.fasta \
	-n ../data/normal/bqsr/n_bqsr.bam \
	-t ../data/tumour/bqsr/t_bqsr.bam \
	-g ../data/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
	-p ../data/somatic-hg38_1000g_pon.hg38.vcf.gz \
	-o ../data/variant/test/bqsr/variant_call \
	-d ~/scratch

echo "Job finished with exit code $? at: `date`"
