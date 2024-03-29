#!/bin/bash

# Function for usage
usage() {
    echo -e "Usage: clean.sh <-r path_to_ref> <-i path_to_aligned_reads> [-s sample] [-d tmp_dir] [-o out_dir] [-b] [-h]\n"
            echo -e "Converts SAM files to BAM and cleans them by removing duplicates and recalibrating base scores.\n"
            echo -e "Options:"    
            echo -e "\t-r path_to_ref\t\t\tReference genome"
            echo -e "\t-i path_to_aligned_reads\tAligned reads"
            echo -e "\t-s sample\t\t\tSample name. Default is the name of the file containing the reads. (Optional)"
            echo -e "\t-d tmp_dir\t\t\tTemporary directory. Default: /tmp. (Optional)"
            echo -e "\t-o out_dir\t\t\tOutput directory. Default: same as input directory. (Optional)"
            echo -e "\t-b\t\t\t\tUse this flag to indicate that base scores should be recalibrated. Default is no recalibration. (Optional)"
            echo -e "\t-h\t\t\t\tHelp"
}

# Default parameters
BQSR=false
TMP=/tmp

while getopts "r:i:s:o:bd:h" option; do
    case $option in
        r)  ref=$OPTARG
            ;;
        i)  in=$OPTARG
            ;;
        s)  sample=$OPTARG
            ;;
        o)  out=$OPTARG
            ;;       
        b)  BQSR=true
            ;;
        d)  TMP=$OPTARG
            ;;
        h) 
            usage
            exit 1
            ;;
    esac
done

# Error message if no arguments are given
if [ $# -eq 0 ]; then
    echo -e "Error: No arguments given\n"
    usage
    exit 1
elif [ -z $ref ]; then
    echo -e "Error: Path to reference genome not specified\n"
    usage
    exit 1
elif [ -z $in ]; then
    echo -e "Error: Path to aligned reads not specified\n"
    usage
    exit 1
fi

# Extracting file path
path=$(dirname $in)

# Extracting sample name if not specified
if [ -z $sample ]; then
    sample=$(basename $in)
fi

# Make output directory if it doesn't exist
if [ ! -d $out ]; then
    mkdir $out
fi

# Update path if output directory is specified
if [ ! -z $out ]; then
    path=$out
fi

# Start cleaning pipeline
echo "Cleaning $sample..."
echo "Converting from SAM to BAM..."

samtools view -bS $in > $path/${sample}_aligned.bam

echo -e "Query-sorting BAM file...\n"

gatk SortSam \
    INPUT=$path/${sample}_aligned.bam \
    OUTPUT=$path/${sample}_qsorted_aligned.bam \
    SORT_ORDER=queryname \
    TMP_DIR=$TMP

gatk AddOrReplaceReadGroups \
    I=$path/${sample}_sorted_aligned.bam \
    O=$path/${sample}_RR_qsorted_aligned.bam \
    RGID=$sample \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=$sample \
    TMP_DIR=$TMP

echo -e "\nMarking and removing duplicates...\n"

gatk MarkDuplicates \
    INPUT=$path/${sample}_RR_qsorted_aligned.bam \
    OUTPUT=$path/${sample}_dup.bam \
    METRICS_FILE=$path/${sample}_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    ASSUME_SORT_ORDER=queryname

echo -e "\nCoordinate-sorting BAM file...\n"

gatk SortSam \
    INPUT=$path/${sample}_dup.bam \
    OUTPUT=$path/${sample}_csorted.bam \
    SORT_ORDER=coordinate \
    TMP_DIR=$TMP

if $BQSR
then
    echo -e "\nRecalibrating base scores...\n"
    echo -e "\nFirst pass...\n"

    gatk BaseRecalibrator \
        -I $path/${sample}_csorted.bam \
        -R $ref \
        -O $path/${sample}_bqsr.table \
        --known-sites ../data/Homo_sapiens_assembly38.dbsnp138.vcf \
        --known-sites ../data/Homo_sapiens_assembly38.known_indels.vcf.gz \
        --known-sites ../data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

    echo -e "\nSecond pass...\n"

    gatk ApplyBQSR \
        -I $path/${sample}_csorted.bam \
        -R $ref \
        --bqsr-recal-file $path/${sample}_bqsr.table \
        -O $path/${sample}_bqsr.bam

    gatk BaseRecalibrator \
        -I $path/${sample}_bqsr.bam \
        -R $ref \
        -O $path/n_post_bqsr.table \
        --known-sites ../data/Homo_sapiens_assembly38.dbsnp138.vcf \
        --known-sites ../data/Homo_sapiens_assembly38.known_indels.vcf.gz \
        --known-sites ../data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

    echo -e "\nGenerating before/after plots...\n"

    gatk AnalyzeCovariates \
                -before $path/${sample}_bqsr.table \
                -after $path/${sample}_post_bqsr.table \
                -plots $path/recalibration_plots.pdf

    echo -e "\nIndexing BAM file...\n"

    samtools index $path/${sample}_bqsr.bam
else
    echo -e "\nIndexing BAM file...\n"

    samtools index $path/${sample}_csorted.bam
fi
;;

echo -e "\nDone!\n"