#!/bin/bash

# Function for usage
usage() {
    echo -e "Usage: alignment.sh <-r path_to_ref> <-n path_to_reads [path_to_paired_reads]> [-p] [-h]\n"
            echo -e "Aligns sequencing data to the reference genome using BWA-MEM. The reference is indexed if not already.\n"
            echo -e "Options:"    
            echo -e "\t-r path_to_ref\t\t\t\tReference genome"
            echo -e "\t-n path_to_reads [path_to_paired_reads]\tSequencing data. If paired-end reads contained in two files, indicate the pair file with RIGHT"
            echo -e "\t-s sample\t\t\t\tSample name. Default is the name of the file containing the reads. (Optional)"
            echo -e "\t-p\t\t\t\t\tUse this flag to indicate that READS contains paired-end reads. Default is single-end reads. (Optional)"
            echo -e "\t-h\t\t\t\t\tHelp"
}

# Function for getting extra arguments
function getopts-extra () {
    declare i=1
    # if the next argument is not an option, then append it to array OPTARG
    while [[ ${OPTIND} -le $# && ${!OPTIND:0:1} != '-' ]]; do
        OPTARG[i]=${!OPTIND}
        let i++ OPTIND++
    done
}

# Default parameter
PAIRED=false

while getopts "r:n:s:ph" option; do
    case $option in
        r)  ref=$OPTARG
            ;;
        n)  getopts-extra "$@"
            args=( "${OPTARG[@]}" )
            left=${args[0]}
            right=${args[1]}
            ;;
        s)  sample=$OPTARG
            ;;
        p)  PAIRED=true
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
elif [ -z $ref ] || [ -z $left ]; then
    echo -e "Error: Reference genome or sample not specified\n"
    usage
    exit 1
elif [ -z $sample ]; then
    sample=$(basename $left)
fi

# Extracting file path
path=$(dirname $left)

# Indexing reference genome if not already
if [ ! -f $ref.bwt ]; then
    echo -e "Indexing reference genome...\n"
    bwa index $ref
    samtools faidx $ref
fi

# Aligning sample to reference genome
echo -e "\nAligning $sample to reference genome...\n"

if [ -z $right ] # only one file, paired ends concatenated
    then
        if $PAIRED
        then
            bwa mem -p -M $ref $left > $path/${sample}_aligned.sam
        else # single end reads
            bwa mem $ref $left > $path/${sample}_aligned.sam
        fi
else # two files
    bwa mem $ref $left $right > $path/${sample}_aligned.sam
fi

echo "Done!"