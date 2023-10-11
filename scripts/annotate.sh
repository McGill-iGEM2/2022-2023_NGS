#!/bin/bash

# Function for usage
usage() {
    echo -e "Usage: annotate.sh <-i path_to_vcf> [-s sample] [-p path_to_pcgr] [-o out_dir] [-g genome_assembly] [-t tumor_site] [-w assay] [-h]\n"
            echo -e "Annotates VCF files using PCGR given a specific tumor site.\n"
            echo -e "Options:"    
            echo -e "\t-i path_to_vcf\t\tVCF file"
            echo -e "\t-s sample\t\tSample name. Default is the name of the file containing the reads. (Optional)"
            echo -e "\t-p path_to_pcgr\t\tPath to PCGR data bundle. Default: ../pcgr. (Optional)"
            echo -e "\t-o out_dir\t\tOutput directory. Default: same as input directory. (Optional)"
            echo -e "\t-g genome_assembly\tGenome assembly, can be grch38 or grch37. Default: grch38. (Optional)"
            echo -e "\t-t tumor_site\t\tTumor site. Default: 19 (pancreatic cancer). (Optional)"
            echo -e "\t-w assay\t\tAssay. Default: WES. (Optional)"
            echo -e "\t-h\t\t\tHelp"
}

# Default parameters
ASSEMBLY=grch38
ASSAY=WES
TSITE=19 # pancreatic cancer
PCGR=../pcgr

while getopts "i:s:p:o:g:t:w:h" option; do
    case $option in
        i)  in=$OPTARG
            ;;
        s)  sample=$OPTARG
            ;;
        p)  PCGR=$OPTARG
            ;;
        o)  out=$OPTARG
            ;;
        g)  ASSEMBLY=$OPTARG
            ;;
        t)  tsite=$OPTARG
            ;;
        w)  ASSAY=$OPTARG
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
elif [ -z $in ]; then
    echo -e "Error: Path to VCF not specified\n"
    usage
    exit 1
elif [ -z $pcgr ]; then
    echo -e "Error: Path to PCGR data bundle not specified\n"
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

# Annotating VCF
echo -e "\nAnnotating $sampple...\n"

pcgr \
    --input_vcf $in \
    --pcgr_dir $pcgr \
    --output_dir $path \
    --genome_assembly $ASSEMBLY \
    --sample_id $sample \
    --tumor_site $tsite \
    --assay $ASSAY

echo -e "\nDone!\n"