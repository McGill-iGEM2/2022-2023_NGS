#!/bin/bash

# Function for usage
usage() {
    echo -e "Usage: variant.sh <-r path_to_ref> <-t path_to_tumor> [-n path_to_normal] [-s sample] [-o out_dir] [-g path_to_germline] [-p path_to_panel] [-d tmp_dir] [-h]\n"
            echo -e "Calls variants on tumor samples using matched normal samples.\n"
            echo -e "Options:"    
            echo -e "\t-r path_to_ref\t\tReference genome"
            echo -e "\t-t path_to_tumor\tTumor sample"
            echo -e "\t-n path_to_normal\tNormal sample. If not specified, tumor-only mode is used. (Optional)"
            echo -e "\t-s sample\t\tSample name. If in matched tumor-normal mode and not specified, default is name of the normal sample. (Optional)"
            echo -e "\t-o out_dir\t\tOutput directory. Default: same as input directory. (Optional)"
            echo -e "\t-g path_to_germline\tPath to germline resource. Default: None. (Optional)"
            echo -e "\t-p path_to_panel\tPath to panel of normals. Default: None. (Optional)"
            echo -e "\t-d tmp_dir\t\tTemporary directory. Default: /tmp. (Optional)"
            echo -e "\t-h\t\t\tHelp"
}

# Default parameters
TMP=/tmp

while getopts "r:t:n:s:o:g:p:d:h" option; do
    case $option in
        r)  ref=$OPTARG
            ;;
        t)  tumor=$OPTARG
            ;;
        n)  normal=$OPTARG
            ;;
        s)  sample=$OPTARG
            ;;
        o)  out=$OPTARG
            ;;
        g)  germline=$OPTARG
            # Index germline if it doesn't exist
            if [ ! -f $germline.tbi ]; then
                echo -e "\nIndexing feature files...\n"
                gatk IndexFeatureFile -I $germline
            fi
            ;;
        p)  panel=$OPTARG
            # Index panel of normals if it doesn't exist
            if [ ! -f $panel.tbi ]; then
                echo -e "\nIndexing feature files...\n"
                gatk IndexFeatureFile -I $panel
            fi
            ;;
        d)  TMP=$OPTARG
            ;;
        h)
            usage
            exit 1
            ;;
    esac
done

# Extracting file path
path=$(dirname $tumor)

# If a normal is present and no sample is given, extract the sample name from the normal
if [ ! -z $normal ] && [ -z $sample ]; then
    sample=$(basename $normal)
fi

# Make output directory if it doesn't exist
if [ ! -d $out ]; then
    mkdir $out
fi

# Update path if output directory is specified
if [ ! -z $out ]; then
    path=$out
fi

# Error message if no arguments are given
if [ $# -eq 0 ]; then
    echo -e "Error: No arguments given\n"
    usage
    exit 1
elif [ -z $ref ]; then
    echo -e "Error: Reference genome not specified\n"
    usage
    exit 1
elif [ -z $tumor ]; then
    echo -e "Error: Path to tumor sample not specified\n"
    usage
    exit 1
fi

# Start variant calling pipeline
echo -e "\nCalling variants on $sample...\n"

if [ ! -z $normal ]; then # Matched tumor-normal mode
    echo -e "\nMatched tumor-normal mode\n"
    
    # Germline and panel
    if [ ! -z $germline ] && [ ! -z $panel ]; then
        gatk Mutect2 \
            --java-options "-Xmx32G" \
            -R $ref \
            -I $normal \
            -I $tumor \
            -normal $sample \
            --germline-resource $germline \
            --panel-of-normals $panel \
            -O $path/${sample}_variant_call.vcf \
            --tmp-dir $TMP
    elif [ ! -z $germline ] && [ -z $panel ]; then # Only germline
        gatk Mutect2 \
            --java-options "-Xmx32G" \
            -R $ref \
            -I $normal \
            -I $tumor \
            -normal $sample \
            --germline-resource $germline \
            -O $path/${sample}_variant_call.vcf \
            --tmp-dir $TMP
    elif [ -z $germline ] && [ ! -z $panel ]; then # Only panel
        gatk Mutect2 \
            --java-options "-Xmx32G" \
            -R $ref \
            -I $normal \
            -I $tumor \
            -normal $sample \
            --panel-of-normals $panel \
            -O $path/${sample}_variant_call.vcf \
            --tmp-dir $TMP
    elif [ -z $germline ] && [ -z $panel ]; then # Neither germline nor panel
        gatk Mutect2 \
            --java-options "-Xmx32G" \
            -R $ref \
            -I $normal \
            -I $tumor \
            -normal $sample \
            -O $path/${sample}_variant_call.vcf \
            --tmp-dir $TMP
    fi
else # Tumor-only mode
    echo -e "\nTumor-only mode\n"

    # Germline and panel
    if [ ! -z $germline ] && [ ! -z $panel ]; then
        gatk Mutect2 \
            --java-options "-Xmx32G" \
            -R $ref \
            -I $tumor \
            --germline-resource $germline \
            --panel-of-normals $panel \
            -O $path/${sample}_variant_call.vcf \
            --tmp-dir $TMP
    elif [ ! -z $germline ] && [ -z $panel ]; then # Only germline
        gatk Mutect2 \
            --java-options "-Xmx32G" \
            -R $ref \
            -I $tumor \
            --germline-resource $germline \
            -O $path/${sample}_variant_call.vcf \
            --tmp-dir $TMP
    elif [ -z $germline ] && [ ! -z $panel ]; then # Only panel
        gatk Mutect2 \
            --java-options "-Xmx32G" \
            -R $ref \
            -I $tumor \
            --panel-of-normals $panel \
            -O $path/${sample}_variant_call.vcf \
            --tmp-dir $TMP
    elif [ -z $germline ] && [ -z $panel ]; then # Neither germline nor panel
        gatk Mutect2 \
            --java-options "-Xmx32G" \
            -R $ref \
            -I $tumor \
            -O $path/${sample}_variant_call.vcf \
            --tmp-dir $TMP
    fi
fi

echo -e "\nFiltering vcf...\n"

gatk FilterMutectCalls \
    -R $ref \
    -V $path/${sample}_variant_call.vcf \
    -O $path/${sample}_variant_call.filtered.vcf

bgzip -c $path/${sample}_variant_call.filtered.vcf > $path/${sample}_variant_call.filtered.vcf.gz
tabix -p vcf $path/${sample}_variant_call.filtered.vcf.gz

echo -e "\nDone!\n"