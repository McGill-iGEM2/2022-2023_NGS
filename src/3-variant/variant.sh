#!/bin/bash

while getopts "r:n:t:" options; do
    case $option in
        r) ref=$OPTARG;;
        n) normal=$OPTARG;;
        t) tumour=$OPTARG;;
    esac
done

gatk Mutect2 \
    -R $ref \
    -I $normal \
    -I $tumour \
    -O variant_call.vcf.gz