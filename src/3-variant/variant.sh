#!/bin/bash

while getopts "r:n:t:o:" option; do
    case $option in
        r) ref=$OPTARG;;
        n) normal=$OPTARG;;
        t) tumour=$OPTARG;;
        o) out=$OPTARG;;
    esac
done

gatk Mutect2 \
    -R $ref \
    -I $normal \
    -I $tumour \
    -normal normal \
    -O $out