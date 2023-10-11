#!/bin/bash

folder=$1
name=$2

pcgr \
    --input_vcf ../data/$folder/variant/bqsr/variant_call.filtered.vcf.gz \
    --pcgr_dir ../pcgr \
    --output_dir ../data/$folder/variant/pcgr \
    --genome_assembly grch38 \
    --sample_id $name \
    --tumor_site 19 \
    --assay WES