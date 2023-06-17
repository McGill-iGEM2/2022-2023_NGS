#!/bin/bash

while getopts "r:n:t:" options; do
    case $option in
        r) ref=$OPTARG;;
        s) normal=$OPTARG;;
        t) tumour=$OPTARG;;
    esac
done

#bwa index $ref
bwa mem $ref $normal > n_aligned.sam
bwa mem $ref $tumour > t_aligned.sam