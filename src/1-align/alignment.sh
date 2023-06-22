#!/bin/bash

while getopts "r:n:t:" option; do
    case $option in
        r) ref=$OPTARG
            if [ ! -f $ref.fai ]; then
                echo -e "Indexing reference genome...\n"
                bwa index $ref
            fi
            ;;
        n) normal=$OPTARG
            n_file_path=$(dirname $normal)
            echo -e "\nAligning normal sample to reference genome...\n"
            bwa mem $ref $normal > $n_file_path/n_aligned.sam
            ;;
        t) tumour=$OPTARG
            t_file_path=$(dirname $tumour)
            echo -e "\nAligning tumour sample to reference genome...\n"
            bwa mem $ref $tumour > $t_file_path/t_aligned.sam
            ;;
    esac
done

echo "Done!"