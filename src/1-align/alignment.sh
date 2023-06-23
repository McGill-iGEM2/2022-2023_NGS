#!/bin/bash

while getopts "r:pn:t:" option; do
    case $option in
        r) ref=$OPTARG
            if [ ! -f $ref.bwt ]; then
                echo -e "Indexing reference genome...\n"
                bwa index $ref
                samtools faidx $ref
            fi
            ;;
        p)  PAIRED=true
            ;;
        n)  set -f
            IFS=' '
            shift $(($OPTIND-2))

            n_file_path=$(dirname $1)
            echo -e "\nAligning normal sample to reference genome...\n"

            if [ -z $2 ] # only one file, paired ends concatenated
            then
                if $PAIRED ; then
                    bwa mem -p $ref $1 > $n_file_path/n_aligned.sam
                else
                    bwa mem $ref $1 > $n_file_path/n_aligned.sam
                fi
            else # two files
                bwa mem $ref $1 $2 > $n_file_path/n_aligned.sam
            fi
            ;;
        t)  set -f
            IFS=' '
            shift $(($OPTIND-2))

            t_file_path=$(dirname $1)
            echo -e "\nAligning tumour sample to reference genome...\n"

            if [ -z $2 ]; then
                if [ $PAIRED ]; then
                    bwa mem -p $ref $1 > $t_file_path/t_aligned.sam
                else
                    bwa mem $ref $1 > $t_file_path/t_aligned.sam
                fi
            else
                bwa mem $ref $1 $2 > $t_file_path/t_aligned.sam
            fi
            ;;
    esac
done

echo "Done!"

# TODO
# picard addorremovereadgroups
# testing mutect2