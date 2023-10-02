#!/bin/bash

function getopts-extra () {
    declare i=1
    # if the next argument is not an option, then append it to array OPTARG
    while [[ ${OPTIND} -le $# && ${!OPTIND:0:1} != '-' ]]; do
        OPTARG[i]=${!OPTIND}
        let i++ OPTIND++
    done
}

# excluding based on LQ scores - can also do after alignment
# q scores for alignment as well - filter out LQ alignments, FP later on

PAIRED=false

while getopts "r:n:t:p" option; do
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
        n)  getopts-extra "$@"
            args=( "${OPTARG[@]}" )
            n1=${args[0]}
            n2=${args[1]}

            n_file_path=$(dirname $n1)
            echo -e "\nAligning normal sample to reference genome...\n"

            if [ -z $n2 ] # only one file, paired ends concatenated
            then
                if $PAIRED
                then
                    bwa mem -p -M $ref $n1 > $n_file_path/n_aligned.sam
                else
                    bwa mem $ref $n1 > $n_file_path/n_aligned.sam
                fi
            else # two files
                bwa mem $ref $n1 $n2 > $n_file_path/n_aligned.sam
            fi
            ;;
        t)  getopts-extra "$@"
            args=( "${OPTARG[@]}" )
            t1=${args[0]}
            t2=${args[1]}

            t_file_path=$(dirname $t1)
            echo -e "\nAligning tumour sample to reference genome...\n"

            if [ -z $t2 ]; then
                if $PAIRED; 
                then
                    bwa mem -p -M $ref $t1 > $t_file_path/t_aligned.sam
                else
                    bwa mem $ref $t1 > $t_file_path/t_aligned.sam
                fi
            else
                bwa mem $ref $t1 $t2 > $t_file_path/t_aligned.sam
            fi
            ;;
    esac
done

echo "Done!"