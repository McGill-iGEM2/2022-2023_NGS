#!/bin/bash

while getopts "n:t:" option; do
    case $option in
        n) normal=$OPTARG
            n_path=$(dirname $normal)
            echo "Cleaning normal sample..."
            echo "Converting from SAM to BAM..."

            samtools view -bS $normal > $n_path/n_aligned.bam

            echo -e "Query-sorting BAM file...\n"

            picard SortSam \
                INPUT=$n_path/n_aligned.bam \
                OUTPUT=$n_path/n_qsorted_aligned.bam \
                SORT_ORDER=queryname

            echo -e "\nMarking and removing duplicates...\n"

            picard MarkDuplicates \
                INPUT=$n_path/n_qsorted_aligned.bam \
                OUTPUT=$n_path/n_dup.bam \
                METRICS_FILE=$n_path/n_dup_metrics.txt \
                REMOVE_DUPLICATES=true

            echo -e "\nCoordinate-sorting BAM file...\n"

            picard SortSam \
                INPUT=$n_path/n_dup.bam \
                OUTPUT=$n_path/n_csorted.bam \
                SORT_ORDER=coordinate

            echo -e "\nIndexing BAM file...\n"

            samtools index $n_path/n_csorted.bam
            ;;
        t) tumour=$OPTARG
            t_path=$(dirname $tumour)
            echo "Cleaning tumour sample..."
            echo "Converting from SAM to BAM..."

            samtools view -bS $tumour > $t_path/t_aligned.bam

            echo -e "\nQuery-sorting BAM file...\n"

            picard SortSam \
                INPUT=$t_path/t_aligned.bam \
                OUTPUT=$t_path/t_qsorted_aligned.bam \
                SORT_ORDER=queryname

            echo -e "\nMarking and removing duplicates...\n"

            picard MarkDuplicates \
                INPUT=$t_path/t_qsorted_aligned.bam \
                OUTPUT=$t_path/t_dup.bam \
                METRICS_FILE=$t_path/t_dup_metrics.txt \
                REMOVE_DUPLICATES=true

            echo -e "\nCoordinate-sorting BAM file...\n"

            picard SortSam \
                INPUT=$t_path/t_dup.bam \
                OUTPUT=$t_path/t_csorted.bam \
                SORT_ORDER=coordinate

            echo -e "\nIndexing BAM file...\n"

            samtools index $t_path/t_csorted.bam
            ;;
    esac
done

echo "Done!"