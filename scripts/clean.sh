#!/bin/bash

clean=false
# gatk - filtering LQ reads or LQ alignments, or just in CL of BWA for this filtering step
# typically filtering done more upstream here?
# maybe hidden here?
while getopts "n:t:d:r:cf:" option; do
    case $option in
    	r) ref=$OPTARG;;
        d) tmp=$OPTARG;;
        c) clean=true;;
        f) dir=$OPTARG;;
        n) normal=$OPTARG

        n_path=$(dirname $normal)/$dir
        echo "Cleaning normal sample..."
        echo "Converting from SAM to BAM..."

        samtools view -bS $normal > $n_path/n_aligned.bam

        echo -e "Query-sorting BAM file...\n"

        gatk SortSam \
            INPUT=$n_path/n_aligned.bam \
            OUTPUT=$n_path/n_qsorted_aligned.bam \
            SORT_ORDER=queryname \
            TMP_DIR=$tmp

        gatk AddOrReplaceReadGroups \
            I=$n_path/n_qsorted_aligned.bam \
            O=$n_path/n_RR_qsorted_aligned.bam \
            RGID=normal \
            RGLB=lib1 \
            RGPL=illumina \
            RGPU=unit1 \
            RGSM=normal \
            TMP_DIR=$tmp

        echo -e "\nMarking and removing duplicates...\n"

        gatk MarkDuplicates \
            INPUT=$n_path/n_RR_qsorted_aligned.bam \
            OUTPUT=$n_path/n_dup.bam \
            METRICS_FILE=$n_path/n_dup_metrics.txt \
            REMOVE_DUPLICATES=true \
            ASSUME_SORT_ORDER=queryname

        echo -e "\nCoordinate-sorting BAM file...\n"

        gatk SortSam \
            INPUT=$n_path/n_dup.bam \
            OUTPUT=$n_path/n_csorted.bam \
            SORT_ORDER=coordinate \
            TMP_DIR=$tmp

        if $clean
        then
            echo -e "\nRecalibrating base scores...\n"
            echo -e "\nFirst pass...\n"

            gatk BaseRecalibrator \
                -I $n_path/n_csorted.bam \
                -R $ref \
                -O $n_path/n_bqsr.table \
                --known-sites ../../data/Homo_sapiens_assembly38.dbsnp138.vcf \
                --known-sites ../../data/Homo_sapiens_assembly38.known_indels.vcf.gz \
                --known-sites ../../data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

            echo -e "\nSecond pass...\n"

            gatk ApplyBQSR \
                -I $n_path/n_csorted.bam \
                -R $ref \
                --bqsr-recal-file $n_path/n_bqsr.table \
                -O $n_path/n_bqsr.bam

            gatk BaseRecalibrator \
                -I $n_path/n_bqsr.bam \
                -R $ref \
                -O $n_path/n_post_bqsr.table \
                --known-sites ../../data/Homo_sapiens_assembly38.dbsnp138.vcf \
                --known-sites ../../data/Homo_sapiens_assembly38.known_indels.vcf.gz \
                --known-sites ../../data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

            echo -e "\nGenerating before/after plots...\n"

            gatk AnalyzeCovariates \
                        -before $n_path/n_bqsr.table \
                        -after $n_path/n_post_bqsr.table \
                        -plots $n_path/recalibration_plots.pdf

            echo -e "\nIndexing BAM file...\n"

            samtools index $n_path/n_bqsr.bam
        else
            echo -e "\nIndexing BAM file...\n"

            samtools index $n_path/n_csorted.bam
        fi
        ;;

        t) tumour=$OPTARG
        
        t_path=$(dirname $tumour)/$dir
        echo "Cleaning tumour sample..."
        echo "Converting from SAM to BAM..."

        samtools view -bS $tumour > $t_path/t_aligned.bam

        echo -e "\nQuery-sorting BAM file...\n"

        gatk SortSam \
            INPUT=$t_path/t_aligned.bam \
            OUTPUT=$t_path/t_qsorted_aligned.bam \
            SORT_ORDER=queryname \
            TMP_DIR=$tmp

        gatk AddOrReplaceReadGroups \
            I=$t_path/t_qsorted_aligned.bam \
            O=$t_path/t_RR_qsorted_aligned.bam \
            RGID=tumour \
            RGLB=lib1 \
            RGPL=illumina \
            RGPU=unit1 \
            RGSM=tumour \
            TMP_DIR=$tmp

        echo -e "\nMarking and removing duplicates...\n"

        gatk MarkDuplicates \
            INPUT=$t_path/t_RR_qsorted_aligned.bam \
            OUTPUT=$t_path/t_dup.bam \
            METRICS_FILE=$t_path/t_dup_metrics.txt \
            REMOVE_DUPLICATES=true \
            ASSUME_SORT_ORDER=queryname

        echo -e "\nCoordinate-sorting BAM file...\n"

        gatk SortSam \
            INPUT=$t_path/t_dup.bam \
            OUTPUT=$t_path/t_csorted.bam \
            SORT_ORDER=coordinate \
            TMP_DIR=$tmp

        if $clean
        then
            echo -e "\nRecalibrating base scores...\n"
            echo -e "\nFirst pass...\n"

            gatk BaseRecalibrator \
                -I $t_path/t_csorted.bam \
                -R $ref \
                -O $t_path/t_bqsr.table \
                --known-sites ../../data/Homo_sapiens_assembly38.dbsnp138.vcf \
                --known-sites ../../data/Homo_sapiens_assembly38.known_indels.vcf.gz \
                --known-sites ../../data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

            echo -e "\nSecond pass...\n"

            gatk ApplyBQSR \
                -I $t_path/t_csorted.bam \
                -R $ref \
                --bqsr-recal-file $t_path/t_bqsr.table \
                -O $t_path/t_bqsr.bam

            gatk BaseRecalibrator \
                -I $t_path/t_bqsr.bam \
                -R $ref \
                -O $t_path/t_post_bqsr.table \
                --known-sites ../../data/Homo_sapiens_assembly38.dbsnp138.vcf \
                --known-sites ../../data/Homo_sapiens_assembly38.known_indels.vcf.gz \
                --known-sites ../../data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

            echo -e "\nGenerating before/after plots...\n"

            gatk AnalyzeCovariates \
                        -before $t_path/t_bqsr.table \
                        -after $t_path/t_post_bqsr.table \
                        -plots $t_path/recalibration_plots.pdf

            echo -e "\nIndexing BAM file...\n"

            samtools index $t_path/t_bqsr.bam
        else
            echo -e "\nIndexing BAM file...\n"

            samtools index $t_path/t_csorted.bam
        fi
        ;;
    esac
done

echo "Done!"
