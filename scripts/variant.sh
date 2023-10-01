#!/bin/bash

while getopts "r:n:t:o:g:p:d:" option; do
    case $option in
        r) ref=$OPTARG;;
        n) normal=$OPTARG;;
        t) tumour=$OPTARG;;
        o) out=$OPTARG;;
        g) germline=$OPTARG
            if [ ! -f $germline.tbi ]; 
            then
                echo -e "\nIndexing feature files...\n"
                gatk IndexFeatureFile -I $germline
            fi
        ;;
        p) panel=$OPTARG
            if [ ! -f $panel.tbi ];
            then
                echo -e "\nIndexing feature files...\n"
                gatk IndexFeatureFile -I $panel
            fi
        ;;
        d) tmp=$OPTARG;;
    esac
done

# job on CC - allocate more memory
# by default, java has max size to its heap - no matter how much memory on node, if it exceeds it stops
# if calling through java, option available to increase size of heap (google)
# help desk CC -> good resource
# gatk way to go
#   most commonly used

# comparing multiple variant callers
#   time consuming, gatk may have some FP
#   intersecting predictions could remove more FP but not best use of time
# FP clinical setting
#   group of researches @ mcgill -> working on this, C3G
#   can contact them
# detecting mutations?
# accurate? which are pathogenic? harder question to validate
#   variants detected by gatk come with confidence values
#   beyond that, people will visualize alignments with IGV to convince themselves that calls look reasonable (not for every variant, for those you are interested in)

# Variant calling
gatk Mutect2 \
    --java-options "-Xmx8G" \
    -R $ref \
    -I $normal \
    -I $tumour \
    -normal normal \
    -O $out.vcf \
    --tmp-dir $tmp

gatk FilterMutectCalls \
    -R $ref \
    -V $out.vcf \
    -O $out.filtered.vcf

bgzip -c $out.filtered.vcf > $out.filtered.vcf.gz
tabix -p vcf $out.filtered.vcf.gz
