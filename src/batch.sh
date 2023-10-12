#!/bin/bash

# To use this script, we've provided an example based on our own workflow
# Ideally, we would like to implement a workflow manager that runs each for loop ONLY when all previous for loop have finished
# But at the time of writing, we have been manually commenting the loops not in use :p

dir=(PAM43_MET_AOR PAM43_PRIM PAM43_MET_PAN)
accession=(SRR12331374 SRR12331371 SRR12331376)
name=(AOR PRIM PAN)

for (( i=0; i<${#dir[@]}; i++ ));
do
    folder=${dir[i]}
    sample=${name[i]}
    id=${accession[i]}

    sbatch --mail-type=ALL --mail-user=(your email here) --job-name "align $sample" ./1_align.sh ../data/$folder/tumour/${id}_1.fastq ../data/$folder/tumour/${id}_2.fastq $sample
done

for (( i=0; i<${#dir[@]}; i++ ));
do
    folder=${dir[i]}
    sample=${name[i]}

    sbatch --mail-type=ALL --mail-user=(your email here) --job-name "clean $sample" ./2_clean.sh ../data/$folder/tumour/${sample}_aligned.sam $sample
done

for (( i=0; i<${#dir[@]}; i++ ));
do
    folder=${dir[i]}
    sample=${name[i]}

    sbatch --mail-type=ALL --mail-user=(your email here) --job-name "variant $sample" ./3_variant.sh ../data/$folder/tumour/${sample}_clean.bam ../data/normal/${sample}_normal.bam $sample ../data/$folder/variant
done

for (( i=0; i<${#dir[@]}; i++ ));
do
    folder=${dir[i]}
    sample=${name[i]}

    sbatch --mail-type=ALL --mail-user=(your email here) --job-name "annotate $sample" ./4_annotate.sh ../data/$folder/variant/${sample}_variant_call.filtered.vcf.gz $sample ../data/$folder/variant/pcgr
done