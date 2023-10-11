#!/bin/bash

dir=(PAM43_MET_AOR PAM43_PRIM PAM43_MET_PAN)
name=(AOR PRIM PAN)

# for (( i=0; i<${#dir[@]}; i++ ));
# do
#     folder=${dir[i]}
#     sample=${name[i]}
#     id=${accession[i]}

#     sbatch --job-name "align $(basename $sample)" ./1_align.sh ../data/$folder/tumour/${id}_1.fastq ../data/$folder/tumour/${id}_2.fastq
# done

# for (( i=0; i<${#dir[@]}; i++ ));
# do
#     folder=${dir[i]}
#     sample=${name[i]}
#     sbatch --job-name "clean $(basename $sample)" ./2_clean.sh $folder
# done

# for (( i=0; i<${#dir[@]}; i++ ));
# do
#     folder=${dir[i]}
#     sample=${name[i]}
#     sbatch --job-name "variant $(basename $sample)" ./3_variant.sh $folder
# done

for (( i=0; i<${#dir[@]}; i++ ));
do
    folder=${dir[i]}
    sample=${name[i]}
    sbatch --job-name "annotate $(basename $sample)" ./4_annotate.sh $folder $sample
done