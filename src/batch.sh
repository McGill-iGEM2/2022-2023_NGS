#!/bin/bash

dir=(PAM43_PRIM PAM43_MET_AOR PAM43_MET_LIV PAM43_MET_LUNG PAM43_MET_MES PAM43_MET_PAN)
name=(PRIM AOR LIV LUNG MES PAN)
accession=(SRR12331374 SRR12331375 SRR12331372 SRR12331373 SRR12331376)

# for (( i=0; i<${#dir[@]}; i++ ));
# do
#     folder=${dir[i]}
#     id=${accession[i]}

#     sbatch --job-name "align $(basename $folder)" ./1_align.sh ../data/$folder/tumour/${id}_1.fastq ../data/$folder/tumour/${id}_2.fastq
# done

for (( i=0; i<${#dir[@]}; i++ ));
do
    folder=${dir[i]}
    sample=${name[i]}
    sbatch --job-name "clean $(basename $sample)" ./2_clean.sh $folder
done

# for (( i=0; i<${#dir[@]}; i++ ));
# do
#     folder=${dir[i]}
#     sbatch --job-name "variant $(basename $folder)" ./3_variant.sh $folder
# done