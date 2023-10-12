#!/bin/bash
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0:30:0

module add mugqic/pcgr/1.2.0

echo "Starting run at: `date`"

cd ../scripts

# Static
PCGR=#path_to_pcgr

# Reading from command line
in=$1         # path_to_input
sample=$2     # sample_name
out=$3        # path_to_output

./annotate.sh \
    -i $in \
    -s $sample \
    -p $PCGR \
    -o $out

echo "Job finished with exit code $? at: `date`"