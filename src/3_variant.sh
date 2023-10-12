#!/bin/bash
#SBATCH --mem-per-cpu=32G
#SBATCH --time=80:0:0

module load gatk/4.2.5.0
module load tabix/0.2.6

echo "Starting run at: `date`"

cd ../scripts

# Static
REF=#path_to_ref
GERM=#path_to_germline_resource
PON=#path_to_PoN

# Reading from command line
tumor=$1	# path_to_tumor
normal=$2	# path_to_normal
sample=$3 	# sample_name
out=$4		# path_to_output

./variant.sh \
	-r $REF \
	-t $tumor \
	-n $normal \
	-s $sample \
	-o $out \
	-g $GERM \
	-p $PON

echo "Job finished with exit code $? at: `date`"
