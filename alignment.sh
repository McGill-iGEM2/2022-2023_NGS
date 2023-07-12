bwa mem \
# invokes the BWA-MEM algorithm
-t 8 \
# Specifies the number of threads (parallel processes) to be used by BWA-MEM during the alignment
-T 0 \
# Sets the score threshold for a primary alignment
-R <read_group> \
/Users/adminelya/Documents/GitHub/NGS/data/reference\ genome/GRCh38ref.gz.sa \
/Users/adminelya/Documents/GitHub/NGS/data/fastq/SRR12331326_GL000219.1.fastq |
samtools view \
-Shb
-o /Users/adminelya/Documents/GitHub/NGS/data/aligned.sam
# transforms to bam

