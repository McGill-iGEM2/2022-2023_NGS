# NGS

## Pipeline
1) Prefetching `fastq` Dataset

  - Program: SRAtoolkit
  
BAM files are converted to fastq. Datasets are downloaded through prefetch commands.

The following steps are performed on each sequence (healthy and cancer) individually.

2) Genome Alignmnet

  - Program: [BWA-MEM](https://github.com/lh3/bwa)
  - Input: fastq file
  - Output: BAM file, aligned reads
  
Read groups are aligned to the reference genome.

3) Alignment Co-Cleaning

  - Program: [GATK](https://software.broadinstitute.org/gatk/) ([IndelRealigner](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php) -> [BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator))
  - Input: BAM file, aligned reads
  - Output: BAM file, harmonized aligned reads
  
Improving alignment quality.

4) Somatic Variant Calling

  - Programs: [DeepVariant](https://github.com/google/deepvariant) or [MuTect2](https://gatk.broadinstitute.org/hc/en-us/articles/360047232772--Notebook-Intro-to-using-Mutect2-for-somatic-data) (Others: MuSE, VarScan2, Pindel)
  - Input: BAM file, aligned reads
  - Output: VCF file, raw simple somatic mutation
  
  - If using MuTect2, it takes in both the healthy and cancer sequence together to output one VCF file
  
Aligned and co-cleaned BAM files are processed through the Somatic Mutation Calling Workflow as tumor-normal pairs. Variant calling is the process by which we identify variants from sequence data.

If using MuTect2, can also annotate and filter.

5) Variant Filtering

  - Program: VarScan (filtering) or MuTect2
  - Input: VCF file, simple
  - Output: VCF file, filtered
  
Removes false positives.
  
6a) Variant Annotation

  - Program: ANNOVAR (annotation) or MuTect2
  - Input: VCF file, filtered
  - Output: VCF file, annotated somatic mutation
  
Annotates location of each mutation, its biological consequences (frameshift/silent mutation) and the affected genes.

6b) Variant Frequency

  - Program: some conversion program
  - Input: VCF file, filtered
  - Output: BAM

To determine the frequency of a variant, we can convert the filtered VCF file to BAM so that it can be processed using the GATK (toolkit).

  - Program: GATK
  - Input: BAM file
  - Output: Variant frequency (hopefully)

Find # of reads.
