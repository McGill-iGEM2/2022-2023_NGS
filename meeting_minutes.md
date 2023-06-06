## Monday, June 5, 2023
Updates:
- Pilar found that marked dupliciates are useful but depends on the data. More useful in case there are issues with the data itself. 
    - A couple ways to do duplicates Mark Duplicate Spark: Reqires a lot of disk operations but optimized to run locally as best as possible (at least 16GB RAM recommended)
        - Compute Canada??!? for those of us who have weak computers
- BQSR (Base recalibration, finetuning quality scores) not 100% necessary but could give more detailed metadata.
- [Useful review paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8506910/) for an evaluation of the whole pipeline - Thanks Albert!
- TCGA has atlas for pancreatic cell lines and expression levels
- SHAD Python test run happening tomorrow

Problems/Tasks:
- Read the above paper depth to see if we're missing any steps or could be incorporating anything to improve results (@pilar, but everybody should read it)
- Try other parts of GATK to see if it's working/we can use it for other things. Ex. marking duplicates, which doesn't require reference file (@bsegall02)
- Continue trying to solve contigs issue for Mutect2. Ideas are renaming them in the fasta file, convertinig BAM to SAM, gzipping, and editing that, or looking for more reference files that have matching names (@besgall02, @jennifertramsu)
- Look through the github of Albert's paper, see if we can steal* parts of  their scripts to run the whole pipeline or parts.

## Thursday, June 1, 2023
Updates:
- Mutect2 should have a built-in function for mutation fractions, allele frequency
- Mutect2 is running without errors; except, there's no output
- Python SHAD drafts to be completed by Friday, test run on Tuesday
- The BAM -> Harmonized BAM step may be unnecessary. Depends on quality of data

Problems/Tasks:
- Figure out how to get the read count/fraction in Mutect2 output (@henry @olivialopardo)
- Solve output problem for Mutect2 (@bsegall02 @jennifertramsu)
- Look into BAM -> Harmonized BAM step, whether it's necessary and what tools to use (@pilar)
- Finish draft SHAD presentation/script, small edits in the activity and tutorial(@owen @bsegall02)
- Will start discussion on gRNA design: what factors to prioritize and what libraries can figure out those things (@albert @bsegall02 @jonas @henry)

## Monday, May 29, 2023
Problems:
- IndelRealigner not available in GATK > v3.6 -> look into alternatives (@olivialopardo)
- Calling Mutect2 on BAM files -> get GATK working then call Mutect2 (@bsegall02 @jennifertramsu)
- How to count # of reeds to get mutation frequency -> either manually or with a tool (@henry @olivialopardo)
    - Does Mutect2 already do this or do we need another tool?
