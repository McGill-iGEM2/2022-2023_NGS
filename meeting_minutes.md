## Thursday, June 22, 2023
Questions:
- Optical duplicates vs amplification duplicates
    - Optical duplicates not marked by `picard` because of RG error, regular duplicates are being detected and removed
    - Optical are due to instrument error. Not sure why they aren't removed normally, they're probably not significant
- Two files produced by `fasterq-dump` -> paired-end reads, can indicate in `bwa mem` that we are inputting paired-end reads. Can ignore technical reads
    - May be aligned or "interleaved," will have to check the first couple lines of the fastq files to see if the reads match
    - In cases of a third file, will have to be aligned separately on its own.
    - ![Image instructions for bwa](images/bwa_paired_ends.png)
- Why are we finding % of aligned reads?
    - Sanity check, making sure the data is good and we can proceed to next steps
    - Good % is context-dependent - when aligning to genome, ~85% is enough. Otherwise, would have to figure out what went wrong with collection, or something else
 - What info do we need for read groups?
    - [Listed here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups), but the naming is probably arbitrary and can work with whatever we call it.
    - Can add either with [picard AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-) or during the alignment with bwa ![Image for bwa flag](images/bwa_RG_flag.png).
  
Things to consider for future (and current) steps:
- Does our variant calling include masking of the germline variants?
- Make sure to add a step comparing to dbSNP to check for common mutations
- If normal data isn't great, can run tumor only and reference against common germline variants
- Indel realignment probably necessary for quality control still

## Monday, June 19, 2023
Assignments:
- Finding % of aligned reads (@elya)
- Getting raw fastq, aligning to reference genome and bringing them through preprocessing/mutect (@jennifertramsu, @bsegall02)
- Script to get expression level of different genes (@albert)

## Thursday, June 15, 2023
gRNA things to consider:
- Cas-9 DNA melting temperature
- mRNA structure
- expression level
- Whether the mutation kills the gene expression
- Frequency of those with mutation vs unmutated (and don't want it to bind to unmutated)
- Off targets are rarely an issue, usually a problem like paralogs since often things similar at the protein level would havae likely diverged at the DNA level
    - Can also check quickly by just blasting or using BWA, bowtie to search for similarities to the chosen gRNA
- More data on TCGA compared to SRA, much better for determing what genes to target vs SRA just to see the whole pipeline working from the start
    - [MAF files](https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22simple%20nucleotide%20variation%22%5D%7D%7D%5D%7D)
    -  Should also look to figure out what driver mutations are  known, try to get a fairly comprehensive list of known ones since that would influence the decision of what gRNA to use, whether the target is a driver mutation
        - [Software for finding driver mutations](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00830-0)
- Some  of the factors are binary, but with those that are not, it's hard to know what weight to give to each of them.. May just have to guess if there isn't enough time to test a bunch of gRNAs in the lab

gRNA how to do it:
- Some packages for cas-9 are based on melting temperature, which should still work for our purposes
- folding, can use RNA fold, mostly 2D folding though: usually pretty good though
- Tools for detecting NMD (Nonsense mediated decay)
    - Found a tool, [NMD classifier](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0174798). Requires a transcript so we would want to feed it every transcript that ref-seq predicts
- Targetting expression level - would be good to target maybe top 50% because it tends to be logorithmic
    - Need to start working on this, shouldn't be too bad. Looking through MAF files Jasmin send and writing a script to get the fraction of samaples where the gene is expressed
 
## Monday, June 12, 2023
Updates:
- Pilar working on the BWA part of the pipeline
- Elya has alignment working starting with a FASTQ file 
    - GCh38 as reference genome from NCBI, results in a SAM equivalent to the downloaded file from the original SRA data site
    - Still does not have the % of aligned reads, further steps to figure out how to get that
    - Can copy paste read in the UCSC Genome Browser to benchmark the read and see how it maps
- Looks like we can't store in compute canada and run locally, but shouldn't be an issue to run job on there. Jasmin's students have rarely had to wait for resources to open up.

Tasks:
- Will create a new markdown file to log intructions when trying to set up and use the programs (@bsegall02 @jennifertramsu)
- Continue work for Mutect2, can use albert's paper as a reference for commands (@pilar @bsegall02 @jennifertramsu)
- % of aligned reads, start logging process on github (@elya)
- Everybody apply for compute canada account (@everybody)
- Get wet lab to make a list of tools that would be useful to have for us to start working on

## Thursday, June 8, 2023
Tasks:
- Running the [pipeline](https://github.com/bioinform/somaticseq)
    - Includes a manual on how to download and use SomaticSeq --> get this working
- From the manual
    1) (Alignment and Cleaning) Run `makeAlignmentScripts.py` with the necessary flags to trim, align, and mark duplicates (@pilar)
        - Input: fastq (tumour, normal), fasta (reference)
        - Output: merged fastq file, trimmed bam file
    2) (Variant Calling - Tumour-Normal Paired Mode) Run `makeSomaticScripts.py` to create scripts for MuText2, SomaticSniper, VarDict, etc. Use the flag `--run-workflow` to execute the scripts in parallel once the scripts are created (@jennifertramsu @bsegall02)
        - Input: bam (tumour, normal)
        - Output: ??
    - We will probably want to run MuTect2 individually. Will also need to find argument to pass to MuTect2 that calculates allele frequency.
- Compute Canada -> storing data on Compute Canada and finding out how to access the data locally (@bsegall02)

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
- Look through the github of Albert's paper, see if we can steal* their scripts to run the whole pipeline or individual parts (@pilar, @bsegall02)

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
