## Monday, May 29, 2023
Problems:
- IndelRealigner not available in GATK > v3.6 -> look into alternatives (@olivialopardo)
- Calling Mutect2 on BAM files -> get GATK working then call Mutect2 (@bsegall02 @jennifertramsu)
- How to count # of reeds to get mutation frequency -> either manually or with a tool (@henry @olivialopardo)
    - Does Mutect2 already do this or do we need another tool?

## Thursday, June 1, 2023
Updates:
- Mutect2 should have a built-in function for mutation fractions, allele frequency
- Mutect2 is running without errors; except, there's no output
- Python SHAD drafts to be completed by Friday, test run on Tuesday
- The BAM -> Harmonized BAM step may be unnecessary. Depends on quality of data

Problems/Assignments:
- Figure out how to get the read count/fraction in Mutect2 output (@henry @olivialopardo)
- Solve output problem for Mutect2 (@bsegall02 @jennifertramsu)
- Look into BAM -> Harmonized BAM step, whether it's necessary and what tools to use (@pilar)
- Finish draft SHAD presentation/script, small edits in the activity and tutorial(@owen @bsegall02)
- Will start discussion on gRNA design: what factors to prioritize and what libraries can figure out those things (@albert @bsegall02 @jonas @henry)
