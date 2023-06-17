# Errors encountered and their solutions
## BWA Mem

## GATK

## Mutect2

## 2) Cleaning

### Section: 2c) Marking and removing duplicates

- Problem: `picard MarkDuplicates` gives the following warning:
    
        AbstractOpticalDuplicateFinderCommandLineProgram A field field parsed out of a read name was expected to contain an integer and did not. Read name: SRR12331371.100000157. Cause: String 'SRR12331371.100000157' did not start with a parsable number.

- Solution: 

### Section: 3b) Variant Calling

- Problem: `gatk Mutect2` throws the following exception:

        java.lang.IllegalArgumentException: samples cannot be empty

- Solution: The BAM files are missing read group headers.

        picard AddOrReplaceReadGroups I=<<input.bam>> O=<<output.bam>> RGID=<<read_group_id>> RGLB=<<read_group_library>> RGPL=<<read_group_platform>> RGPU=<<read_group_platform_unit>> RGSM=<<read_group_sample>>

- Problem: `gatk Mutect2` throws the following exception:

        java.lang.IllegalArgumentException: Reference name for '91553841' not found in sequence dictionary.

- Solution: 