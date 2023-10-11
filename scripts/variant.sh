#!/bin/bash

while getopts "r:n:t:o:g:p:d:" option; do
    case $option in
        r) ref=$OPTARG;;
        n) normal=$OPTARG;;
        t) tumour=$OPTARG;;
        o) out=$OPTARG;;
        g) germline=$OPTARG
            if [ ! -f $germline.tbi ]; 
            then
                echo -e "\nIndexing feature files...\n"
                gatk IndexFeatureFile -I $germline
            fi
        ;;
        p) panel=$OPTARG
            if [ ! -f $panel.tbi ];
            then
                echo -e "\nIndexing feature files...\n"
                gatk IndexFeatureFile -I $panel
            fi
        ;;
        d) tmp=$OPTARG;;
    esac
done

# job on CC - allocate more memory
# by default, java has max size to its heap - no matter how much memory on node, if it exceeds it stops
# if calling through java, option available to increase size of heap (google)
# help desk CC -> good resource
# gatk way to go
#   most commonly used

# comparing multiple variant callers
#   time consuming, gatk may have some FP
#   intersecting predictions could remove more FP but not best use of time
# FP clinical setting
#   group of researches @ mcgill -> working on this, C3G
#   can contact them
# detecting mutations?
# accurate? which are pathogenic? harder question to validate
#   variants detected by gatk come with confidence values
#   beyond that, people will visualize alignments with IGV to convince themselves that calls look reasonable (not for every variant, for those you are interested in)

# Variant calling
gatk Mutect2 \
    --java-options "-Xmx32G" \
    -R $ref \
    -I $normal \
    -I $tumour \
    -normal normal \
    --germline-resource $germline \
    --panel-of-normals $panel \
    -O $out.vcf \
    --tmp-dir $tmp

echo -e "\nFiltering vcf...\n"

gatk FilterMutectCalls \
    -R $ref \
    -V $out.vcf \
    -O $out.filtered.vcf

bgzip -c $out.filtered.vcf > $out.filtered.vcf.gz
tabix -p vcf $out.filtered.vcf.gz

# oct 2
# statistical model..? like sensitivity or precision, is this needed
# using software like VEP to annotate pathogenicity
# driver mutations, how do we find them? based on freuqency

# challenging question of identifying driver mutations
# - programs that take list of mutations and assign them severity score
# - SIFT, polyphen (both old)
# - determine whetehr mutation in coding region would change an amino acid in protein or silent
# - whether AA change is likely to have consequence on protein function -> score or ranked list (!)
# - not super accurate
# - don't know interested in certain cancer, don't know which genes are associated with the cancer
#     - how to navigate that?
#     - talk to C3G!!
#     - unsure, maybe there are programs out there
#     - make a manual list of known 
#     - if many pairs, can do analysis separete for each pair, then is same mutation found across samples or diff mnutatoin in same gene
#         - that suggests this is driver mmutation
#     - with only one pair, every mutatoin in  mprinciple could be cause in cancer
#         - one thing you can do is 
#             - suppose if you find 20 mutations with significent severity -> look up those 20 genes see whta's known about them in role in cancer and try to convince yourself that it's a driver mutation
#             - this is manual and time consuming, error prone -> esp if no ecisting data in ltierature talks about it
# - it just looks at general consequence
# - obvious, itnroduced stop codon in middle of sequence -> totally kill function of protein
# - substitution may have cons or may not

# testing diff programs
# - varainta nnotation is where can benefit from testing different programs!
# - many differences between software
#     - take intersection to narrow down list is good idea

# oct 3, 2023
# - driver mutations
#     - look for mutations that are reoccurring across multiple tumours, not just random passneger mutation
# - challenge
#     - unlikely that two tumours will have the exact mutation at the same position
#     - diff tumours may have diff tumours int he same gene (all disable gene but in the diff ways)
#     - more important o look at disabling of a aprticular gene
#     - not his area
#     - tehre are tools for this analysis
#         - predicting effect on gene tool
#         - another takes input as mutation found in multiple tumours and then identify genes hit by mutations across these mutations
#         - identifying recurrent mutations or genes that are ercurrently being disabled by these mutations
# - statistical model
#     - if developing those tools, important to have statistics
#     - if a user, then less of a burden to demonstrate that
#     - can't take output for granted, can't assume that it's correct
#         - certain tasks like identifying recurrent mutations
#         - suppose you look at 10 pairs, identify mutation in gene in 3/10 tumours
#             - is this indicative that gene is involved? or just random chance?
#                 - this is the question that needs to be addressed
#                 - some kind of p-value
#                 - statistics that he would expect a tool would have
#         - can't prove beyond all doubts with just a computational analysis
#             - best result to hope for is that all 10 tumours all have mutation in a given gene
#             - often people will follow up discovery with lab work
#                 - what happens if you take cell with this mutation (add or remove)
# - finding off-target effects when designing gRNAS
#     - input on best way to go about that
#     - we were thinking of BLAST or BLAT to find closely matching sequences
#         - not his field :')
#         - would expect that there are other tools for off-target effects for gRNAs
#             - can look for these tools, can save downstream tools
#     - existing tools we've looked at are for the Cas9 system but we are using Cas7-11
#         - a longer sequence, 31 nt and does not have a PAM sequence
#         - main thing is that our sysetm cleaves RNA, not DNA
#             - If no explicit tool, then might have to do it with our own code
#             - Not sure what level of complementary required to trigger
#     - main diff blast blat
#         - blat is much more faster but will only find matches for grna that are highly similar
#         - for 10 grna, both tools are fast enough (runtime)
#         - bigger question is what level of complementarity is needed between grna and mrna to trigger
#         - might be published somewhere
#         - if you run blast with given grna, will find hundreds of matches
#             - some with high complemetnarity, some with lower
#             - but doesn't mean all those alignments would all generate off-target effects
#                 - maybe only the most complementary would work
#         - we know about seed region, but not much else in terms if complementarity
#             - might be useful to assess the off-target effects experimentally
#             - but PI might be in better place to make suggestions
#             - once we've designed gRNA, assess specificity experimentally
#     - what can we do computationally
#         - how to design the grna to target which portion of the target gene
#         - in such a way that they have little predicted off-target effects as possible
