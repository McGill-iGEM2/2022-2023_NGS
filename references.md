# References
## Data
Reference genome: https://www.ncbi.nlm.nih.gov/genome/guide/human/ 

Tumour: https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR12331371&display=alignment

Normal: https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR12331408&display=download

dbSNP (Mutect 2): https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/
## Software
Fasterq-dump
  https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump/a58c2a07022dc7b0b043529ee01879b5381fe862

RNA Structure Prediction Tool (mfold)
  http://www.unafold.org/mfold/software/download-mfold.php

GATK Workshop GDrive
  https://drive.google.com/drive/folders/1y7q0gJ-ohNDhKG85UTRTwW1Jkq4HJ5M3

NGS Pipeline Tutorials
  https://hbctraining.github.io/variant_analysis/lessons/06_alignment_file_processing.html
  
## Readings
A software tool to help prioritize potential driver mutations (since these are more likely to be specific and universal to the tumor).
  https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00830-0

A relatively recent and well (enough) cited predictor of non-sense-mediated decay 
(to avoid targeting mutations that result in transcript degradation. remember you only need to run this 
for transcript that have premature stop codons, either because the mutation introduces a stop codon or a frame-shift).
  https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0174798
  
Best practices paper with SomaticSeq pipeline.
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8506910/
  
You will need to guess whether the mRNA that you want to target will be expressed (at a sufficient level) in the tumor. 
I suggest that for each gene, you count the fraction of tumor samples where the gene is unexpressed or lowly expressed 
(say in the bottom 50% of genes) and avoid targeting genes with a high fraction. 
After that, you will also need to predict whether the mutation(s) in the gene would cause it to be unexpressed 
(mostly through non-sense-mediated decay), only for those with frameshifts or pre-mature stop codons.
  https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22content%22%3A%7B%22field%22%3A%22files.cases.primary_site%22%2C%22value%22%3A%5B%22pancreas%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22transcriptome%20profiling%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Gene%20Expression%20Quantification%22%5D%7D%7D%5D%7D
  
Using the sra toolkit.
  https://erilu.github.io/python-fastq-downloader/#downloading-fastq-files-using-the-sra-toolkit

For the guide designs, I suggest you start straight from the TCGA maf files.
  https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22simple%20nucleotide%20variation%22%5D%7D%7D%5D%7D
 
This is the link to all the whole exome sequencing (WES for future reference) data from pancreatic cancer and the matched normals.
  https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP273725&o=acc_s%3Aa
  
  To find a sample that has its matched normal, look for runs with the same isolate (i.e. the 5LTS isolate has 5LTS_Met1, metastatic sample, 
  and 5LTS_N, its matched normal sample)

This paper seems to suggest Mutect does keep track of mutation fractions (i.e. total # of reads w a certain mutation). 
From the paper: “Variant detection: Variants in the tumor are identified by analyzing the data at each site under two alternative models: 
(i) a reference model, M0, which assumes there is no variant at the site and any observed non-reference bases are due to random sequencing errors; 
and (ii) a variant model, Mmf, which assumes the site contains a true variant allele m at allele fraction f in addition to sequencing errors. 
The allele fraction f is unknown and is estimated as the fraction of tumor reads that support m.”
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3833702/

Other papers:
  https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00830-0
