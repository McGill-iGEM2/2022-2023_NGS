#!/bin/bash
bwa index GRCh38ref.gz
bwa mem GRCh38ref.gz SRR12331326_GL000219.1.fastq > aligned.sam