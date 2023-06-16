#!/bin/bash
ref=$2
sample=$3

bwa index $ref
bwa mem $ref $sample > aligned.sam