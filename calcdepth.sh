#! /bin/bash

files=(*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _merged_aligned_PE_realigned_rmd.bam`

samtools depth ${name} > 
