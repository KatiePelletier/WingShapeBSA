#!/bin/bash

set -e
set -u
set -o pipefail

files=$1

#a check to make sure this is the files you think it is
echo "the files are" $files

outname=$2

genome=/home/katie/flyGenome/dmel_r6.23/bwa/dmel-all-chromosome-r6.23.fasta


#the first step is to create the mpileup file and then pipes this to call variants 

samtools mpileup -B -Q 20 -f $genome $files > ${outname}.mpileup

