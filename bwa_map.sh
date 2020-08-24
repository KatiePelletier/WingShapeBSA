#!/bin/bash

#project directory (to keep each run separate)
run=PhO_controls

# directory of processed sequences with trimmomatic 
trim_dir=/home/katie/dacsProject/Test

# variable for the reference genome
refGenome=/home/katie/flyGenome/dmel_r6.23/bwa/dmel-all-chromosome-r6.23.fasta 

# make output directory from mapping outputs
output=//home/katie/dacsProject/Test

# make BWA directory path
bwa_dir=/usr/local/bwa/0.7.8

cd ${bwa_dir}

#list all files to be read (this selects the left end from each PE pair)
files=(${trim_dir}/*_R1_PE.fastq.gz)

#echo ${files[@]}

#For loop over every file
for file in ${files[@]} 
do
name=${file}
base=`basename ${name} _Head.fastq`

#echo base is
#echo ${base}
#echo name is 
#echo ${name}

bwa mem -t 8 -M ${refGenome} ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz > ${output}/${base}_bwa_PE.SAM
done


