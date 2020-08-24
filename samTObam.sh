#!/bin/bash

#specify run folder 
run=

#Specify input directory 
sam_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/SAM_files/{run}/

#Specify output directory 
bam_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/{run}/


#Loop over each file to convert it 
files=(${sam_dir}/*.SAM)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .SAM`
samtools view -b -q -@5 ${sam_dir}/${base}.SAM | samtools sort - ${bam_dir}/${base}
done

