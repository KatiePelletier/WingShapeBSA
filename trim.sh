#!/bin/bash
# Project directory variable:
# project_dir=/home/katie/dachsProject/WingShapeBSA/20141205_A_DNASeq_PE

# make variable for trimmomatic program location
trim=/usr/local/trimmomatic/trimmomatic-0.36.jar

# make input directory for raw reads
raw_dir=/home/katie/dachsProject/WingShapeBSA/20141205_A_DNASeq_PE

# make output directory from trimmomatic outputs
trim_dir=/2/scratch/Katie/WingShapeBSA/trimmomatic/20141205_A_DNASeq_PE

# make path to adapter sequences (to be used with ILLUMINACLIP)
adapter=/home/katie/adapterSequences/TruSeq3-PE.fa

#list all files to be read (all raw data)
files=(${raw_dir}/*_R1_001.fastq.gz)

#For loop over every file
for file in ${files[@]} 
do
name=${file}
base=`basename ${name}_R1_001.fastq.gz`

java -jar ${trim} PE -phred33 \
  -trimlog ${trim_dir}/trimlog.txt \
  ${raw_dir}/${base}_R1_001.fastq.gz \
  ${raw_dir}/${base}_R2_001.fastq.gz \
  ${trim_dir}/${base}_R1_PE.fastq.gz \
  ${trim_dir}/${base}_R1_SE.fastq.gz \
  ${trim_dir}/${base}_R2_PE.fastq.gz \
  ${trim_dir}/${base}_R2_SE.fastq.gz \
  ILLUMINACLIP:${adapter} \
  LEADING:3 \
  TRAILING:3 \
  MAXINFO:40:0.5 \
  MINLEN:36
  
done