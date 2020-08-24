#!/bin/bash

#Variable for project:
#project_dir=/home/sarahm/cvl/storage

#Path to input directory
final_bam=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge/

#Path to output directory
gatk_dir=/2/scratch/Katie/WingShapeBSA/temp/gatk/index/

#Variable for reference genome (non-zipped)
#index_dir=/home/sarahm/cvl/index_dir
ref_genome=/home/katie/flyGenome/dmel_r6.23/gatk/dmel-all-chromosome-r6.23.fasta

#Path to GATK
gatk=/usr/local/gatk/GenomeAnalysisTK.jar


files=(${final_bam}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`

java -Xmx32g -jar ${gatk} -I ${final_bam}/${base}_RG.bam -R ${ref_genome} \
  -T RealignerTargetCreator \
  -o ${gatk_dir}/${base}.intervals

done