#! /bin/bash

#Variable for project:
project_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge

#Path to Picard
picdir=/usr/local/picard-tools/picard.jar


files=(${project_dir}/*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`

java -jar ${pic} AddOrReplaceReadGroups I=${final}/${base}.bam \
  O=${final}/${base}_RG.bam \
  RGID=L001_L002 \
  RGLB=library1 \
  RGPL=illumina \
  RGPU=None \
  RGSM=${base}

done

