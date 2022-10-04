#! /bin/bash 


#Variable to put indexed files in 
index_dir=/2/scratch/Katie/WingShapeBSA/temp/gatk/index

#Path to input directory
input=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge

files=(${input}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`

samtools index ${input}/${base}_RG.bam 

done

