#!/bin/bash
#Only works if all lanes are L001/L002

run1=
run2=
bam_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/
merged_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge
files=(${bam_dir}/${run1}/*_L001_aligned_pe.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_aligned_pe.bam
	samtools merge ${merged_dir}/${base}_merged_aligned_PE.bam \
	${bam_dir}/${run1}/${base}_L001_aligned_pe.bam \
	${bam_dir}/${run1}/${base}_L002_aligned_pe.bam \
	${bam_dir}/${run2}/${base}_L001_aligned_pe.bam \
	${bam_dir}/${run2}/${base}_L002_aligned_pe.bam \
done





