#! /bin/bash

mapped_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge/gatkindel/
outdir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge/gatkindel/wildpop/
picdir=/usr/local/picard-tools/picard.jar

files=(${mapped_dir}/C*.bam)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`
#echo ${name}
java -Xmx2g -jar ${picdir} MarkDuplicates I=${mapped_dir}/${base}.bam O=${outdir}/${base}_rmd.bam M=${outdir}/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
done



