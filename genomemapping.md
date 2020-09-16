# Mapping and preparing sequences for analysis 

Trimming adapter sequences from files.  I used a very large adapter file because subsets of this were leaving adapter contamination in my sequences 
Set minimum length to 36 bases. 

[trimming script](https://github.com/KatiePelletier/WingShapeBSA/blob/master/trim.sh)


Mapped with BWA to the Drosophila melanogaster v6.23 genome

[mapping script](https://github.com/KatiePelletier/WingShapeBSA/blob/master/bwa_map.sh)

I then merged reads from the same biological sample. Samples were sequenced twice to increase read depth. 

[merging biological replicates](https://github.com/KatiePelletier/WingShapeBSA/blob/master/merge.sh)

Following mapping, I converted SAM to BAM files (compressing) using SAMtools. This step also sorts the file.
Will also filter for mapping quality less than 20. 

[convert to BAM and sort](https://github.com/KatiePelletier/WingShapeBSA/blob/master/samTObam.sh)

# Realigning around indels using GATK 

First I have to add read group information. In the future, I should do this earlier in the pipeline 

[add read groups](https://github.com/KatiePelletier/WingShapeBSA/blob/master/addreplacegroups.sh)

Indexing the files 

[GATK index](https://github.com/KatiePelletier/WingShapeBSA/blob/master/gatkindex.sh)

This step identifies the indels in the files 

[find intrivals](https://github.com/KatiePelletier/WingShapeBSA/blob/master/gatkintravals.sh)

Realignment step

[GATK reallign](https://github.com/KatiePelletier/WingShapeBSA/blob/master/gatkalign.txt) 

Finally, I removed duplicate reads using Picard 

[remove duplicates](https://github.com/KatiePelletier/WingShapeBSA/blob/master/dedup.sh)

At this point, I can calculate the read depth for each sample

[calculate read depth](https://github.com/KatiePelletier/WingShapeBSA/blob/master/calcdepth.sh)

I also merged treatments together at this point for the artificial selection experiment 
[merging ee lineages](https://github.com/KatiePelletier/WingShapeBSA/blob/master/ee_merging.txt)

At this point, I was able to move on to creating pileup/mpileup files and population genetics analysis. 








