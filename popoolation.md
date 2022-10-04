# Make mpileup 

Use the final bam files, using merged files for the artifical selection and unmerged files for wild 

[make mpileup script](https://github.com/KatiePelletier/WingShapeBSA/blob/master/GenomicMapping/make_mpileup.sh)

# Masking Repeats

Used the repetitive regions identified by Repeat Masker and the script from Popoolation2

```
perl /home/katie/bin/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --gtf /home/katie/flyGenome/dmel_r6.23/dmel-all-chromosome-r6.23.fasta.out.gff --input ee_merge.mpileup --output ee_merge_repeatmask.mpileup &
perl /home/katie/bin/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --gtf /home/katie/flyGenome/dmel_r6.23/dmel-all-chromosome-r6.23.fasta.out.gff --input wild_nomerge.mpileup --output wild_nomerge_repeatmask.mpileup
```

# Masking InDels 

First I identified indels using the popoolation2 script, used a 5bp window around indels  

```
perl /home/katie/bin/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input ee_merge.mpileup --output ee_merge_indels.gtf --indel-window 5 &
perl /home/katie/bin/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input wild_nomerge.mpileup --output wild_nomerge_indels.gtf --indel-window 5 
```

Then mask indels in mpileup file 

```
perl /home/katie/bin/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --gtf ee_merge_indels.gtf --input ee_merge_repeatmask.mpileup --output ee_merge_repeatmask_indelmask.mpileup &
perl /home/katie/bin/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --gtf wild_nomerge_indels.gtf --input wild_nomerge_repeatmask.mpileup --output wild_nomerge_repeatmask_indelmask.mpileup
```
# Make sync file 

```
java -ea -Xmx32g -jar /usr/local/popoolation/mpileup2sync.jar --input ee_merge_repeatmask_indelmask.mpileup --output ee_masked_merged.sync &
java -ea -Xmx32g -jar /usr/local/popoolation/mpileup2sync.jar --input wild_nomerge_repeatmask_indelmask.mpileup --output wild_masked_nomerge.sync
```

#Calculate Fst in sliding windows. 

```
perl /home/katie/bin/popoolation2_1201/fst-sliding.pl --input ee_masked_merged.sync --output ../population/ee_merged_5000windows.fst --suppress-noninformative --min-count 2 --min-coverage 10 --max-coverage 600 --min-covered-fraction 0.2 --window-size 5000 --step-size 5000 --pool-size 75 &
perl /home/katie/bin/popoolation2_1201/fst-sliding.pl --input wild_masked_nomerge.sync --output ../population/wild_nomerge_100windows.fst --suppress-noninformative --min-count 2 --min-coverage 10 --max-coverage 600 --min-covered-fraction 0.2 --window-size 100 --step-size 5000 --pool-size 75 
```

For the wild popualtions, I want to pull out only the comparisons I am interested in. 
```
awk '{print $1, $2, $3, $4, $5, $6, $36, $57, $71}' wild_nosync_100window.fst > wild_lvr_nosync_100window.fst &
awk '{print $1, $2, $3, $4, $5, $7, $37, $58, $69}' wild_nosync_100window.fst > wild_lvc_nosync_100window.fst &
awk '{print $1, $2, $3, $4, $5, $17, $44, $62, $70}' wild_nosync_100window.fst > wild_rvc_nosync_100window.fst 
```

For the artifical selection population the order of comparisons are: control vs down, control vs up, down vs up. 



