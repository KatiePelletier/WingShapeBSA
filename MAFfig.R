##Going to look at the allele F in ds and emc estimated from selected lines. 

#DO NOT DO
#library(pegas)
#ds.vcf <- read.vcf("../Data/parents2L.vcf.gz", which.loci = 240013:714983)

#library(vcfR)

#chr2L <- read.vcfR("../Data/parents2L.vcf.gz")
library(tidyr)
library(data.table)


#This is the dgrp vcf file but subsetted for the 30 lines used in this experiment. 
chr3L <- fread("../Data/parents3L.vcf")

crap <- chr3L[chr3L$POS >= 749400,]
emconly <- crap[crap$POS <= 753492,]

#Now I am going to extract the genotypes. A 1 == alt allele and 0 == ref allele. there are . == no data in here as well. 

#this is a really bad way to do this but it works for now. 

for(i in 10:39){
  emconly[[i]] <- gsub("1/1", "1", emconly[[i]])
}
for(i in 10:39){
  emconly[[i]] <- gsub("0/0", "0", emconly[[i]])
}
for(i in 10:39){
  emconly[[i]] <- gsub("./.", "NA", emconly[[i]])
}


for(i in 10:39){
  emconly[[i]] <- as.numeric(emconly[[i]])
}


#now to get allele frequency. 

emconly$aallele <- rowSums(emconly[,10:39], na.rm = TRUE)
#number of lines without data avalible. 
emconly$na.count <- rowSums(is.na(emconly[,10:39]))

emconly$freq <- emconly$aallele/(30 - emconly$na.count)

#This is the alternate allele compared to the refrence genome. 
#MAF is also weird because this is called with the whole data set and then I subsetted. 
hist(emconly$freq)

#Now I want to also calculate pi using:
# h = n/n-1(1-sum(pi^2)) where pi represents the freq of the ith allele at any site 
#then pi = sum(hj) where hj is the heterozygosity at the jth site 

#first, I only want the sites that are variant 
emc.var1 <- emconly[emconly$freq > 0,]
#there are also indels here? So I want only SNPs
emc.var <- emc.var1[grep("SNP", emc.var1$ID),]

#now to write out how I would get h at each site in peices. 

#30 lines in sample 
emc.var$n <- 30 - emc.var$na.count
  
emc.var$p1 <- emc.var$freq^2
emc.var$p2 <- (1 - emc.var$freq)^2

emc.var$h <- (emc.var$n/(emc.var$n-1))*(1- (emc.var$p1 +emc.var$p2))
hist(emc.var$h)

emclength <- 753492 - 749400 


#0.0051
emc.pi <- sum(emc.var$h)/emclength


#This is the dgrp vcf file but subsetted for the 30 lines used in this experiment. 
chr2L <- fread("../Data/parents2L.vcf")

crap <- chr2L[chr2L$POS >= 240013,]
dsonly <- crap[crap$POS <= 714983,]

#Now I am going to extract the genotypes. A 1 == alt allele and 0 == ref allele. there are . == no data in here as well. 

#this is a really bad way to do this but it works for now. 

for(i in 10:39){
  dsonly[[i]] <- gsub("1/1", "1", dsonly[[i]])
}
for(i in 10:39){
  dsonly[[i]] <- gsub("0/0", "0", dsonly[[i]])
}
for(i in 10:39){
  dsonly[[i]] <- gsub("./.", "NA", dsonly[[i]])
}


for(i in 10:39){
  dsonly[[i]] <- as.numeric(dsonly[[i]])
}

class(dsonly[[39]])
is.na(dsonly[[39]])

#now to get allele frequency. 

dsonly$aallele <- rowSums(dsonly[,10:39], na.rm = TRUE)
#number of lines without data avalible. 
dsonly$na.count <- rowSums(is.na(dsonly[,10:39]))

dsonly$freq <- dsonly$aallele/(30 - dsonly$na.count)

#This is the alternate allele compared to the refrence genome. 
#MAF is also weird because this is called with the whole data set and then I subsetted. 
hist(dsonly$freq)


#Now I want to also calculate pi using:
# h = n/n-1(1-sum(pi^2)) where pi represents the freq of the ith allele at any site 
#then pi = sum(hj) where hj is the heterozygosity at the jth site 

#first, I only want the sites that are variant 
ds.var1 <- dsonly[dsonly$freq > 0,]
#there are also indels here? So I want only SNPs
ds.var <- ds.var1[grep("SNP", ds.var1$ID),]

#now to write out how I would get h at each site in peices. 

#30 lines in sample 
ds.var$n <- 30 - ds.var$na.count

ds.var$p1 <- ds.var$freq^2
ds.var$p2 <- (1 - ds.var$freq)^2

ds.var$h <- (ds.var$n/(ds.var$n-1))*(1- (ds.var$p1 +ds.var$p2))
hist(ds.var$h)

dslength <- 714983 - 240013


#0.0039273
ds.pi <- sum(ds.var$h)/dslength
ds.pi

#Here, the frequency are on two diffrent axis because there are more alleles in ds than emc. but I don't think that matters?

#Going to try making this in ggplot. 

library(ggplot2)

ggplot(dsonly, aes(x = freq)) + 
  geom_histogram() + 
  xlim(0, 1) + 
  ylim(0, 3000)


ggplot(emconly, aes(x = freq)) + 
  geom_histogram() + 
  xlim(0, 1) + 
  ylim(0, 30)

#but I want to plot both together

emcfreq <- data.frame(emconly$freq)
emcfreq$gene <- "emc"
colnames(emcfreq)[1] <- "freq"

dsfreq <- data.frame(dsonly$freq)
dsfreq$gene <- "ds"
colnames(dsfreq)[1] <- "freq"

genefreq <- bind_rows(emcfreq, dsfreq)

png("../Figures/AltAlleleFreq_plottedTogether.png")
ggplot(genefreq, aes(x = freq, fill = gene)) + 
  geom_histogram() + 
  xlab("Alt Allele Freq") + 
  theme_classic()
dev.off()

library(cowplot)

dsplot <- ggplot(dsfreq, aes(x = freq)) + 
  geom_histogram() + 
  xlab("Alt Allele Freq") + 
  theme_classic() + 
  xlim(0, 1) + 
  ylim(0, 3000)


emcplot <- ggplot(emcfreq, aes(x = freq)) + 
  geom_histogram() + 
  xlab("Alt Allele Freq") + 
  theme_classic()+ 
  xlim(0, 1) + 
  ylim(0, 30)

allPlot <- plot_grid(dsplot, emcplot, labels = c("ds", "emc"), label_size = 10)

png("../Figures/AltAlleleFreq_plottedSeperate.png")
allPlot
dev.off()
