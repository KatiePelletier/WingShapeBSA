library(poolSeq)
library(ACER)

args <- commandArgs(trailingOnly = TRUE)

#read in the file
reps <- c(1:12)
gen <- rep(0,12)
sync <- read.sync("../mpileup/wild_masked_coregenome_spaces.sync",
                  gen=gen, repl=reps,
                  polarization = "minor",
                  keepOnlyBiallelic = TRUE)


sync <- read.sync("../mpileup/neur_new_coregenome.sync", 
                  gen=gen, repl=reps, 
                  polarization = "minor", 
                  keepOnlyBiallelic = TRUE)

pops <-c('CMO.L', 'CMO.R',
            'FVW12.L', 'FVW12.R',
            'FVW14.L', 'FVW14.R',
            'PHO.L', 'PHO.R')
pops <-c('CMO.L',  'CMO.R', 'CMO.C',
         'FVW12.L',  'FVW12.R', 'FVW12.C',
         'FVW14.L', 'FVW14.R', 'FVW14.C',
         'PHO.C', 'PHO.L', 'PHO.R')
       

af.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 12)
colnames(af.mat) <- pops

for (i in 1:ncol(af.mat)){
  tempdat <- af(sync, repl = i, gen = 0)
  af.mat[,i] <- as.matrix(tempdat)
}

af.mat <- na.omit(af.mat)
head(af.mat)
dim(af.mat)

af.mat2 <- af.mat[,c(1,2,4,5,7,8,11,12)]
dim(af.mat2)

no.pho.af <- af.mat2[,1:6]
head(no.pho.af)
    
#now to make a coverage one. 

cov.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 12)
colnames(cov.mat) <- pops

for (i in 1:ncol(cov.mat)){
  tempdat <- coverage(sync, repl = i, gen = 0)
  cov.mat[,i] <- as.matrix(tempdat)
}

crap <- data.frame(cov.mat, sync@alleles[,1:2])
crap[crap==0] <- NA
crap2 <- na.omit(crap)
location <- crap2[,13:14]

cov.mat[cov.mat == 0] <- NA
cov.mat <- na.omit(cov.mat)

dim(cov.mat)

head(cov.mat)

cov.mat2 <- cov.mat[,c(1,2,4,5,7,8,11,12)]
dim(cov.mat2)

no.pho.cov <- cov.mat[,c(1,2,4,5,7,8)]
head(no.pho.cov)

#Now I want to estimate Ne to use below. 
# ne <- estimateNe(p0=af.mat2[,"CMO.L"], pt=af.mat2[,"CMO.R"], 
#            cov0=cov.mat2[,"CMO.L"], covt=cov.mat2[,"CMO.R"], 
#            t=0, method = "P.planII", poolSize=c(150, 150))

#Creating the vars for the CMH test 
rep<-c(1,1,2,2,3,3,4,4) #Number of replicates
Ne<-rep(10e6, 4) #Using a pretty general estimate. Should this be estimated from the unselected pools?
tp<-rep(0,4) #Generations of evolution for each sample
ps<-rep(75, 8) #Pool size

#for 3 populations (leaving pho out) 
rep2<-c(1,1,2,2,3,3) #Number of replicates
Ne2<-rep(10e6, 3) #Using a pretty general estimate. Should this be estimated from the unselected pools?
tp2<-rep(0,3) #Generations of evolution for each sample
ps2<-rep(75, 6) #Pool size

#cmh test 
#I think I need to include the random pools as a gen 0?

pval <- adapted.cmh.test(freq=af.mat2, coverage=cov.mat2, 
                         Ne=Ne, gen=tp, repl=rep, poolSize=ps)

pval.no.pho <- adapted.cmh.test(freq=no.pho.af, coverage=no.pho.cov, 
                         Ne=Ne2, gen=tp2, repl=rep2, poolSize=ps2)


# Warning messages:
#   1: In adapted.cmh.test(freq = af.mat, coverage = cov.mat, Ne = Ne,  :
#       Ne value(s) which are not integer are converted to integer
#I took care of this (warning 2) by setting left to gen 0 and right to gen 1 
#   2: In adapted.cmh.test(freq = af.mat, coverage = cov.mat, Ne = Ne,  :
#       Value of 'Ne' will be ignored because no random genetic drift is assumed.
#   3: In adapted.cmh.test(freq = af.mat, coverage = cov.mat, Ne = Ne,  :
#       The counts that equal 0 or equal the coverage in all replicates are changed to 1 or to coverage-1 respectively.


#these are all 1. So everything goes away when we account for drift?
#padj <- p.adjust(pval, "fdr")

padj <- p.adjust(pval, "bonferroni")

padj <- p.adjust(pval, "fdr")

padj.no.pho <-  p.adjust(pval.no.pho, "fdr")


#This doesn't line up with the locations from the original file because we droped those NA lines. 
afdat <- cbind(af.mat, pval, padj)

afdat2 <- cbind(afdat, location)

dat.no.pho <- cbind(no.pho.af, pval.no.pho, padj.no.pho, location)

#3 sites.... seems low. but cool (with 1 gen of selection)
#15 sites for 0 gen. 
library(dplyr)
sig <- filter(afdat2, padj <= 0.05)

sig.no.pho <- filter(dat.no.pho, padj.no.pho <= 0.05)


write.csv(afdat2, "ds_wild_ACER_zeroGen_fdr.csv", quote = FALSE, row.names = FALSE)
# write.csv(afdat2, "ds_wild_ACER_significantSites_zeroGen.csv", quote = FALSE, row.names = FALSE)

write.csv(dat.no.pho, "ds_wild_ACER_zeroGen_fdr_noPho.csv", quote = FALSE, row.names = FALSE)


write.table(sig, "ds_wild_ACER_significantSites_zeroGen_fdr.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sig.no.pho, "ds_wild_ACER_significantSites_zeroGen_fdr_noPho.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
