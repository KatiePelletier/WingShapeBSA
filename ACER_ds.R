library(poolSeq)
library(ACER)

#args <- commandArgs(trailingOnly = TRUE)

#read in the file
reps <- c(1:8)
gen <- rep(1,8)
sync <- read.sync(file="../mpileup/LRpoolsOnly_wild_masked_coregenome_75subsample.sync", 
                  gen=gen, repl=reps, 
                  polarization = "minor", 
                  keepOnlyBiallelic = TRUE)

head(sync@alleles[,1:2])

pops <-c('CMO.L', 'CMO.R', 
             'FVW12.L', 'FVW12.R',
             'FVW14.L', 'FVW14.R',
             'PHO.L', 'PHO.R')
# pops <-c('CMO.L', 'CMO.C', 'CMO.R', 
#                 'FVW12.L', 'FVW12.C', 'FVW12.R',
#                 'FVW14.L', 'FVW14.R',
#                 'PHO.L', 'PHO.R')
       

af.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 8)
colnames(af.mat) <- pops

for (i in 1:ncol(af.mat)){
  tempdat <- af(sync, repl = i, gen = 1)
  af.mat[,i] <- as.matrix(tempdat)
}

af.mat <- na.omit(af.mat)
head(af.mat)
dim(af.mat)

dim(sync@alleles)

    
#now to make a coverage one. 

cov.mat <- matrix(NA, nrow = nrow(sync@alleles), ncol = 8)
# cov.mat[,1:2] <- sync@alleles[1,]
colnames(cov.mat) <- pops

for (i in 1:ncol(cov.mat)){
  tempdat <- coverage(sync, repl = i, gen = 1)
  cov.mat[,i] <- as.matrix(tempdat)
}

crap <- data.frame(cov.mat, sync@alleles[,1:2])
crap[crap==0] <- NA
crap2 <- na.omit(crap)
location <- crap2[,9:10]

cov.mat[cov.mat==0] <- NA
cov.mat <- na.omit(cov.mat)

dim(cov.mat)

head(cov.mat)


#Now I want to estimate Ne to use below. 
ne <- estimateNe(p0=af.mat[,"CMO.L"], pt=af.mat[,"CMO.R"], 
           cov0=cov.mat[,"CMO.L"], covt=cov.mat[,"CMO.R"], 
           t=1, method = "P.planII", poolSize=c(75, 75))

#Creating the vars for the CMH test 
rep<-c(1,1,2,2,3,3,4,4) #Number of replicates
Ne<-rep(ne, 4) #Using a pretty general estimate. Should this be estimated from the unselected pools?
tp<-rep(rep(c(0,1)),4) #Generations of evolution for each sample
ps<-rep(75, 8) #Pool size

#cmh test 
#I think I need to include the random pools as a gen 0?

pval <- adapted.cmh.test(freq=af.mat, coverage=cov.mat, 
                         Ne=Ne, gen=tp, repl=rep, poolSize=ps)


# Warning messages:
#   1: In adapted.cmh.test(freq = af.mat, coverage = cov.mat, Ne = Ne,  :
#       Ne value(s) which are not integer are converted to integer
#I took care of this (warning 2) by setting left to gen 0 and right to gen 1 
#   2: In adapted.cmh.test(freq = af.mat, coverage = cov.mat, Ne = Ne,  :
#       Value of 'Ne' will be ignored because no random genetic drift is assumed.
#   3: In adapted.cmh.test(freq = af.mat, coverage = cov.mat, Ne = Ne,  :
#       The counts that equal 0 or equal the coverage in all replicates are changed to 1 or to coverage-1 respectively.


#these are all 1. So everything goes away when we account for drift?
padj <- p.adjust(pval, "fdr")
#This doesn't line up with the locations from the original file because we droped those NA lines. 
afdat <- cbind(location, af.mat, pval, padj)

write.csv(afdat, "subsample75cov_correctPoolSize_dsWild_ACER.csv")

#write.csv("neur_ACER.csv")