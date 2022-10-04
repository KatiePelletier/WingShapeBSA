#doing the log regression on sites in my top 100 windows for fst. This is not super helpful but whatever. 

require(tidyr)
require(dplyr)
require(data.table)

#First, I need to take the sync and turn it into ref and alt counts. 

dat <- read.delim("../mpileup/wild_controlPop_nomerge.sync", header = F)

#Now I need to split the middle cols. 

colnames(dat) <- c("chr", "pos", "ref", paste("pop", rep(1:(ncol(dat)-3)), sep = ""))

#Now I want to split each of the cols into counts for each 
#Going to use a function to do this then apply it over all the cols.

alleles <- c("A", "T", "C", "G", ".", "NA")

pops <- paste("pop", rep(1:(ncol(dat)-3)), sep = "")


#Now I want to loop over each popualtion and extract the ref and alternate counts 

paste(pops, rep(c("ref", "alt")), sep = "_")
################
crap <- (dat %>%
           filter(ref == "A") )

adf <- crap[,1:3]
      
#loop over for each popualtion and record.    
for (x in 1:length(pops)) {
  pop <- pops[x]
  eachA <- (crap %>%
              dplyr::select(pop)%>%
              separate(pop, into = alleles, convert = TRUE)) 
  
  Aref <- eachA$A
  Aalt <- apply(eachA[,2:6], 1, max)
  
  fuck <- data.frame(Aref, Aalt)
  names(fuck) <- paste(pop, rep(c("ref", "alt"), 1), sep = "_")
  
  adf <- data.frame(adf, fuck)
    
}


############################################
crap <- (dat %>%
           filter(ref == "T") )

tdf <- crap[,1:3]

#loop over for each popualtion and record.    
for (x in 1:length(pops)) {
  pop <- pops[x]
  eachT <- (crap %>%
              dplyr::select(pop)%>%
              separate(pop, into = alleles, convert = TRUE)) 
  
  Tref <- eachT$T
  Talt <- apply(eachT[,c(1,3:6)], 1, max)
  
  fuck <- data.frame(Tref, Talt)
  names(fuck) <- paste(pop, rep(c("ref", "alt"), 1), sep = "_")
  
  tdf <- data.frame(tdf, fuck)
  
}
#######################################

crap <- (dat %>%
           filter(ref == "C") )

cdf <- crap[,1:3]

#loop over for each popualtion and record.    
for (x in 1:length(pops)) {
  pop <- pops[x]
  eachC <- (crap %>%
              dplyr::select(pop)%>%
              separate(pop, into = alleles, convert = TRUE)) 
  
  Cref <- eachC$C
  Calt <- apply(eachC[,c(1:2, 4:6)], 1, max)
  
  fuck <- data.frame(Cref, Calt)
  names(fuck) <- paste(pop, rep(c("ref", "alt"), 1), sep = "_")
  
  cdf <- data.frame(cdf, fuck)
  
}

############################
crap <- (dat %>%
           filter(ref == "G") )

gdf <- crap[,1:3]


for (x in 1:length(pops)) {
  pop <- pops[x]
  eachG <- (crap %>%
              dplyr::select(pop)%>%
              separate(pop, into = alleles, convert = TRUE)) 
  
  Gref <- eachG$G
  Galt <- apply(eachG[,c(1:3, 5:6)], 1, max)
  
  fuck <- data.frame(Gref, Galt)
  names(fuck) <- paste(pop, rep(c("ref", "alt"), 1), sep = "_")
  
  gdf <- data.frame(gdf, fuck)
  
}

#all back together. 
alleleCount <- rbind(adf, tdf, cdf, gdf)

#Now I want to reorder this by chromosome and position 
sortedCounts <- alleleCount[with(alleleCount, order(chr, pos)),]

#need to pull out the core genome. 

coreGenome <- filter(sortedCounts, chr == "X" | chr == "2L"| chr == "2R"| chr == "3L"| chr == "3R"| chr == "4")

#only alt counts for SNP file for SMARTPCA
altCounts <- coreGenome[,c(5,7,9,11)]

#I also want only variant sites because only variants contribute to variation 
#removing rows with only 4 variants because I don't trust those
snpCounts <- altCounts[rowSums(altCounts) > 4, ]

#first I want only variant sites. Going to write this specifically for 4 popualtions (wild data) using the greater than 4 alternate counts to include the data. 

#2567303 sites kept for wild data. 
pass.filter <- coreGenome[rowSums(coreGenome[,c(5,7,9,11)]) >= 4,]

#I need to get allele frequency for each population. 
af.dat <- data.frame(pass.filter[,1:3])

af.dat$CMO <- pass.filter$pop1_alt/pass.filter$pop1_ref

af.dat$FVVW13 <- pass.filter$pop2_alt/pass.filter$pop2_ref

af.dat$FVVW14 <- pass.filter$pop3_alt/pass.filter$pop3_ref

af.dat$PHO <- pass.filter$pop4_alt/pass.filter$pop4_ref


write.csv(af.dat, "wild_controlPop_altAlleleFreq.csv", quote = FALSE, row.names = FALSE)

write.csv(coreGenome, "wild_controlPop_cleanedCounts_coreGenome.csv", quote = FALSE, row.names = FALSE)


data("nancycats")
cats <- nancycats
