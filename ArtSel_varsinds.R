#Looking at the 10 ds variants (and later all variants) in the artifical selection lineages. These are calls with varscan then I turned into a table and annotated with SNPeff
library(tidyverse)

dsArt <- read_tsv("../Data/dsArtSel_subsetfordspositions.nohead.ann.table")

#changing the sample names to something meaningful. There is probably a way to give this to varscan but I need to look into that and work into the script 
#there is a faster way to do this, but I am not a good enough programmer. 
colnames(dsArt) <- sub("Sample1", "CNA", colnames(dsArt))
colnames(dsArt) <- sub("Sample2", "CNB", colnames(dsArt))
colnames(dsArt) <- sub("Sample3", "CNC", colnames(dsArt))
colnames(dsArt) <- sub("Sample4", "DNA", colnames(dsArt))
colnames(dsArt) <- sub("Sample5", "DNB", colnames(dsArt))
colnames(dsArt) <- sub("Sample6", "DNC", colnames(dsArt))
colnames(dsArt) <- sub("Sample7", "UPA", colnames(dsArt))
colnames(dsArt) <- sub("Sample8", "UPB", colnames(dsArt))
colnames(dsArt) <- sub("Sample9", "UPC", colnames(dsArt))

#for some reason pos is read as 0?
colnames(dsArt)[2] <- "pos"

#Looks good. 
(names(dsArt))

emc <- read_tsv("../Data/emc_subsetfordspositions.nohead.ann.table")
colnames(emc) <- sub("Sample1", "DNA", colnames(emc))
colnames(emc) <- sub("Sample2", "DNB", colnames(emc))
colnames(emc) <- sub("Sample3", "DNC", colnames(emc))
colnames(emc) <- sub("Sample4", "UPA", colnames(emc))
colnames(emc) <- sub("Sample5", "UPB", colnames(emc))
colnames(emc) <- sub("Sample6", "UPC", colnames(emc))


colnames(emc)[2] <- "pos"
#looks good. 
names(emc)

#These are the SNPs in Pitchers 2019 (subset for the ds SNPs)
snps <- read.csv("../Data/SNPsFromPitchers.csv", header = TRUE)

#Subsetting the variant calls to get the 10 SNPs 

#only the INDEL? 
snps.ds <- filter(dsArt, dsArt$pos %in% snps$FB6)

#Everything else must be masked? 
snps.emc <- filter(emc, emc$pos %in% snps$FB6)

#So the indel looks really promising in the ds data but the alt allele is fixed (or nearly fixed) in all of the emc data 

#for now, I want to pull the ds freq to see if I can estimate magnitude and sign. 

ds.pops <- c("CNA", "CNB", "CNC", "DNA", "DNB", "DNC", "UPA", "UPB", "UPC")
ds.freq <- snps.ds[c(18,26,34,42,50,58,66,74,82)]

#stripping off the % and making this a number. Why this isn't the defult beats me. 
ds.freq <- as.numeric(sapply(ds.freq, function(x) {x <- gsub("%", "", x)}))/100

#now I can add this all together. 

ds.indel.freq <- data.frame(ds.pops, ds.freq)
str(ds.indel.freq)
ds.indel.freq$selection <- c(rep("CN", 3), rep("DN", 3), rep("UP", 3))    
ds.indel.freq$line <- rep(c("A", "B", "C"), 3)
ds.indel.freq


#This should really have a random line effect. This is imperfect because the As are unrelated              
#But the sign and the effect is there 

#indel.mod <- lmer(ds.freq ~ selection + (selection|line), dat = ds.indel.freq)
ds.indel.freq$selection <- as.factor(ds.indel.freq$selection)
ds.indel.freq$line <- as.factor(ds.indel.freq$line)


indel.mod <- lm(ds.freq ~ selection + (selection|line), dat = ds.indel.freq)


summary(indel.mod)
Anova(indel.mod)

#Imperfect but I think this is something? Talk to Ian? 
plot(emmeans(indel.mod, ~selection))
                    
                      