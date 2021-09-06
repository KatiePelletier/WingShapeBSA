#looking at the SFS in ds in the artifical selection lineages. 


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

#Subsetting for only frequencies and shit I care about 
ds.dat <- dsArt[,c(1:5,8,18,26,34,42,50, 58, 66, 74, 82 )]

ds.dat[,7:15] <- sapply(ds.dat[,7:15], function(x) {x <- as.numeric(gsub("%", "", x))})

#converting to decimal to be in alt. allele frequency. 
ds.dat[,7:15] <- ds.dat[,7:15]/100

#Lots at interemediate freq in the Dn and control linages but lots driven to fix (or near) in up... at least for A as a quick check 
hist(ds.dat$CNA.FREQ)
hist(ds.dat$DNA.FREQ)
hist(ds.dat$UPA.FREQ)

#converting to MAF (if AF >50 taking the inverse) 
ds.maf <- ds.dat


ds.maf[,7:15] <- lapply(ds.maf[,7:15], function(x) {ifelse(x > 0.5, 1 - x, x)})

#This looks better 
hist(ds.maf$CNA.FREQ)
hist(ds.maf$DNA.FREQ)
hist(ds.maf$UPA.FREQ)

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

emc.dat <- emc[,c(1:5,8,18,26,34,42,50, 58)]

emc.dat[,7:12] <- sapply(emc.dat[,7:12], function(x) {x <- as.numeric(gsub("%", "", x))})

#converting to decimal to be in alt. allele frequency. 
emc.dat[,7:12] <- emc.dat[,7:12]/100

#ds control is here because I think it was the same? 
#much less clear in this case. Again makes sense. 
hist(ds.maf$CNA.FREQ)
hist(emc.dat$DNA.FREQ)
hist(emc.dat$UPA.FREQ)

#converting to MAF (if AF >50 taking the inverse) 
emc.maf <- emc.dat


emc.maf[,7:12] <- lapply(emc.maf[,7:12], function(x) {ifelse(x > 0.5, 1 - x, x)})

#Way less clear of a picture. Something is weird with the down data here. 
hist(ds.maf$CNA.FREQ)
hist(emc.maf$DNA.FREQ)
hist(emc.maf$UPA.FREQ)


#Now to draw something a little nicer. 

ggplot(ds.maf, aes(x = CNA.FREQ)) + 
  geom_histogram(binwidth = 0.05, fill = "black") + 
  geom_histogram(aes(x = DNA.FREQ), binwidth = 0.05, 
                 fill = "grey", alpha = 0.9) + 
  geom_histogram(aes(x = UPA.FREQ), binwidth = 0.05,
                 fill = "red", alpha = 0.6) + 
  xlab("MAF") + 
  theme_classic()


#Now with just Up and contol 
ggplot(ds.maf, aes(x = CNA.FREQ)) + 
  geom_histogram(binwidth = 0.05, fill = "black") + 
  geom_histogram(aes(x = UPA.FREQ), binwidth = 0.05,
                 fill = "red", alpha = 0.6) + 
  xlab("MAF")+ 
  theme_classic()


#now averaging over all three replicates (probably a better thing to do here?)

ds.maf$up.avg <- (ds.maf$UPA.FREQ + ds.maf$UPB.FREQ + ds.maf$UPC.FREQ)/3
ds.maf$dn.avg <- (ds.maf$DNA.FREQ + ds.maf$DNB.FREQ + ds.maf$DNC.FREQ)/3
ds.maf$cn.avg <- (ds.maf$CNA.FREQ + ds.maf$CNB.FREQ + ds.maf$CNC.FREQ)/3


#Now to draw something a little nicer. 

ggplot(ds.maf, aes(x = cn.avg)) + 
  geom_histogram(binwidth = 0.05, fill = "black") + 
  geom_histogram(aes(x = dn.avg), binwidth = 0.05, 
                 fill = "grey", alpha = 0.9) + 
  geom_histogram(aes(x = up.avg), binwidth = 0.05,
                 fill = "red", alpha = 0.6) + 
  xlab("MAF") + 
  theme_classic()


#Now with just Up and contol 
png("MAF_atds_ArtSel.png")

ggplot(ds.maf, aes(x = cn.avg)) + 
  geom_histogram(binwidth = 0.05, fill = "black") + 
  geom_histogram(aes(x = up.avg), binwidth = 0.05,
                 fill = "red", alpha = 0.6) + 
  xlab("MAF")+ 
  theme_classic()

dev.off()

#####A nicer emc plot 


emc.maf$up.avg <- (emc.maf$UPA.FREQ + emc.maf$UPB.FREQ + emc.maf$UPC.FREQ)/3
emc.maf$dn.avg <- (emc.maf$DNA.FREQ + emc.maf$DNB.FREQ + emc.maf$DNC.FREQ)/3



#Now to draw something a little nicer. 
#Just down and up here. Have to think more about how to add control. 
ggplot(emc.maf, aes(x = dn.avg)) + 
  geom_histogram(binwidth = 0.05, fill = "black") + 
  geom_histogram(aes(x = up.avg), binwidth = 0.05,
                 fill = "red", alpha = 0.6) + 
  xlab("MAF")+ 
  theme_classic()

