#a second shot at doing the PCoA for the wild populations for the paper. 
#Using ape -> has a pcoa function 
library(ape)
library(data.table)
library(tidyverse)
library(vegan)

#allele counts by popualtion 
#pop1 = cmo, 2 = fvw13, 3 = fvw14, 4 = pho
dat <- fread("wild_controlPop_cleanedCounts_coreGenome.csv")

colnames(dat)[4:11] <- c("CMO.ref", "CMO.alt", "FVW13.ref", "FVW13.alt", "FVW14.ref", "FVW14.alt", "PHO.ref", "PHO.alt")

#Now I also want to filter for a min read count here. 
#I want a min read depth of 20 for all reads. 
dat$cmo.count <- dat$CMO.ref +  dat$CMO.alt
dat$fvw13.count <- dat$FVW13.ref + dat$FVW13.alt
dat$fvw14.count <- dat$FVW14.ref + dat$FVW14.alt
dat$pho.count <- dat$PHO.ref + dat$PHO.alt

dat2 <- filter(dat, cmo.count >= 20 & fvw13.count >= 20 & fvw14.count >= 20 & pho.count >= 20) 

#alternate alleles with at least 5 alterate counts
alt.counts.pass <- dat2[rowSums(dat2[,c(5,7,9,11)]) >= 10,]


alt.counts.pass$cmo.altF <- alt.counts.pass$CMO.alt/alt.counts.pass$cmo.count

alt.counts.pass$pho.altF <- alt.counts.pass$PHO.alt/alt.counts.pass$pho.count

alt.counts.pass$fvw13.altF <- alt.counts.pass$FVW13.alt/alt.counts.pass$fvw13.count

alt.counts.pass$fvw14.altF <- alt.counts.pass$FVW14.alt/alt.counts.pass$fvw14.count

#now I am going to try to make the distance matrix for this. (vegan and ape are already on the cluster) 

#using allele freq

af.D <- vegdist(t(alt.counts.pass[,16:19]), "bray")

af.D.dist <- dist(t(alt.counts.pass[,16:19]))


af.pcoa.dist <- pcoa(af.D.dist)
biplot(af.pcoa.dist)

af.pcoa <- pcoa(af.D)
biplot(af.pcoa)
vecs <- data.frame(af.pcoa$vectors)
vecs$population <- rownames(vecs)



pcoa1.2plot <- ggplot(vecs, aes(x = Axis.1, y = Axis.2, label = population)) + 
  geom_point() + 
  geom_text(vjust="inward",hjust="inward")

pcoa2.3plot <- ggplot(vecs, aes(x = Axis.2, y = Axis.3, label = population)) + 
  geom_point() + 
  geom_text(vjust="inward",hjust="inward")

library(cowplot)

png("../Figures/geneticDistance_pcoa_goodCopy.png", res = 100, units = "px")
plot_grid(pcoa1.2plot, pcoa2.3plot, labels = c("A", "B"))
dev.off()

#saved on cluster and moved locally for better plotting options. 

load("../Data/genetic_distances.Rda")


# png("../Figures/genetic_dist_pcoa_distFunction.png")
# biplot(af.pcoa.dist)
# dev.off()


#I want to make a figure plotting the genetic and phenotypic distances. 
#first going to take the distance matrix and turn it into a df with pairwise distance. This is not the most elegant way to make it but it works. 

geno_mat <- as.matrix(af.D)
af.D

#first col is the phenotypic distances beteween populations. 
phenotype_dist <- read.csv("../Tables/WildPopulation_pairwiseTest.csv")

colnames(phenotype_dist)[1] <- "comparisons"
phenotype_dist$comparisons <- gsub("fvw12", "fvw13", phenotype_dist$comparisons)
colnames(phenotype_dist)[2] <- "phenotypic_distance"


dist.df <- data.frame(phenotype_dist[,1])
colnames(dist.df)[1] <- "comparisons"

dist.df$genetic_distance <- c(geno_mat[1,3], geno_mat[1,4], geno_mat[1,2], geno_mat[3,4], geno_mat[2,3], geno_mat[2,4])

dist.df$phenotype_dist <- phenotype_dist[,2]

png("../Figures/genetic_phenotypic_distance.png", res = 100, units = "px")
ggplot(dist.df, aes(x = genetic_distance, y = phenotype_dist, label = comparisons)) + 
  geom_point() + 
  geom_text(vjust = "inward", hjust="inward") + 
  xlab("Genetic Distance (Bray's Distance)") + 
  ylab("Phenotypic Distance (Procrustes Distance)")

dev.off()



