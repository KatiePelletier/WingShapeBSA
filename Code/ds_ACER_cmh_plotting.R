library(data.table)
library(tidyverse)

source('KP_genomescan_source.R')

ds.cmh <- fread("../Data/ds_wild_ACER_zeroGen_fdr.csv")

#reordering the genome so the X is first and adding numbers for plotting. 

ds.cmh.fix <- chrNumbering(ds.cmh)

ds.cmh.middle <- middleChr(ds.cmh)

#need -logp for plotting.

ds.cmh.fix$logp <- -log(ds.cmh.fix$pval)


#filtering the sig sites for plotting. 
hits <- filter(ds.cmh.fix, padj <= 0.05)
hits
#15 sites. 
nrow(hits)


cmhPlot <- ggplot(data = ds.cmh.fix, aes(x=number, y=logp, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("-log(p-val)") +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  geom_point(data = hits, aes(x=number, y=logp), col = "red")+
  scale_x_discrete(limits=c(ds.cmh.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))

png("../Figures/ds_wild_cmh_ACER_zeroGen_withFDR.png", width=1060,height=412,units="px")
cmhPlot
dev.off()

#Now subsetting out ds alone. 

ds <- (ds.cmh.fix %>% 
         filter(chr == "2L" ) %>%
         filter(pos >= 640013) %>%
         filter(pos <= 714983)
)


dsplot <- ggplot(data = ds, aes(x=pos, y=logp)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("postition") +
  ylab("-log(p-val)") +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))

png("../Figures/wild_ds_ACRER_zeroGen_cmh.png", width=1060,height=412,units="px")
dsplot
dev.off()

################################################
ds.cmh.noPHO <- fread("../Data/ds_wild_ACER_zeroGen_fdr_noPHO.csv")

#reordering the genome so the X is first and adding numbers for plotting. 

ds.cmh.noPHO.fix <- chrNumbering(ds.cmh.noPHO)

ds.cmh.middle <- middleChr(ds.cmh.noPHO)

#need -logp for plotting.

ds.cmh.noPHO.fix$logp <- -log(ds.cmh.noPHO.fix$pval.no.pho)


#filtering the sig sites for plotting. 
hits.noPHO <- filter(ds.cmh.noPHO.fix, padj.no.pho <= 0.05)
hits.noPHO
#158 sites. 
nrow(hits.noPHO)


cmhPlot <- ggplot(data = ds.cmh.noPHO.fix, aes(x=number, y=logp, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("-log(p-val)") +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  geom_point(data = hits.noPHO, aes(x=number, y=logp), col = "red")+
  scale_x_discrete(limits=c(ds.cmh.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))

png("../Figures/ds_wild_cmh_ACER_zeroGen_withFDR_NoPHO.png", width=1060,height=412,units="px")
cmhPlot
dev.off()

#Now subsetting out ds alone. 

ds.noPHO <- (ds.cmh.noPHO.fix %>% 
         filter(chr == "2L" ) %>%
         filter(pos >= 640013) %>%
         filter(pos <= 714983)
)


dsplot <- ggplot(data = ds.noPHO, aes(x=pos, y=logp)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("postition") +
  ylab("-log(p-val)") +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))

png("../Figures/wild_ds_ACRER_zeroGen_cmh_noPHO.png", width=1060,height=412,units="px")
dsplot
dev.off()



#####Quick GO analysis using TopGO 

library(topGO)
library(org.Dm.eg.db)

#Reading in the annotated sig sites table 
#this adds an extra col with nothing but the data looks just fine. 
ds.noPHO.sig <- fread("../Data/ds_wild_significantSites_zeroGen_genes_fdr_noPho.table")

ds.sig.genes <- ds.noPHO.sig$FBID
#vector of genes. Hopefully the duplicates don't mess this up. I dont think they will. 
ds.sig.genes

#loading in the GO index we made 
gene_GO <- readMappings("fly_to_GO.delim")


#I need to make a list of all possible genes (all genes in fly) and then classify if these are in my data set or not. 
allgenes <- data.frame(GenomicFeatures::genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene))


#this is all the genes that are near diff sites 
allgenes$diff <- as.numeric(allgenes$gene_id %in% ds.sig.genes)


all.coded <- as.factor(allgenes$diff)
names(all.coded) <- allgenes$gene_id
head(all.coded)

diff.sites <- as.factor(rep(1, length(ds.sig.genes)))
names(diff.sites) <- ds.sig.genes
diff.sites

#test for sites 
gene_filter <- function(allScore){
  return(allScore == 1)
}

#testing for enrichment
allgenes <- new("topGOdata",
                ontology = "BP", 
                allGenes = all.coded,
                annotationFun = annFUN.gene2GO, 
                gene2GO = gene_GO
)
#Specific hippo test 
hippogenes <- genesInTerm(allgenes, "GO:0035329")[[1]]
all <- names(all.coded)
diff.genes <- sigGenes(allgenes)

hippotest <- new("classicCount", testStatistic = GOFisherTest, 
                 name = "fisher",
                 allMembers = all, groupMembers = hippogenes,
                 sigMembers = diff.genes)

contTable(hippotest)
#There is an enrichment for hippo I guess? p = 0.02
runTest(hippotest)

termStat(allgenes, "GO:0035329")

resultFisher <- runTest(allgenes, algorithm = "classic", statistic = "fisher")

#Mostly WNT terms
noPHOtermsTop20 <-GenTable(allgenes, classic = resultFisher, ranksOf = "classic", topNodes = 20)

write.csv(noPHOtermsTop20, file = "../Tables/GOanalysisTop20_cmhSites_fdr_noPHO_ds_wild.csv")

noPHOtermsTop50 <-GenTable(allgenes, classic = resultFisher, ranksOf = "classic", topNodes = 50)
write.csv(noPHOtermsTop50, file = "../Tables/GOanalysisTop50_cmhSites_fdr_noPHO_ds_wild.csv")
