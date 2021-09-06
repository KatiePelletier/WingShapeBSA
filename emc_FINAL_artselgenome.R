#This is me trying to finish this part of the project finally. 
library(data.table)
library(tidyverse)

#my plotting functions
source("KP_genomescan_source.R")

#Starting with what might be the easiest one, the ee populations. 
#This is in 1000bp windows with a 1000bp step. 
#Pretty much I want to just draw this and identify mean + 2sd to use as my GO cutoff. 
#This file has all of the comparisons included to make things easier. 

# art.sel <- fread("../Data/emc_merge_1000window.fst")
# 
# 
# art.sel <- populaionFst_cleanup(art.sel, x = c('cndn', 'cnup', 'dnup'))
# 
# #Reordering and numbering the chr for plotting. 
# fixed_art.sel <- chrNumbering(art.sel)
# 
# chrlabel1000 <- middleChr(fixed_art.sel)
# 
# #Make sure nothing crazy happened. 
# #wow this actually looks so much better. look at those peaks! 
# genomePlot <- ggplot(data = fixed_art.sel, aes(x=number, y=dnup, color=chr)) + 
#   geom_point(size=1, show.legend = F, alpha = 0.6) + 
#   theme(panel.background = element_blank()) +
#   xlab("Chromosome") +
#   ylab("meanFst") +
#   scale_x_discrete(limits=c(chrlabel1000),
#                    labels = c("X","2L", "2R", '3L', '3R', '4')) +
#   scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
#   theme(text = element_text(size=15),
#         axis.text.x= element_text(size=12), 
#         axis.text.y= element_text(size=12), 
#         panel.border = element_rect(colour = "black",
#                                     fill=NA, size=0.5))
# 
# 
# #hippo non polarity
# 
# #2L:640,013..714,983
# ds <- subset(fixed_art.sel, chr == "2L" & window == 640500)
# #2L:4,198,404..4,221,796
# ft <- subset(fixed_art.sel, chr == '2L' & window == 4198500)
# #2L:16,466,512..16,485,998
# dac <- subset(fixed_art.sel, chr == '2L' & window == 1646500)
# #2R:18,232,761..18,236,311
# fj <- subset(fixed_art.sel, chr == "2R" & window == 18233500)
# #3R:31,054,085..31,061,210
# dco <- subset(fixed_art.sel, chr == "3R" & window == 31054500)
# #2L:10,364,471..10,366,410 
# lft <- subset(fixed_art.sel, chr == "2L" & window == 10364500)
# #	3L:12,206,238..12,266,841
# app <- subset(fixed_art.sel, chr == "3L" & window == 12206500)
# #2L:9,839..21,376 
# #this is l(2)gl
# lgl <- subset(fixed_art.sel, chr == "2L" & window == 10500)
# #3R:26,536,350..26,604,051
# scrib <- subset(fixed_art.sel, chr == "3R" & window == 26536500)
# #	2R:24,065,975..24,068,485
# yki <- subset(fixed_art.sel, chr == "2R" & window == 24066500)
# #X:15,804,370..15,827,682
# sd <- subset(fixed_art.sel, chr == "X" & window == 15804500)
# #3R:10,507,561..10,639,568 
# hth <- subset(fixed_art.sel, chr == "3R" & window == 10508500)
# #2L:21,828,593..21,837,011
# tsh <- subset(fixed_art.sel, chr == "2L" & window == 21839500)
# #3R:30,789,625..30,806,619
# warts <- subset(fixed_art.sel, chr == "3R" & window == 30790500)
# #3R:22,421,161..22,423,302
# mats <- subset(fixed_art.sel, chr == "3R" & window == 22421500)
# #2R:19,493,996..19,496,856
# hippo <- subset(fixed_art.sel, chr == "2R" & window == 19494500)
# #3R:23,058,609..23,061,304 
# sav <- subset(fixed_art.sel, chr == "3R" & window == 23059500)
# #3R:14,697,689..14,724,261
# kibra <- subset(fixed_art.sel, chr == "3R" & window == 14698500)
# #2L:431,227..448,701
# ex <- subset(fixed_art.sel, chr == "2L" & window == 431500)
# #X:19,689,697..19,693,500
# mer <- subset(fixed_art.sel, chr == "X" & window == 19690500)
# #3R:24,295,078..24,314,541
# crb <- subset(fixed_art.sel, chr == "3R" & window == 24295500)
# #3R:23,114,014..23,117,931
# rassf <- subset(fixed_art.sel, chr == "3R" & window == 23114500)
# #X:13,826,041..13,830,317
# jub <- subset(fixed_art.sel, chr == "X" & window == 13826500)
# 
# dsgenes.number <- rep(NA, 17) 
# dsgenes.number <- c(ds$number, ft$number, dac$number, fj$number, dco$number, lft$number, app$number, lgl$number, scrib$number, yki$number, sd$number, hth$number, tsh$number, warts$number, mats$number, hippo$number, sav$number, kibra$number, ex$number, mer$number, crb$number, rassf$number, jub$number)
# 
# #hippo genes not in fig 
# 
# #2L:4,477,462..4,614,300 
# dpy <- subset(fixed_art.sel, chr == '2L' & window == 4477500)
# #4:1,057,365..1,065,001
# zyx <- subset(fixed_art.sel, chr == '4' & window == 1057500)
# 
# #not hippo genes 
# 
# #2R:21,522,420..21,559,977
# egfr <- subset(fixed_art.sel, chr == '2R' & window == 21522500)
# #3L:749,400..753,492
# emc <- subset(fixed_art.sel, chr == '3L' & window == 749500)
# #3R:9,020,348..9,039,471
# neur <- subset(fixed_art.sel, chr == '3R' & window == 9020348)
# 
# nothippo.number <- rep(NA, 3)
# nothippo.number <- c(egfr$number, emc$number, neur$number)
# 
# m2 <- mean(fixed_art.sel$dnup)
# s2 <- sd(fixed_art.sel$dnup)
# cutoff2 <- m2 + 2*s2
# 
# 
# fixed_art.sel$outlier <- ifelse(fixed_art.sel$dnup > cutoff2, "yes", "no")
# 
# 
# #I sort of think I want to keep 5000 bp windows. 
# linescan <- ggplot(data = fixed_art.sel, aes(x=number, y=dnup, color=chr, alpha = outlier)) + 
#   geom_point(size=1, show.legend = F) + 
#   scale_alpha_discrete(range = c(0.1,0.5)) +
#   theme(panel.background = element_blank()) +
#   #scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1)) +
#   xlab("Chromosome") +
#   ylab(expression(F[ST])) +
#   scale_x_discrete(limits=c(chrlabel1000),
#                    labels = c("X","2L", "2R", '3L', '3R', '4')) +
#   scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
#   #geom_segment(aes(x = egfrline$number, y = 0, xend = egfrline$number, yend = egfrline$meanFst, color = 'red', alpha = 0.5)) +
#   #geom_vline(xintercept = dsgenes.number, size = 1, color = 'red', alpha = 0.5) +
#   #geom_vline(xintercept = egfr$number, size = 1, color = 'blue', alpha = 0.5) +
#   geom_vline(xintercept = emc$number, size = 1, color = 'purple', alpha = 0.5) +
#   geom_vline(xintercept = ds$number, size = 1, color = 'red', alpha = 0.5) +
#   geom_hline(yintercept = cutoff2, alpha = 0.5) +
#   theme(text = element_text(size=15),
#         axis.text.x= element_text(size=12), 
#         axis.text.y= element_text(size=12), 
#         panel.border = element_rect(colour = "black", 
#                                     fill=NA, size=0.5))
# 
# png(file = "../Output/emc_1000windows_dsemc.png",width=1060,height=412,units="px")
# linescan
# dev.off()



#I want to look at this with 5000 bp windows. 
artsel.5000 <- fread("../Data/emc_merge_5000window.fst")

artsel.5000 <-  populaionFst_cleanup(artsel.5000, x = c('cndn', 'cnup', 'dnup'))


fixed_artsel.5000 <- chrNumbering(artsel.5000)

chrlabel <- middleChr(fixed_artsel.5000)

#This looks really nice
genomePlot <- ggplot(data = fixed_artsel.5000, aes(x=number, y=dnup, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("meanFst") +
  scale_x_discrete(limits=c(chrlabel),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))


#2L:640,013..714,983
ds <- subset(fixed_artsel.5000, chr == "2L" & window == 642500)
#2L:4,198,404..4,221,796
ft <- subset(fixed_artsel.5000, chr == '2L' & window == 4202500)
#2L:16,466,512..16,485,998
dac <- subset(fixed_artsel.5000, chr == '2L' & window == 1647500)
#2R:18,232,761..18,236,311
fj <- subset(fixed_artsel.5000, chr == "2R" & window == 18237500)
#3R:31,054,085..31,061,210
dco <- subset(fixed_artsel.5000, chr == "3R" & window == 31057500)
#2L:10,364,471..10,366,410 
lft <- subset(fixed_artsel.5000, chr == "2L" & window == 10367500)
#	3L:12,206,238..12,266,841
app <- subset(fixed_artsel.5000, chr == "3L" & window == 12207500)
#2L:9,839..21,376 
#this is l(2)gl
lgl <- subset(fixed_artsel.5000, chr == "2L" & window == 12500)
#3R:26,536,350..26,604,051
scrib <- subset(fixed_artsel.5000, chr == "3R" & window == 26537500)
#	2R:24,065,975..24,068,485
yki <- subset(fixed_artsel.5000, chr == "2R" & window == 24067500)
#X:15,804,370..15,827,682
sd <- subset(fixed_artsel.5000, chr == "X" & window == 15807500)
#3R:10,507,561..10,639,568 
hth <- subset(fixed_artsel.5000, chr == "3R" & window == 10512500)
#2L:21,828,593..21,837,011
tsh <- subset(fixed_artsel.5000, chr == "2L" & window == 21842500)
#3R:30,789,625..30,806,619
warts <- subset(fixed_artsel.5000, chr == "3R" & window == 30792500)
#3R:22,421,161..22,423,302
mats <- subset(fixed_artsel.5000, chr == "3R" & window == 22422500)
#2R:19,493,996..19,496,856
hippo <- subset(fixed_artsel.5000, chr == "2R" & window == 19497500)
#3R:23,058,609..23,061,304 
sav <- subset(fixed_artsel.5000, chr == "3R" & window == 23062500)
#3R:14,697,689..14,724,261
kibra <- subset(fixed_artsel.5000, chr == "3R" & window == 14702500)
#2L:431,227..448,701
ex <- subset(fixed_artsel.5000, chr == "2L" & window == 432500)
#X:19,689,697..19,693,500
mer <- subset(fixed_artsel.5000, chr == "X" & window == 19692500)
#3R:24,295,078..24,314,541
crb <- subset(fixed_artsel.5000, chr == "3R" & window == 24297500)
#3R:23,114,014..23,117,931
rassf <- subset(fixed_artsel.5000, chr == "3R" & window == 23117500)
#X:13,826,041..13,830,317
jub <- subset(fixed_artsel.5000, chr == "X" & window == 13827500)

dsgenes.number <- rep(NA, 17) 
dsgenes.number <- c(ds$number, ft$number, dac$number, fj$number, dco$number, lft$number, app$number, lgl$number, scrib$number, yki$number, sd$number, hth$number, tsh$number, warts$number, mats$number, hippo$number, sav$number, kibra$number, ex$number, mer$number, crb$number, rassf$number, jub$number)

#hippo genes not in fig 

#2L:4,477,462..4,614,300 
dpy <- subset(fixed_artsel.5000, chr == '2L' & window == 4477500)
#4:1,057,365..1,065,001
zyx <- subset(fixed_artsel.5000, chr == '4' & window == 1057500)

#not hippo genes 

#2R:21,522,420..21,559,977
egfr <- subset(fixed_artsel.5000, chr == '2R' & window == 21522500)
#3L:749,400..753,492
emc <- subset(fixed_artsel.5000, chr == '3L' & window == 752500)
#3R:9,020,348..9,039,471
neur <- subset(fixed_artsel.5000, chr == '3R' & window == 9025000)

nothippo.number <- rep(NA, 3)
nothippo.number <- c(egfr$number, emc$number, neur$number)


#I also want to find where the 2sd cutoff is. this is based on a recomendation in Cutter pop gen textbook.

m <- mean(fixed_artsel.5000$dnup)
s <- sd(fixed_artsel.5000$dnup)

cutoff <- m + 3*s
cutoff
cutoff2 <-  + 2*s

#similar in 1000 bp data? 0.22 vs 0.24 so close-ish? Shouldn't be identical 

fixed_artsel.5000$outlier <- ifelse(fixed_artsel.5000$dnup > cutoff, "yes", "no")

#I sort of think I want to keep 5000 bp windows. 
linescan <- ggplot(data = fixed_artsel.5000, aes(x=number, y=dnup, color=chr, alpha = outlier)) + 
  geom_point(size=1, show.legend = F) + 
  scale_alpha_discrete(range = c(0.2,0.7)) +
  theme(panel.background = element_blank()) +
  #scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1)) +
  xlab("Chromosome") +
  ylab(expression(F[ST])) +
  scale_x_discrete(limits=c(chrlabel),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  #geom_segment(aes(x = egfrline$number, y = 0, xend = egfrline$number, yend = egfrline$meanFst, color = 'red', alpha = 0.5)) +
  #geom_vline(xintercept = dsgenes.number, size = 1, color = 'red', alpha = 0.5) +
  #geom_vline(xintercept = egfr$number, size = 1, color = 'blue', alpha = 0.5) +
  geom_vline(xintercept = emc$number, size = 1, color = 'purple', alpha = 0.5) +
  geom_vline(xintercept = ds$number, size = 1, color = 'red', alpha = 0.5) + 
  geom_hline(yintercept = cutoff2, alpha = 0.5) +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12), 
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=0.5))

png(file = "../Output/emc_5000windows_dsemc_3sd.png",width=1060,height=412,units="px")
linescan
dev.off()

#############GO Analysis######################
library(bumphunter)

#The first thing I want to do is take what I want to call a 'peak' 
# I need a good cutoff for what this Fst value should be... for now I'm just going to pick one. 0.3 seems conservitive but thats ok with me for now. 



peaks <- filter(fixed_artsel.5000, dnup >= cutoff2)
#Now I need to actaully make this a range of positions covered. Popoolation gives the middle point for the window as the position. So the range covered is +/- 2500 bp from there. 

peaks$Lwindow <- peaks$window - 2500
peaks$Rwindow <- peaks$window + 2500



#now I want to combine these together into chunks of the genome, this is done with bumphunter. 
#This makes an indexing table for each line in the peaks file, to organize them into clusters

test <- fixed_artsel.5000

#This max gap term may need to be adjusted, still lots of singletons. 
c1 <- clusterMaker(test$chr, test$window, maxGap = 10000)
table(c1)

#now that its indexed, I need to find the segments in bumps... 

segs <- getSegments(test$dnup, c1, cutoff=cutoff2)


#Now I can make a table of these regions. 
tab <- regionFinder(test$dnup, test$chr, test$window, c1, cutoff=cutoff)
head(tab)

#a annotation database for the fly genome 
#BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#BiocManager::install("org.Dm.eg.db")
library(org.Dm.eg.db)


genes <- annotateTranscripts(TxDb.Dmelanogaster.UCSC.dm6.ensGene, 
                             #annotationPackage = org.Dm.eg.db,
                             by = "gene", codingOnly = FALSE)

#This didn't match the Grange object so adding chr at the start 
tab$chr <- paste0("chr", tab$chr)

#now going to annotate the genes in my data set. 

#Need to trun table into Grange object 
#first creating a list of all the genomic intrivals. 
kp.tab <- rep(NA, nrow(tab))
kp.tab <- paste(tab$chr, tab$start, tab$end, sep = ":")

#now coverting this into a Grange object. 

myranges <- makeGRangesFromDataFrame(tab, keep.extra.columns = FALSE, 
                                     ignore.strand = TRUE, 
                                     start.field = "start", 
                                     end.field = "end")

#finding the overlap between genes and my ranges. 
#here are the genes. 
#genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

myresults <- subsetByOverlaps(genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene), myranges)

#the metadata col contains all the genes. 
myresults
mygenes <- myresults$gene_id

#write.csv(mygenes, "../Tables/emc_3sd_outliergenelist.csv")

#going to try to do the GO analysis with TopGO using my 722 tutorial 
library(topGO)

#loading in the GO index we made 
gene_GO <- readMappings("fly_to_GO.delim")


#I need to make a list of all possible genes (all genes in fly) and then classify if these are in my data set or not. 
allgenes <- data.frame(GenomicFeatures::genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene))

allgenes$peak <- as.numeric(allgenes$gene_id %in% mygenes)

#topGO wants the gene names as rownames
test <- data.frame(allgenes$peak)
rownames(test) <- allgenes$gene_id

test <- as.factor(allgenes$peak)
names(test) <- allgenes$gene_id

interesting <- as.factor(rep(1, length(mygenes)))
names(interesting) <- mygenes

gene_filter <- function(allScore){
  return(allScore == 1)
}

allgenes <- new("topGOdata",
                ontology = "BP", 
                allGenes = test,
                annotationFun = annFUN.gene2GO, 
                gene2GO = gene_GO
)

#All the fly genes
head(genes(allgenes))
numGenes(allgenes)
numSigGenes(allgenes)

#This drops a lot more genes than the ds one.  
length(inpeaks)
length(mygenes)

length(base::intersect(inpeaks, mygenes))

#going to see what happens?
#GO term for hippo 
#Given my other results, I would expect this.  
hippogenes <- genesInTerm(allgenes, "GO:0035329")[[1]]
all <- genes(allgenes)
inpeaks <- sigGenes(allgenes)

#This doesn't make sense with what I saw before. 
hippotest <- new("classicCount", testStatistic = GOFisherTest, 
                 name = "fisher",
                 allMembers = all, groupMembers = hippogenes,
                 sigMembers = inpeaks)

contTable(hippotest)
runTest(hippotest)

#I guess this counts? 
termStat(allgenes, "GO:0035329")


#what are all the enriched terms?

resultFisher <- runTest(allgenes, algorithm = "classic", statistic = "fisher")

allRes20_2sd <- GenTable(allgenes, classic = resultFisher, ranksOf = "classic", topNodes = 20)

write.csv(allRes20_2sd, file = "../Output/emc_artselGOtop20_2sd.csv")

allRes50_2sd <- GenTable(allgenes, classic = resultFisher, ranksOf = "classic", topNodes = 50)

write.csv(allRes50_2sd, file = "../Output/emc_artselGOtop50_2sd.csv")


allRes100_2sd <- GenTable(allgenes, classic = resultFisher, ranksOf = "classic", topNodes = 100)

write.csv(allRes100_2sd, file = "../Output/emc_artselGOtop100_2sd.csv")



