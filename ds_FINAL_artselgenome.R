#This is me trying to finish this part of the project finally. 
library(data.table)
library(tidyverse)

#my plotting functions
source("KP_genomescan_source.R")


#a annotation database for the fly genome 
#BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#BiocManager::install("org.Dm.eg.db")
library(org.Dm.eg.db)
# library(bumphunter)
# 
# #if below code is run, this fails. 
# genes <- annotateTranscripts(TxDb.Dmelanogaster.UCSC.dm6.ensGene, 
#                              #annotationPackage = org.Dm.eg.db,
#                              by = "gene", codingOnly = FALSE)

library(topGO)
gene_GO <- readMappings("fly_to_GO.delim")

#Starting with what might be the easiest one, the ee populations.
#This is in 1000bp windows with a 1000bp step.
#Pretty much I want to just draw this and identify mean + 2sd to use as my GO cutoff.
#This file has all of the comparisons included to make things easier.

# art.sel <- fread("../Data/ee_merged_1000windows.fst")
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
#   scale_x_discrete(limits=c(chrlabel),
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
# fixed_art.sel$outlier <- ifelse(fixed_art.sel$dnup > cutoff2, "yes", "no")


#I sort of think I want to keep 5000 bp windows.
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
#   #geom_vline(xintercept = emc$number, size = 1, color = 'purple', alpha = 0.5) +
#   geom_vline(xintercept = ds$number, size = 1, color = 'red', alpha = 0.5) +
#   geom_hline(yintercept = cutoff2, alpha = 0.5) +
#   theme(text = element_text(size=15),
#         axis.text.x= element_text(size=12),
#         axis.text.y= element_text(size=12),
#         panel.border = element_rect(colour = "black",
#                                     fill=NA, size=0.5))

# png(file = "../Output/1000windows_dsonly",width=1060,height=412,units="px")
# linescan
# dev.off()

#I need to get the mean from all sites 
#this actually gives the same answer. 
all.sites <- fread("../Data/ds_all_wild_ee.mask_5000windows_allsites.fst")

#0s become NAs
all.sites <-  populaionFst_cleanup(all.sites, x = c('cndn', 'cnup', 'dnup'))

to.zero <- function(x){
  gsub(x, NA, 0)
}

sum(is.na((all.sites$dnup)))

all.sites$dnup[is.na(all.sites$dnup)] <- 0

sum(is.na((all.sites$dnup)))

m <- mean(all.sites$dnup)
s <- sd(all.sites$dnup, na.rm = TRUE)

#I want to look at this with 5000 bp windows. 
artsel.5000 <- fread("../Data/ee_merged_5000windows.fst")

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

hipponods.number <- rep(NA, 16)
hipponods.number <- c(ft$number, dac$number, fj$number, dco$number, lft$number, app$number, lgl$number, scrib$number, yki$number, sd$number, hth$number, tsh$number, warts$number, mats$number, hippo$number, sav$number, kibra$number, ex$number, mer$number, crb$number, rassf$number, jub$number)

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

cutoff <- m + 3*s # 0.3449325

cutoff2 <- m + 2*s
#now to define values above cutoff for better plotting 

#similar in 1000 bp data? 0.277 vs 0.263 so close-ish? Shouldn't be identical 
# m2 <- mean(fixed_art.sel$dnup)
# s2 <- sd(fixed_art.sel$dnup)
# cutoff2 <- m2 + 2*s2

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
  #geom_vline(xintercept = hipponods.number, size = 1, color = 'blue', alpha = 0.5) +
  #geom_vline(xintercept = egfr$number, size = 1, color = 'blue', alpha = 0.5) +
  #geom_vline(xintercept = emc$number, size = 1, color = 'purple', alpha = 0.5) +
  geom_vline(xintercept = ds$number, size = 1, color = 'red', alpha = 0.7) +
  geom_hline(yintercept = cutoff, alpha = 0.5) +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12), 
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=0.5))

#png(file = "../Figures/5000windows_dsmarked_3sd_allsites.png",width=1060,height=412,units="px")
linescan
#dev.off()
#
#Now with all the lines 
ds_all_lines <- ggplot(data = fixed_artsel.5000, aes(x=number, y=dnup, color=chr, alpha = outlier)) + 
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
  geom_vline(xintercept = dsgenes.number, size = 1, color = 'red', alpha = 0.5) +
  #geom_vline(xintercept = hipponods.number, size = 1, color = 'blue', alpha = 0.5) +
  #geom_vline(xintercept = egfr$number, size = 1, color = 'blue', alpha = 0.5) +
  geom_vline(xintercept = emc$number, size = 1, color = 'purple', alpha = 0.5) +
  geom_vline(xintercept = ds$number, size = 1, color = 'red', alpha = 0.7) +
  geom_hline(yintercept = cutoff, alpha = 0.5) +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12), 
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=0.5))

#png(file = "../Figures/5000windows_allmarked_3sd_allsites.png",width=1060,height=412,units="px")
ds_all_lines
#dev.off()


########################top 100000 sites########################
#for the logistic regression, we want to take the top 10000 windows
#This is way too many for the number of observations. remember that each site is *5000 for the regression

# sitesforreg <- top_frac(fixed_artsel.5000, 0.01, dnup)
# #top 1% of sites is still over 10000000 so I dont want that. 
# nrow(sitesforreg)
# 241*5000

#going to start with the top 100 windows and we can always scale up later. 
sitesforreg <- top_n(fixed_artsel.5000, 100, dnup)

hist(sitesforreg$dnup)
hist(fixed_artsel.5000$dnup)

quantile(sitesforreg$dnup)
quantile(fixed_artsel.5000$dnup)

#Now I need to write this out in a form that will be useful. 
#I want this to be a gtf/gff2 because that works with population for filtering. 
#these are tab seperated with these cols:
# seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. 
# source - name of the program that generated this feature, or the data source (database or project name)
# feature - feature type name, e.g. Gene, Variation, Similarity
# start - Start position* of the feature, with sequence numbering starting at 1.
# end - End position* of the feature, with sequence numbering starting at 1.
# score - A floating point value.
# strand - defined as + (forward) or - (reverse).
# frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
# attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

#Now to build this into a data frame. 
#I'll start with an empty matrix and fill it 

gff <- data.frame(sitesforreg$chr, #chromosome
           rep("popoolationFst", 100), #program
           rep("highDiff", 100), #feature
           sitesforreg$window - 2500, #start
           sitesforreg$window + 2500, #end
           sitesforreg$dnup, #score
           rep("+", 100),#strand, doesnt matter 
           rep(0, 100), #frame, doesnt matter
           rep("dnVup", 100)#attribute, place holder
           )

str(gff)
dim(gff)


#write_tsv(gff, "../Data/top100windows_ds_5000windowFstCalc.gff", 
 #         col_names = FALSE)

#####################################################################

head(fixed_artsel.5000)
fixed_artsel.5000$Lwindow <- fixed_artsel.5000$window - 2500
fixed_artsel.5000$Rwindow <- fixed_artsel.5000$window + 2500

sig.sites2sd <- fixed_artsel.5000[fixed_artsel.5000$dnup >= cutoff2,]
sig.sites2sd$chr <- paste0("chr", sig.sites2sd$chr)

twosdRanges <- makeGRangesFromDataFrame(sig.sites2sd,
                                     keep.extra.columns = FALSE,
                                     ignore.strand = TRUE,
                                     start.field = "Lwindow",
                                     end.field = "Rwindow")

twosdResults <- subsetByOverlaps(GenomicFeatures::genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene), twosdRanges)

twosdgenes <- twosdResults$gene_id

allgenes <- data.frame(GenomicFeatures::genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene))

allgenes$twosd <- as.numeric(allgenes$gene_id %in% twosdgenes)

twosd <- data.frame(allgenes$twosd)
rownames(twosd) <- allgenes$gene_id
twosd <- as.factor(allgenes$twosd)
names(twosd) <- allgenes$gene_id

genes.2sd <- as.factor(rep(1, length(twosdgenes)))
names(genes.2sd) <- twosdgenes

gene_filter <- function(allScore){
  return(allScore == 1)
}

go2sd <- new("topGOdata",
                ontology = "BP", 
                allGenes = twosd,
                annotationFun = annFUN.gene2GO, 
                geneSelectionFun = gene_filter1,
                gene2GO = gene_GO
)

hippogenes <- genesInTerm(go2sd, "GO:0035329")[[1]]
all <- genes(go2sd)
topGenes <- sigGenes(go2sd)

length(base::intersect(topGenes, twosdgenes))

hippotest.2sd <- new("classicCount", testStatistic = GOFisherTest,
                 name = "fisher",
                 allMembers = all, groupMembers = hippogenes,
                 sigMembers = twosdgenes)

contTable(hippotest.2sd)
runTest(hippotest.2sd)

termStat(go2sd, "GO:0035329")

resultFisher2sd <- runTest(go2sd, algorithm = "classic", statistic = "fisher")

allRes20_2sd <- GenTable(go2sd, classic = resultFisher2sd, ranksOf = "classic", topNodes = 20)

write.csv(allRes20_2sd, file = "../Tables/ds_artselGOtop20_2sd.csv")

allRes50_2sd <- GenTable(go2sd, classic = resultFisher2sd, ranksOf = "classic", topNodes = 50)

write.csv(allRes50_2sd, file = "../Tables/ds_artselGOtop50_2sd.csv")

allRes100_2sd <- GenTable(go2sd, classic = resultFisher2sd, ranksOf = "classic", topNodes = 100)

write.csv(allRes100_2sd, file = "../Tables/ds_artselGOtop100_2sd.csv")

obs.SigHippo2sd <- as.numeric(termStat(go2sd, "GO:0035329")[2])
obs.HippoRatio2sd <- as.numeric(termStat(go2sd, "GO:0035329")[2])/as.numeric(termStat(go2sd, "GO:0035329")[3])

obs.SigNegHippo2sd <- as.numeric(termStat(go2sd, "GO:0035331")[2])
obs.NegHippoRatio2sd <- as.numeric(termStat(go2sd, "GO:0035331")[2])/as.numeric(termStat(go2sd, "GO:0035331")[3])

#############3sd####################

sig.sites3sd <- fixed_artsel.5000[fixed_artsel.5000$dnup >= cutoff,]
sig.sites3sd$chr <- paste0("chr", sig.sites3sd$chr)

threesdRanges <- makeGRangesFromDataFrame(sig.sites3sd,
                                          keep.extra.columns = FALSE,
                                          ignore.strand = TRUE,
                                          start.field = "Lwindow",
                                          end.field = "Rwindow")

threesdResults <- subsetByOverlaps(GenomicFeatures::genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene), threesdRanges)

threesdgenes <- threesdResults$gene_id

allgenes <- data.frame(GenomicFeatures::genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene))

allgenes$threesd <- as.numeric(allgenes$gene_id %in% threesdgenes)

threesd <- data.frame(allgenes$threesd)
rownames(threesd) <- allgenes$gene_id
threesd <- as.factor(allgenes$threesd)
names(threesd) <- allgenes$gene_id

genes.3sd <- as.factor(rep(1, length(threesdgenes)))
names(genes.3sd) <- threesdgenes

gene_filter <- function(allScore){
  return(allScore == 1)
}

go3sd <- new("topGOdata",
             ontology = "BP", 
             allGenes = threesd,
             annotationFun = annFUN.gene2GO, 
             geneSelectionFun = gene_filter1,
             gene2GO = gene_GO
)

hippogenes <- genesInTerm(go3sd, "GO:0035329")[[1]]
all <- genes(go3sd)
topGenes <- sigGenes(go3sd)

length(base::intersect(topGenes, threesdgenes))

hippotest.3sd <- new("classicCount", testStatistic = GOFisherTest,
                     name = "fisher",
                     allMembers = all, groupMembers = hippogenes,
                     sigMembers = threesdgenes)

contTable(hippotest.3sd)
runTest(hippotest.3sd)

termStat(go3sd, "GO:0035329")

resultFisher3sd <- runTest(go3sd, algorithm = "classic", statistic = "fisher")

allRes20_3sd <- GenTable(go3sd, classic = resultFisher3sd, ranksOf = "classic", topNodes = 20)

write.csv(allRes20_3sd, file = "../Tables/ds_artselGOtop20_3sd.csv")

allRes50_3sd <- GenTable(go3sd, classic = resultFisher3sd, ranksOf = "classic", topNodes = 50)

write.csv(allRes50_3sd, file = "../Tables/ds_artselGOtop50_3sd.csv")

allRes100_3sd <- GenTable(go3sd, classic = resultFisher3sd, ranksOf = "classic", topNodes = 100)

write.csv(allRes100_3sd, file = "../Tables/ds_artselGOtop100_3sd.csv")

obs.SigHippo3sd <- as.numeric(termStat(go3sd, "GO:0035329")[2])
obs.HippoRatio3sd <- as.numeric(termStat(go3sd, "GO:0035329")[2])/as.numeric(termStat(go3sd, "GO:0035329")[3])

obs.SigNegHippo3sd <- as.numeric(termStat(go3sd, "GO:0035331")[2])
obs.NegHippoRatio3sd <- as.numeric(termStat(go3sd, "GO:0035331")[2])/as.numeric(termStat(go3sd, "GO:0035331")[3])

###################purmutation testing##############
#the goal of this is to ask how common this over rep of ds actually is. 

#For the permutation. 

#I also want to streamline this dataset. This is all the avalible windows in my dataset 
permgenome <- fixed_artsel.5000[,1:2]
#fixing for later. 
permgenome$chr <- paste0("chr", permgenome$chr)
permgenome$Lwindow <- permgenome$window - 2500
permgenome$Rwindow <- permgenome$window + 2500

#using 3sdfirst 
nwindows3sd <- dim(filter(fixed_artsel.5000, dnup >= cutoff))[1]
nwindows2sd <- dim(filter(fixed_artsel.5000, dnup >= cutoff2))[1]

###############Loop over the genome######################
n <- 1000
#make vectors to store 
sig.hippo.count3sd <- rep(NA, n)
ratio.vec3sd <- rep(NA, n)

SigHippo3sd <- rep(NA, n)
HippoRatio3sd <- rep(NA, n)

SigNegHippo3sd <- rep(NA, n)
NegHippoRatio3sd <- rep(NA, n)

#commented out so it doesnt run forever. 
# for(i in 1:n){
# 
# #step 1: sample the genome to get nwidows significant sites
# perm.sig.sites <- permgenome[sample(nrow(permgenome), nwindows3sd),]
# 
# #Assuming these are independent sites and going to make a grange object.
# permRanges <- makeGRangesFromDataFrame(perm.sig.sites,
#                                      keep.extra.columns = FALSE,
#                                      ignore.strand = TRUE,
#                                      start.field = "Lwindow",
#                                      end.field = "Rwindow")
# 
# 
# 
# #getting the genes in the peaks.
# permResults <- subsetByOverlaps(GenomicFeatures::genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene), permRanges)
# 
# permgenes <- permResults$gene_id
# 
# allgenes$perm <- as.numeric(allgenes$gene_id %in% permgenes)
# 
# test <- data.frame(allgenes$perm)
# rownames(test) <- allgenes$gene_id
# 
# test <- as.factor(allgenes$perm)
# names(test) <- allgenes$gene_id
# 
# interesting <- as.factor(rep(1, length(permgenes)))
# names(interesting) <- permgenes
# 
# permGO <- new("topGOdata",
#                 ontology = "BP",
#                 allGenes = test,
#                 annotationFun = annFUN.gene2GO,
#                 geneSelectionFun = gene_filter,
#                 gene2GO = gene_GO
# )
# 
# SigHippo3sd[i] <- as.numeric(termStat(permGO, "GO:0035329")[2])
# HippoRatio3sd[i] <- as.numeric(termStat(permGO, "GO:0035329")[2])/as.numeric(termStat(permGO, "GO:0035329")[3])
# 
# SigNegHippo3sd[i] <- as.numeric(termStat(permGO, "GO:0035331")[2])
# NegHippoRatio3sd[i] <- as.numeric(termStat(permGO, "GO:0035331")[2])/as.numeric(termStat(permGO, "GO:0035331")[3])
# 
# }

#save(SigHippo3sd, HippoRatio3sd, SigNegHippo3sd, NegHippoRatio3sd, obs.SigHippo3sd, obs.HippoRatio3sd, obs.SigNegHippo3sd, obs.NegHippoRatio3sd, file = "../Data/dsGOpermResults3sd.Rda")

load("../Data/dsGOpermResults3sd.Rda")

hist(SigHippo3sd)
abline(v = obs.SigHippo3sd, col = "red")

sum(obs.SigHippo3sd <= SigHippo3sd, na.rm = TRUE)

hist(HippoRatio3sd)
abline(v = obs.HippoRatio3sd, col = "red")

sum(obs.HippoRatio3sd <= HippoRatio3sd, na.rm = TRUE)


hist(SigNegHippo3sd)
abline(v = obs.SigNegHippo3sd, col = "red")

sum(obs.SigNegHippo3sd <= SigNegHippo3sd, na.rm = TRUE)


hist(NegHippoRatio3sd)
abline(v = obs.NegHippoRatio3sd, col = "red")

sum(obs.NegHippoRatio3sd <= NegHippoRatio3sd, na.rm = TRUE)


#drawing the plots quickly in ggplot because it is just easier 
library(cowplot)

perm.df <- data.frame(NegHippoRatio3sd, HippoRatio3sd)

hippoHist <- ggplot(perm.df, aes(x = HippoRatio3sd)) + 
  geom_histogram(binwidth = 0.5) + 
  xlab("Ratio of significant:expected genes") + 
  geom_vline(xintercept = obs.HippoRatio3sd, col = "red") +
  scale_x_continuous(breaks = seq(0, 7, 1))

negHipoHist <- ggplot(perm.df, aes(x = NegHippoRatio3sd)) + 
  geom_histogram(binwidth = 0.5) + 
  xlab("Ratio of significant:expected genes") + 
  geom_vline(xintercept = obs.NegHippoRatio3sd, col = "red") +
  scale_x_continuous(breaks = seq(0, 9, 1))

histPlot <- plot_grid(hippoHist, negHipoHist, nrow = 1, labels = c("Hippo Signaling, GO:0035329", "Negitive regulator of hippo signalling, GO:0035331"), hjust = -0.3)


png("../Figures/GOtermPerm.png", width = 1000, height = 500, units = "px")
histPlot
dev.off()




