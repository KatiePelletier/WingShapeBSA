#plotting by population fst plots for wild pools (left vs right)
library(data.table)
library(tidyverse)

#my plotting functions
source("KP_genomescan_source.R")

ds.wild <- fread("../Data/wild_lvr_nosync_100window.fst")


ds.wild <- populaionFst_cleanup(ds.wild, x = c('CMO', 'FVW13', 'FVW14', 'PHO'))

#Reordering and numbering the chr for plotting.
fixed_ds.wild <- chrNumbering(ds.wild)

chrlabel.ds <- middleChr(fixed_ds.wild)

#Quick Plots 
#ID suggested adding either the mean fst to each chr or a smoothing spline and increasing alpha 
cmo.xMean <- mean(fixed_ds.wild[chr == "X", CMO])
cmo.2LMean <- mean(fixed_ds.wild[chr == "2L", CMO])
cmo.2RMean <- mean(fixed_ds.wild[chr == "2R", CMO])
cmo.3LMean <- mean(fixed_ds.wild[chr == "3L", CMO])
cmo.3RMean <- mean(fixed_ds.wild[chr == "3R", CMO])
cmo.4Mean <- mean(fixed_ds.wild[chr == "4", CMO])


#can resuse these all for all of them (because the numbering is actually common)
cmo.xstart <- min(fixed_ds.wild[chr == "X", number])
cmo.2Lstart <- min(fixed_ds.wild[chr == "2L", number])
cmo.2Rstart <- min(fixed_ds.wild[chr == "2R", number])
cmo.3Lstart <- min(fixed_ds.wild[chr == "3L", number])
cmo.3Rstart <- min(fixed_ds.wild[chr == "3R", number])
cmo.4start <- min(fixed_ds.wild[chr == "4", number])


CMO <- ggplot(data = fixed_ds.wild, aes(x=number, y=CMO, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.2) +
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("meanFst") +
  scale_x_discrete(limits=c(chrlabel.ds),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  geom_segment(x = cmo.xstart, xend = cmo.2Lstart, y = cmo.xMean, yend = cmo.xMean , col = "red") + 
  geom_segment(x = cmo.2Lstart, xend = cmo.2Rstart, y = cmo.2LMean, yend = cmo.2LMean , col = "red") + 
  geom_segment(x = cmo.2Rstart, xend = cmo.3Lstart, y = cmo.2RMean, yend = cmo.2RMean , col = "red") + 
  geom_segment(x = cmo.3Lstart, xend = cmo.3Rstart, y = cmo.3LMean, yend = cmo.3LMean , col = "red") + 
  geom_segment(x = cmo.3Rstart, xend = cmo.4start, y = cmo.3RMean, yend = cmo.3RMean , col = "red") + 
  geom_segment(x = cmo.4start, xend = max(fixed_ds.wild$number), y = cmo.4Mean, yend = cmo.4Mean , col = "red") + 
  ylim(0,0.75) + 
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12),
        axis.text.y= element_text(size=12),
        panel.border = element_rect(colour = "black",
                                    fill=NA, size=0.5))

f13.xMean <- mean(fixed_ds.wild[chr == "X", FVW13])
f13.2LMean <- mean(fixed_ds.wild[chr == "2L", FVW13])
f13.2RMean <- mean(fixed_ds.wild[chr == "2R", FVW13])
f13.3LMean <- mean(fixed_ds.wild[chr == "3L", FVW13])
f13.3RMean <- mean(fixed_ds.wild[chr == "3R", FVW13])
f13.4Mean <- mean(fixed_ds.wild[chr == "4", FVW13])

FVW13 <- ggplot(data = fixed_ds.wild, aes(x=number, y=FVW13, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.2) +
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("meanFst") +
  scale_x_discrete(limits=c(chrlabel.ds),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  ylim(0,0.75) + 
  geom_segment(x = cmo.xstart, xend = cmo.2Lstart, y = f13.xMean, yend = f13.xMean , col = "red") + 
  geom_segment(x = cmo.2Lstart, xend = cmo.2Rstart, y = f13.2LMean, yend = f13.2LMean , col = "red") + 
  geom_segment(x = cmo.2Rstart, xend = cmo.3Lstart, y = f13.2RMean, yend = f13.2RMean , col = "red") + 
  geom_segment(x = cmo.3Lstart, xend = cmo.3Rstart, y = f13.3LMean, yend = f13.3LMean , col = "red") + 
  geom_segment(x = cmo.3Rstart, xend = cmo.4start, y = f13.3RMean, yend = f13.3RMean , col = "red") + 
  geom_segment(x = cmo.4start, xend = max(fixed_ds.wild$number), y = f13.4Mean, yend = f13.4Mean , col = "red") + 
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12),
        axis.text.y= element_text(size=12),
        panel.border = element_rect(colour = "black",
                                    fill=NA, size=0.5))

f14.xMean <- mean(fixed_ds.wild[chr == "X", FVW14])
f14.2LMean <- mean(fixed_ds.wild[chr == "2L", FVW14])
f14.2RMean <- mean(fixed_ds.wild[chr == "2R", FVW14])
f14.3LMean <- mean(fixed_ds.wild[chr == "3L", FVW14])
f14.3RMean <- mean(fixed_ds.wild[chr == "3R", FVW14])
f14.4Mean <- mean(fixed_ds.wild[chr == "4", FVW14])

FVW14 <- ggplot(data = fixed_ds.wild, aes(x=number, y=FVW14, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.2) +
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("meanFst") +
  scale_x_discrete(limits=c(chrlabel.ds),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  geom_segment(x = cmo.xstart, xend = cmo.2Lstart, y = f14.xMean, yend = f14.xMean , col = "red") + 
  geom_segment(x = cmo.2Lstart, xend = cmo.2Rstart, y = f14.2LMean, yend = f14.2LMean , col = "red") + 
  geom_segment(x = cmo.2Rstart, xend = cmo.3Lstart, y = f14.2RMean, yend = f14.2RMean , col = "red") + 
  geom_segment(x = cmo.3Lstart, xend = cmo.3Rstart, y = f14.3LMean, yend = f14.3LMean , col = "red") + 
  geom_segment(x = cmo.3Rstart, xend = cmo.4start, y = f14.3RMean, yend = f14.3RMean , col = "red") + 
  geom_segment(x = cmo.4start, xend = max(fixed_ds.wild$number), y = f14.4Mean, yend = f14.4Mean , col = "red") + 
  ylim(0,0.75) + 
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12),
        axis.text.y= element_text(size=12),
        panel.border = element_rect(colour = "black",
                                    fill=NA, size=0.5))


pho.xMean <- mean(fixed_ds.wild[chr == "X", PHO])
pho.2LMean <- mean(fixed_ds.wild[chr == "2L", PHO])
pho.2RMean <- mean(fixed_ds.wild[chr == "2R", PHO])
pho.3LMean <- mean(fixed_ds.wild[chr == "3L", PHO])
pho.3RMean <- mean(fixed_ds.wild[chr == "3R", PHO])
pho.4Mean <- mean(fixed_ds.wild[chr == "4", PHO])

PHO <- ggplot(data = fixed_ds.wild, aes(x=number, y=PHO, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.2) +
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("meanFst") +
  scale_x_discrete(limits=c(chrlabel.ds),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  ylim(0,0.75) + 
  geom_segment(x = cmo.xstart, xend = cmo.2Lstart, y = pho.xMean, yend = pho.xMean , col = "red") + 
  geom_segment(x = cmo.2Lstart, xend = cmo.2Rstart, y = pho.2LMean, yend = pho.2LMean , col = "red") + 
  geom_segment(x = cmo.2Rstart, xend = cmo.3Lstart, y = pho.2RMean, yend = pho.2RMean , col = "red") + 
  geom_segment(x = cmo.3Lstart, xend = cmo.3Rstart, y = pho.3LMean, yend = pho.3LMean , col = "red") + 
  geom_segment(x = cmo.3Rstart, xend = cmo.4start, y = pho.3RMean, yend = pho.3RMean , col = "red") + 
  geom_segment(x = cmo.4start, xend = max(fixed_ds.wild$number), y = pho.4Mean, yend = pho.4Mean , col = "red") + 
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12),
        axis.text.y= element_text(size=12),
        panel.border = element_rect(colour = "black",
                                    fill=NA, size=0.5))

library(cowplot)

png("../Figures/wild_ds_byPoolFst.png")
plot_grid(FVW13, FVW14, CMO, PHO, nrow = 2, labels = c("FVW13", "FVW14", "CMO", "PHO"))
dev.off()

mean(fixed_ds.wild$CMO)
mean(fixed_ds.wild$PHO)
mean(fixed_ds.wild$FVW13)
mean(fixed_ds.wild$FVW14)


###Maybe I want to look at some of those higher windows? I dont know? The couple hits in fvw13 look sketch.


#taking some of the high hits in fvw13 quickly. 
#these all have pretty low coverage. 
f13.top <- filter(fixed_ds.wild, FVW13 >= 0.45)
dim(f13.top)



