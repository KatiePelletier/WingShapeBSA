#This is looking at the Fst between popualtions. 

library(data.table)
library(tidyverse)
library(cowplot)

#my plotting functions
source("KP_genomescan_source.R")

#This has all the populations calculated in 1000 bp windows. 
# populations in order: 1) CMO 2) F12 3) F13 4) PHO 
wild1000 <- fread("../Data/allpools_1000window.fst")

wild1000 <- populaionFst_cleanup(wild1000, x = c('cmo.f12', 'cmo.f13', 'cmo.pho', 'f12.f13', 'f12.pho', 'f13.pho'))

#Reordering and numbering the chr for plotting. 
fixed_wild1000 <- chrNumbering(wild1000)

chrlabel1000 <- middleChr(fixed_wild1000)

#plotting 
cmo.f12 <- ggplot(data = fixed_wild1000, aes(x=number, y=cmo.f12, color=chr)) + 
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


cmo.f13 <- ggplot(data = fixed_wild1000, aes(x=number, y=cmo.f13, color=chr)) + 
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

cmo.pho <- ggplot(data = fixed_wild1000, aes(x=number, y=cmo.pho, color=chr)) + 
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

f12.f13 <- ggplot(data = fixed_wild1000, aes(x=number, y=f12.f13, color=chr)) + 
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


f12.pho <- ggplot(data = fixed_wild1000, aes(x=number, y=f12.pho, color=chr)) + 
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

f13.pho <- ggplot(data = fixed_wild1000, aes(x=number, y=f13.pho, color=chr)) + 
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

quartz()
png('../Output/betweenpop_fst1000window.png', height = 1000, width = 2000 )
plot_grid(cmo.f12, cmo.f13, cmo.pho, f12.f13, f12.pho, f13.pho, labels = c('CMO F12', 'CMO F13', 'CMO PHO', 'F12 F13', 'F12 PHO', 'F13 PHO' ), label_x = 0)
dev.off()

#I also just want to get average fst for each pair 
#super low maybe f12 is a little more diffrent?
avg <- colMeans(fixed_wild1000[,6:11])

