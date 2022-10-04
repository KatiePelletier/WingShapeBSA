library(data.table)
library(tidyverse)

source('KP_genomescan_source.R')

ds.cmh <- fread("../Data/subsample75cov_correctPoolSize_dsWild_ACER.csv")

#reordering the genome so the X is first and adding numbers for plotting. 

ds.cmh.fix <- chrNumbering(ds.cmh)

ds.cmh.middle <- middleChr(ds.cmh)

#need -logp for plotting.

ds.cmh.fix$logp <- -log(ds.cmh.fix$pval)


#filtering the sig sites for plotting. 
hits <- filter(ds.cmh.fix, padj <= 0.05)
hits
#More sites than the 15 with all the data. 
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

png("../Figures/ds_wild_cmh_ACER_downsample75x_withFDR_rightPopSize.png", width=1060,height=412,units="px")
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

png("../Figures/wild_ds_ACRER_downsample75x_cmh_betterPoolSize.png", width=1060,height=412,units="px")
dsplot
dev.off()

#printing hit sites for annotation
dim(hits)

write_delim(hits, file = "../Data//wild_ds_ACER_downsample75x_cmh_hits.csv", delim = "\t")
