library(data.table)
library(tidyverse)

source('KP_genomescan_source.R')

neur.cmh <- fread("../Data/neur_wild_ACER_zeroGen_fdr.csv")

#reordering the genome so the X is first and adding numbers for plotting. 

neur.cmh.fix <- chrNumbering(neur.cmh)

neur.cmh.middle <- middleChr(neur.cmh)

#need -logp for plotting.

neur.cmh.fix$logp <- -log(neur.cmh.fix$pval)


#filtering the sig sites for plotting. 
hits <- filter(neur.cmh.fix, padj <= 0.05)
hits
#5 sites. 
nrow(hits)


cmhPlot <- ggplot(data = neur.cmh.fix, aes(x=number, y=logp, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("-log(p-val)") +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  geom_point(data = hits, aes(x=number, y=logp), col = "red")+
  scale_x_discrete(limits=c(neur.cmh.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))

png("../Figures/neur_wild_cmh_ACER_zeroGen_withFDR.png", width=1060,height=412,units="px")
cmhPlot
dev.off()

#Now subsetting out neur alone. 

neur <- (neur.cmh.fix %>% 
           filter(chr == "3R" ) %>%
           filter(pos >= 9020348) %>%
           filter(pos <= 9039471)
)


neurplot <- ggplot(data = neur, aes(x=pos, y=logp)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("postition") +
  ylab("-log(p-val)") +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))

png("../Figures/wild_neur_ACRER_zeroGen_cmh.png", width=1060,height=412,units="px")
neurplot
dev.off()



