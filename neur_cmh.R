#This script is for plotting the cmh test stat across the genome in the wild wing popuation 

library(data.table)
library(tidyverse)

#my plotting functions
source("KP_genomescan_source.R")

#############################################
#data. I already took out the middle cols. 
wildneur.cmh <- fread("../Data/neur_new.cmh")

wildneur.cmh <- wildneur.cmh[,c(1:3, 12)]

colnames(wildneur.cmh) <- c("chr", "pos", "ref", "pval")

#There are NA in the dataset, these are non-informative sites that I will remove. 
length(is.na(wildneur.cmh$pval))

#remove non-informative points. 
wildneur.cmh <- filter(wildneur.cmh, pval != "NaN")

#negative log of p-val for plotting 
wildneur.cmh$logpval <- -log(wildneur.cmh$pval)

#Putting the X first and adding numbers for plotting 
fixed_wildneur_cmh <- chrNumbering(wildneur.cmh)

#going to save this cleaned up for later/other scripts. 
#write_csv(fixed_wildneur_cmh, "wildneur_cmh.csv")

#I also now need to pull out the middle of each chr for plotting. 
chrlabel <- middleChr(fixed_wildneur_cmh)


cmhPlot <- ggplot(data = fixed_wildneur_cmh, aes(x=number, y=logpval, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("-logpval") +
  scale_x_discrete(limits=c(chrlabel),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))

png("../Figures//neur_wild_cmh.png", width=1060,height=412,units="px")
cmhPlot
dev.off()

#Now I want to do a FDR correction for the p-value
fixed_wildneur_cmh$FDRp <- p.adjust(fixed_wildneur_cmh$pval, method = "bonferroni")


hits <- filter(fixed_wildneur_cmh, FDRp <= 0.05)
hits
#short list. 
nrow(hits)

#writing to use for subsetting 
write.table(hits, file = "../Tables/neur_sigCMHloci_missingX.table", row.names = FALSE, quote = FALSE, sep = "\t")

#hits positions only 
loc.only <- dplyr::select(hits, chr, pos)
write.table(loc.only, file = "../Tables/neur_location_hitsOnly.table", row.names = FALSE, quote = FALSE, sep = "\t")

#now to highlight these. 
cmhPlot <- ggplot(data = fixed_wildneur_cmh, aes(x=number, y=logpval, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("-logpval") +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  geom_point(data = hits, aes(x=number, y=logpval), col = "red")+
  scale_x_discrete(limits=c(chrlabel),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))

png("../Figures/neur_wild_cmh_withFDR.png", width=1060,height=412,units="px")
cmhPlot
dev.off()

#now I want to pull out neur specifically 

neur <- (fixed_wildneur_cmh %>% 
         filter(chr == "3R" ) %>%
         filter(pos >= 9020348) %>%
         filter(pos <= 9039471)
)


neur_FDR <- (hits %>% 
             filter(chr == "3R" ) %>%
             filter(pos >= 9020348) %>%
             filter(pos <= 9039471)
)

