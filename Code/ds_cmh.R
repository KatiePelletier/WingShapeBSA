#This script is for plotting the cmh test stat across the genome in the wild wing popuation 

library(data.table)
library(tidyverse)

#my plotting functions
source("KP_genomescan_source.R")

#############################################
#data. I already took out the middle cols. 
wildds.cmh <- fread("../Data/wild_nomerge_tabs.cmh")


wildds.cmh <- wildds.cmh[,c(1:3, 16)]

colnames(wildds.cmh) <- c("chr", "pos", "ref", "pval")

#There are NA in the dataset, these are non-informative sites that I will remove. 
length(is.na(wildds.cmh$pval))

#remove non-informative points. 
wildds.cmh <- filter(wildds.cmh, pval != "NaN")

#negative log of p-val for plotting 
wildds.cmh$logpval <- -log(wildds.cmh$pval)


#Putting the X first and adding numbers for plotting 
fixed_wildds_cmh <- chrNumbering(wildds.cmh)

#going to save this cleaned up for later/other scripts. 
#write_csv(fixed_wildds_cmh, "wildds_cmh.csv")

#I also now need to pull out the middle of each chr for plotting. 
chrlabel <- middleChr(fixed_wildds_cmh)


cmhPlot <- ggplot(data = fixed_wildds_cmh, aes(x=number, y=logpval, color=chr)) + 
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

png("../Figures/ds_wild_cmh.png", width=1060,height=412,units="px")
cmhPlot
dev.off()

#Now I want to do a FDR correction for the p-value
fixed_wildds_cmh$FDRp <- p.adjust(fixed_wildds_cmh$pval, method = "bonferroni")


hits <- filter(fixed_wildds_cmh, FDRp <= 0.05)
hits
#short list. 
nrow(hits)

#writing to use for subsetting 
write.table(hits, file = "../Tables/ds_sigCMHloci.table", row.names = FALSE, quote = FALSE, sep ="\t")

#hits positions only 
loc.only <- dplyr::select(hits, chr, pos)
write.table(loc.only, file = "../Tables/ds_location_hitsOnly.table", row.names = FALSE, quote = FALSE, sep ="\t") 



#now to highlight these. 
cmhPlot <- ggplot(data = fixed_wildds_cmh, aes(x=number, y=logpval, color=chr)) + 
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

png("../Figures/ds_wild_cmh_withFDR.png", width=1060,height=412,units="px")
cmhPlot
dev.off()

#Looking specifically at ds. 
ds <- (fixed_wildds_cmh %>% 
         filter(chr == "2L" ) %>%
         filter(pos >= 640013) %>%
         filter(pos <= 714983)
)


ds_FDR <- (hits %>% 
             filter(chr == "2L" ) %>%
             filter(pos >= 640013) %>%
             filter(pos <= 714983)
)

dsplot <- ggplot(data = ds, aes(x=pos, y=logpval)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("postition") +
  ylab("-logpval") +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))

png("../Figures/wild_ds_cmh.png", width=1060,height=412,units="px")
dsplot
dev.off()

#Another gene?  
ft <- (fixed_wildds_cmh %>% 
         filter(chr == "2L" ) %>%
         filter(pos >= 4198404) %>%
         filter(pos <= 4221796)
)


ft_FDR <- (hits %>% 
             filter(chr == "2L" ) %>%
             filter(pos >= 4198404) %>%
             filter(pos <= 4221796)
)

#or emc? 
emc <- (fixed_wildds_cmh %>% 
         filter(chr == "3L" ) %>%
         filter(pos >= 749400) %>%
         filter(pos <= 753492)
)

emc_FDR <- (hits %>% 
             filter(chr == "2L" ) %>%
             filter(pos >= 749400) %>%
             filter(pos <= 753492)
)




