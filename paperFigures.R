###This script has the protocol for making the figures for WingShapeBSA paper. 
library(cowplot)
library(magick)
library(gridGraphics)

######
source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R" )
source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R" )
#####
#need some wings to plot against. 
wildwings <- read.csv("../Data/BSA_all_wings.csv")

WingPlot(wildwings[,5:100])



##Figure 1. DGRP allignment and wings
rm(list = ls())

source("DGRPfig.R")

png(file = "../Figures/dgrp_all_corplot_bothlabs.png", width =2000, height = 2000, units = "px",res = 300)
pairs(line.PC[, c(154:156, 98:100)], lower.panel = panel.cor)
dev.off()

#wing with landmarks
landmarks.raw <- image_read("../Data/landmarks.tif")
landmarks.crop <- image_trim(landmarks.raw)
landmarks <- ggdraw() + draw_image(landmarks.crop)

#ds effect 

dseffect.raw <- image_read("../Figures/ds_KDeffectPlot_geomorph_1x.png")
dseffect.crop <- image_trim(dseffect.raw)
dseffect <- ggdraw() + draw_image(dseffect.crop)

dgrp.raw <- image_read("../Figures/dgrp_all_corplot_bothlabs.png")
dgrp.crop <- image_trim(dgrp.raw)
dgrp <- ggdraw() + draw_image(dgrp.crop)


shape.row <- plot_grid(dseffect, landmarks,
                       ncol = 2, nrow = 1, 
                       labels = c("B", "C"), 
                       scale = c(0.8, 0.9),
                       label_size = 12 
                       )

figure1 <-plot_grid(dgrp, shape.row, 
          ncol = 1, nrow = 2, 
          scale = c(1, 1),
          rel_heights = c(2, 1),
          labels = c("A", ""),
          label_size = 12, 
          greedy = FALSE
          )

png("../Figures/fig1_DGRP.png", 
    width =1500, height = 2000, 
    units = "px",res = 300)
figure1
dev.off()


#selvec <- read.csv("../Data/seldict_vectors.csv")
#str(selvec)
#vec2 <- read.csv("../Data/newdict_vectors_uncorrelated_to_ds.csv")
#str(vec2)
#neur <- vec2[5,2:97]

# dseff <- WingEffect( meanshape=as.matrix( colMeans(wildwings[,5:100] ) ), effectplus=as.matrix( selvec[1,3:98]/100 ),
#                      effectminus=as.matrix( selvec[1,3:98]/100 ), winglabel=paste(selvec[1,1]),
#                      scale.factor=1, wingcol=c("black", "blue", "red"), wingframe = FALSE )

###Figure 2. Wild PC figure -> in shape code to make individual populations. 

rm(list = ls())
source("NewShapeCode.R")

#choosing to use cmo as the representitive image for the text of the paper. 
png("../Figures/figure2_cmoExample.png", width =2000, height = 2000, units = "px",res = 300)
pairs( cmo[,207:212], lower.panel=panel.cor )
dev.off()


#Figure 3. ds art sel genome scan, and response to selection. 
rm(list = ls())

#for now I am going to use the fst plot. 


#wing shape and selection 
source("dsreponsetosel.R")
#this is the response to seleciton following mutliple

f.response <- ggplot(female_means, aes(x = gen, y = ds, shape = line, col = rep)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.2) + 
  scale_colour_manual(values=c('black', 'grey46', 'grey57')) +
  geom_line(data = linegen_effect2, col = "red") + 
  theme_classic() + 
  ylab(expression(paste(italic("ds"), " shape score")))+
  xlab("Generation") + 
  theme(legend.position = "none")  + 
  scale_x_continuous(breaks = 1:7)  

#because it is near impossible to combine base plots with ggplots I am going to do what Maria does and make those pngs and then covert into ggobject with magik 


# ds.sel.control.raw <- image_read("../Figures/ds_control_wing_gen1v7selection.png")
# ds.sel.control.crop <- image_trim(ds.sel.control.raw)
# ds.sel.control <- ggdraw() + draw_image(ds.sel.control.crop)
# 
# ds.sel.up.raw <- image_read("../Figures/ds_up_wing_gen1v7selection.png")
# ds.sel.up.crop <- image_trim(ds.sel.up.raw)
# ds.sel.up <- ggdraw() + draw_image(ds.sel.up.crop)
# 
# ds.sel.down.raw <- image_read("../Figures/ds_down_wing_gen1v7selection.png")
# ds.sel.down.crop <- image_trim(ds.sel.down.raw)
# ds.sel.down <- ggdraw() + draw_image(ds.sel.down.crop)
# 
# png("../Figures/ds_c_wing_selEffect.png",width =1000, height = 1000, units = "px", res = 300,bg = "transparent")
# 
# WingEffect(c_reflect, c_diff, c_diff,
#            wingcol=c("black", "black", "red"),
#            scale.factor = 0.75,
#            scale.display = FALSE,
#            wingframe = FALSE,
#            winglabel = "No Selection",
#            winglwd=c(1, 1, 1))
# dev.off()
# 
# 
# png("../Figures/ds_up_wing_selEffect.png",width =1000, height = 1000, units = "px", res = 300,bg = "transparent")
# 
# WingEffect(up_reflect, up_diff, up_diff,
#            wingcol=c("black", "black", "red"),
#            scale.factor = 0.75,
#            scale.display = FALSE,
#            wingframe = FALSE,
#            winglabel = "Up Selection",
#            winglwd=c(1, 1, 1))
# dev.off()
# 
# png("../Figures/ds_dn_wing_selEffect.png",width =1000, height = 1000, units = "px", res = 300,bg = "transparent")
# 
# WingEffect(down_reflect, down_diff, down_diff,
#            wingcol=c("black", "black", "red"),
#            scale.factor = 0.75,
#            scale.display = FALSE,
#            wingframe = FALSE,
#            winglabel = "Down Selection",
#            winglwd=c(1, 1, 1))
# dev.off()

ds.sel.control.raw <- image_read("../Figures/ds_selection_shapechange_CN_F1v7_2x.png")
ds.sel.control.crop <- image_trim(ds.sel.control.raw)
ds.sel.control.crop.r <- image_flop(ds.sel.control.crop)
ds.sel.control <- ggdraw() + draw_image(ds.sel.control.crop.r)

ds.sel.up.raw <- image_read("../Figures/ds_selection_shapechange_UP_F1v7_2x.png")
ds.sel.up.crop <- image_trim(ds.sel.up.raw)
ds.sel.up.crop.r <- image_flop(ds.sel.up.crop)
ds.sel.up <- ggdraw() + draw_image(ds.sel.up.crop.r)

ds.sel.down.raw <- image_read("../Figures/ds_selection_shapechange_DN_F1v7_2x.png")
ds.sel.down.crop <- image_trim(ds.sel.down.raw)
ds.sel.down.crop.r <- image_flop(ds.sel.down.crop)
ds.sel.down <- ggdraw() + draw_image(ds.sel.down.crop.r)



wing.col <- plot_grid(ds.sel.up, ds.sel.control, ds.sel.down,
                      ncol = 1, nrow = 3,
                      scale = 0.75, 
                      label_size = 8, 
                      labels = c("Up Selection", "No Selection", "Down Selection"))


figure_3_top <- plot_grid(f.response, wing.col, 
                      ncol = 2, nrow = 1, scale = c(1, 0.9), 
                      rel_widths = c(1, 0.6))

# png("../Figures/fig3_phenotypicSel.png", width =2000, height = 2000, units = "px",res = 300)
# figure_3_top
# dev.off()

#now for the genomics part! 
source("ds_FINAL_artselgenome.R")

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
  geom_vline(xintercept = ds$number, size = 1, color = 'red', alpha = 0.7) +
  geom_hline(yintercept = cutoff, alpha = 0.5) +
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10), 
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=0.5))

# f.response <- ggplot(female_means, aes(x = gen, y = emc, shape = treat, col = rep)) +
#   geom_point(alpha = 0.5) + 
#   geom_line(alpha = 0.2) + 
#   scale_colour_manual(values=c('black', 'grey46', 'grey57')) +
#   geom_line(data = linegen_effect2, col = "red") + 
#   theme_classic() + 
#   ylab(expression(paste(italic("ds"), " shape score")))+
#   xlab("Generation") + 
#   theme(legend.position = "none")  + 
#   scale_x_continuous(breaks = 1:7)  


figure_3 <- plot_grid(figure_3_top, linescan, 
                      ncol = 1, nrow = 2, scale = c(1, 0.95),
                      rel_heights = c(1, 0.6),
                      labels = c("A", "B"))

png("../Figures/fig3_dsSelection_new.png", width =2000, height = 2000, units = "px",res = 300)
figure_3
dev.off()


#Figure 4. emc art sel and response to selection. 

rm(list = ls())

source("emc_responsetosel.R")
source("emc_FINAL_artselgenome.R")

linescan <- ggplot(data = fixed_artsel.5000, aes(x=number, y=dnup, color=chr, alpha = outlier)) + 
  geom_point(size=0.25, show.legend = F) + 
  scale_alpha_discrete(range = c(0.2,0.7)) +
  theme(panel.background = element_blank()) +
  #scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1)) +
  xlab("Chromosome") +
  ylab(expression(F[ST])) +
  scale_x_discrete(limits=c(chrlabel),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  geom_vline(xintercept = emc$number, size = 1, color = 'purple', alpha = 0.5) +
  geom_vline(xintercept = ds$number, size = 1, color = 'red', alpha = 0.5) + 
  geom_hline(yintercept = cutoff, alpha = 0.5) +
  theme(text = element_text(size=10),
        axis.text.x= element_text(size=8), 
        axis.text.y= element_text(size=8), 
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=0.5))

f.response <- ggplot(female_means, aes(x = gen, y = emc, shape = treat, col = rep)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.2) + 
  scale_colour_manual(values=c('black', 'grey46', 'grey57')) +
  geom_line(data = linegen_effect2[linegen_effect2$Sex =="F",], col = "red") + 
  theme_classic() + 
  ylab(expression(paste(italic("emc"), " shape score")))+
  xlab("Generation") + 
  theme(legend.position = "none")  + 
  scale_x_continuous(breaks = 1:7)  

emc.sel.control.raw <- image_read("../Figures/emc_selection_shapechange_CR_F1v7_10xmag.png")
emc.sel.control.crop <- image_trim(emc.sel.control.raw)
emc.sel.control.crop.r <- image_flop(emc.sel.control.crop)
emc.sel.control <- ggdraw() + draw_image(emc.sel.control.crop.r)

emc.sel.up.raw <- image_read("../Figures/emc_selection_shapechange_UP_F1v7_10xmag.png")
emc.sel.up.crop <- image_trim(emc.sel.up.raw)
emc.sel.up.crop.r <- image_flop(emc.sel.up.crop)
emc.sel.up <- ggdraw() + draw_image(emc.sel.up.crop.r)

emc.sel.down.raw <- image_read("../Figures/emc_selection_shapechange_DN_F1v7_10xmag.png")
emc.sel.down.crop <- image_trim(emc.sel.down.raw)
emc.sel.down.crop.r <- image_flop(emc.sel.down.crop)
emc.sel.down <- ggdraw() + draw_image(emc.sel.down.crop.r)



emc.wing.col <- plot_grid(emc.sel.up, emc.sel.control, emc.sel.down,
                      ncol = 1, nrow = 3,
                      scale = 0.75, 
                      label_size = 8, 
                      labels = c("Up Selection", "No Selection", "Down Selection"))


figure_4_top <- plot_grid(f.response, emc.wing.col, 
                          ncol = 2, nrow = 1, scale = c(1, 0.9), 
                          rel_widths = c(1, 0.6))



figure_4 <- plot_grid(figure_4_top, linescan, 
                      ncol = 1, nrow = 2, scale = c(1, 0.95),
                      rel_heights = c(1, 0.6),
                      labels = c("A", "B"))



png("../Figures/fig4_emcSel_new.png", width =1000, height = 1000, units = "px",res = 200)
figure_4
dev.off()

#Figure 5. Fst for both wild groups. 

rm(list = ls())

source("poolmeans.R")

# png("../Figures/ds_wildpools_shapechange_cmo.png")
# WingPlot(cmo14.left, wingcol="black")
# WingPlot(cmo14.right, wingcol="red", add=T)
# dev.off()

ds.pool.groups.raw <- image_read("../Figures/ds_wildpools_shapechange_cmo_mag2.png")
ds.pool.groups.crop <- image_trim(ds.pool.groups.raw)
ds.pool.groups.crop.r <- image_flop(ds.pool.groups.crop)
ds.pools.groups <- ggdraw() + draw_image(ds.pool.groups.crop.r)

source("ds_ACER_cmh_plotting.R")

cmhPlot <- ggplot(data = ds.cmh.fix, aes(x=number, y=logp, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("-log p-value") +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  geom_point(data = hits, aes(x=number, y=logp), col = "red")+
  scale_x_discrete(limits=c(ds.cmh.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  theme(text = element_text(size=10),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10))


dsplot <- ggplot(data = ds, aes(x=pos, y=logp)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  ylim(0, 20) + 
  #geom_point(data = ds_FDR, aes(x=number, y=logp), col = "red")+
  xlab(expression(paste("postition ", italic("ds")))) +
  ylab("-log p-value") +
  theme(text = element_text(size=8),
        axis.text.x= element_text(size=8), 
        axis.text.y= element_text(size=8))

# ftplot <- ggplot(data = ft, aes(x=pos, y=logpval)) + 
#   geom_point(size=1, show.legend = F, alpha = 0.6) + 
#   theme(panel.background = element_blank()) +
#   ylim(0, 20) +
#   geom_point(data = ft_FDR, aes(x=number, y=logpval), col = "red")+
#   xlab(expression(paste("postition ", italic("ft")))) +
#   ylab("") +
#   theme(text = element_text(size=8),
#         axis.text.x= element_text(size=8), 
#         axis.text.y= element_text(size=8))


#now I want to put these together using cowplot
ds.cmh.shape <- plot_grid(dsplot, ds.pools.groups, 
                          ncol = 2, nrow = 1,
                          scale = c(0.9, 0.9),
                          rel_widths = (c(1.2, 0.8)),
                          labels = c("B", "C"))


png("../Figures/fig5_dsonly.png", 
    width =2000, height = 1500, 
    units = "px",res = 300)

plot_grid(cmhPlot, ds.cmh.shape, ncol = 1, nrow = 2, 
          scale = c(1, 0.90), 
          labels = c("A", ""), 
          rel_heights = (c(0.8, 1)))

dev.off()


###########################SUPPPS#####################

###########genomes-all###############

ds.genome.all.lines <- image_read("../Figures/5000windows_allmarked_3sd_allsites.png")
ds.genome.all.lines.crop <- image_trim(ds.genome.all.lines)
ds.genome.all.good <- ggdraw() + draw_image(ds.genome.all.lines.crop)

emc.genome.all.lines <- image_read("../Figures/emc_5000windows_all_lines_3sd.png")
emc.genome.all.lines.crop <- image_trim(emc.genome.all.lines)
emc.genome.all.good <- ggdraw() + draw_image(emc.genome.all.lines.crop)

png("../Figures/allLines_artSelGenome.png", 
    width =2000, height = 1500, 
    units = "px",res = 300)

plot_grid(ds.genome.all.good, emc.genome.all.good , 
          ncol = 1, nrow = 2, 
          labels = c("A", "B"))
dev.off()


#########Neur Wild###################
rm(list = ls())

source("poolmeans.R")

# png("../Figures/cmo_neurLvRpoolshapes_mag2.png")
# WingPlot(neur.cmo14.left, wingcol="black")
# WingPlot(neur.cmo14.right, wingcol="red", add=T)
# dev.off()

neur.pool.groups.raw <- image_read("../Figures/wildpools_shapechange_cmo_mag2.png")
neur.pool.groups.crop <- image_trim(neur.pool.groups.raw)
neur.pool.groups.crop.r <- image_flop(neur.pool.groups.crop)
neur.pools.groups <- ggdraw() + draw_image(neur.pool.groups.crop.r)

source("neur_ACER_cmh_plotting.R")

neur_cmh_plot <- ggplot(data = neur.cmh.fix, aes(x=number, y=logp, color=chr)) + 
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

neurplot <- ggplot(data = neur, aes(x=pos, y=logp)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("postition") +
  ylab("-log(p-val)") +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))


neur.cmh.shape <- plot_grid(neurplot, neur.pools.groups, 
                          ncol = 2, nrow = 1,
                          scale = c(0.9, 0.9),
                          rel_widths = (c(1.5, 0.75)),
                          labels = c("B", "C"))


png("../Figures/sup_neurWild.png", 
    width =2000, height = 1500, 
    units = "px",res = 300)

plot_grid(neur_cmh_plot, neur.cmh.shape, ncol = 1, nrow = 2, 
          scale = c(0.9, 0.90), 
          labels = c("A", ""), 
          rel_heights = (c(0.8, 1)))

dev.off()


#####Supp Figure - WildWings selection  

fvw12_tails <- proj_tails( vector= dx, wings=as.matrix( wildwings[ wildwings$pop_yr == "fvw12" ,5:100] ), trunk_size=50, ID=wildwings$new_ID[ wildwings$pop_yr == "fvw12"] )
fvw14_tails <- proj_tails( vector= dx, wings=as.matrix( wildwings[ wildwings$pop_yr == "fvw14" ,5:100] ), trunk_size=50, ID=wildwings$new_ID[ wildwings$pop_yr == "fvw14"] )
pho14_tails <- proj_tails( vector= dx, wings=as.matrix( wildwings[ wildwings$pop_yr == "pho14" ,5:100] ), trunk_size=50, ID=wildwings$new_ID[ wildwings$pop_yr == "pho14"] )
cmo14_tails <- proj_tails( vector= dx, wings=as.matrix( wildwings[ wildwings$pop_yr == "cmo14" ,5:100] ), trunk_size=50, ID=wildwings$new_ID[ wildwings$pop_yr == "cmo14"] )




png("selectionpoolsdist.png")
par(mfrow=c(2, 2))

par(mar = c(2,2,2,2))
hist( wildwings$ds_proj[ wildwings$pop_yr == "fvw12"], main="FVW 2012", breaks=20, col="grey", ylab = "" )

abline( v=mean( wildwings$ds_proj[ wildwings$new_ID %in% fvw12_tails$right ] ), col='red', lwd=1.5 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% fvw12_tails$right ], 0.025 ), col='red', lwd=1.5, lty=2 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% fvw12_tails$right ], 0.975 ), col='red', lwd=1.5, lty=2 )

abline( v=mean( wildwings$ds_proj[ wildwings$new_ID %in% fvw12_tails$left ] ), col='red', lwd=1.5 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% fvw12_tails$left ], 0.025 ), col='red', lwd=1.5, lty=2 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% fvw12_tails$left ], 0.975 ), col='red', lwd=1.5, lty=2 )



par(mar = c(2,2,2,2))
hist( wildwings$ds_proj[ wildwings$pop_yr == "fvw14"], main="FVW 2014", breaks=20, col="grey", ylab = "" )

abline( v=mean( wildwings$ds_proj[ wildwings$new_ID %in% fvw14_tails$right ] ), col='red', lwd=1.5 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% fvw14_tails$right ], 0.025 ), col='red', lwd=1.5, lty=2 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% fvw14_tails$right ], 0.975 ), col='red', lwd=1.5, lty=2 )

abline( v=mean( wildwings$ds_proj[ wildwings$new_ID %in% fvw14_tails$left ] ), col='red', lwd=1.5 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% fvw14_tails$left ], 0.025 ), col='red', lwd=1.5, lty=2 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% fvw14_tails$left ], 0.975 ), col='red', lwd=1.5, lty=2 )



par(mar = c(2,2,2,2))
hist( wildwings$ds_proj[ wildwings$pop_yr == "cmo14"], main="CMO 2014", breaks=20, col="grey" , ylab = "")

abline( v=mean( wildwings$ds_proj[ wildwings$new_ID %in% cmo14_tails$right ] ), col='red', lwd=1.5 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% cmo14_tails$right ], 0.025 ), col='red', lwd=1.5, lty=2 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% cmo14_tails$right ], 0.975 ), col='red', lwd=1.5, lty=2 )

abline( v=mean( wildwings$ds_proj[ wildwings$new_ID %in% cmo14_tails$left ] ), col='red', lwd=1.5 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% cmo14_tails$left ], 0.025 ), col='red', lwd=1.5, lty=2 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% cmo14_tails$left ], 0.975 ), col='red', lwd=1.5, lty=2 )



par(mar = c(2,2,2,2))
hist( wildwings$ds_proj[ wildwings$pop_yr == "pho14"], main="PHO 2014", breaks=20, col="grey", ylab = "" )

abline( v=mean( wildwings$ds_proj[ wildwings$new_ID %in% pho14_tails$right ] ), col='red', lwd=1.5 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% pho14_tails$right ], 0.025 ), col='red', lwd=1.5, lty=2 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% pho14_tails$right ], 0.975 ), col='red', lwd=1.5, lty=2 )

abline( v=mean( wildwings$ds_proj[ wildwings$new_ID %in% pho14_tails$left ] ), col='red', lwd=1.5 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% pho14_tails$left ], 0.025 ), col='red', lwd=1.5, lty=2 )
abline( v=quantile( wildwings$ds_proj[ wildwings$new_ID %in% pho14_tails$left ], 0.975 ), col='red', lwd=1.5, lty=2 )

dev.off()



########SUPP STUFF I NEED TO CLEAN #############

pho.raw <- image_read("../Figures/dsneur_wildwingsproj_pho_PC1to3.png")
pho.crop <- image_trim(pho.raw)
pho <- ggdraw() + draw_image(pho.crop)
cmo.raw <- image_read("../Figures/dsneur_wildwingsproj_cmo_PC1to3.png")
cmo.crop <- image_trim(cmo.raw)
cmo <- ggdraw() + draw_image(cmo.crop)
f12.raw <- image_read("../Figures/dsneur_wildwingsproj_f12_PC1to3.png")
f12.crop <- image_trim(f12.raw)
f12 <- ggdraw() + draw_image(f12.crop)
f14.raw <- image_read("../Figures/dsneur_wildwingsproj_f14_PC1to3.png")
f14.crop <- image_trim(f14.raw)
f14 <- ggdraw() + draw_image(f14.crop)

pop.corr <- plot_grid(pho, cmo, f12, f14, labels = c("PHO", "CMO", "FVW12", "FVW14"), label_size = 8)

png("../Figures/All_populationCorr.png")
pop.corr
dev.off()


f14all.raw <- image_read("../Figures/dsneur_wildwingsproj_bothSex_f14_PC1to3.png")
f14all.crop <- image_trim(f14all.raw)
f14all <- ggdraw() + draw_image(f14all.crop)

#####################

bothsex <- plot_grid(f14, f14all, labels = c("A", "B"), label_size = 10)


png("../Figures/corrPlot_females_supp.png")
bothsex
dev.off()


######Shape change supplemental######

dseffect.raw <- image_read('../Figures/ds_KDeffectPlot_geomorph_1x.png')
dseffect.crop <- image_trim(dseffect.raw)
dseffect <- ggdraw() + draw_image(dseffect.crop)


emceffect.raw <- image_read("../Figures/emc_KDeffectPlot_geomorph_10x.png")
emceffect.crop <- image_trim(emceffect.raw)
emceffect <- ggdraw() + draw_image(emceffect.crop)

neureffect.raw <- image_read("../Figures/neur_KDeffectPlot_geomorph_2x.png")
neureffect.crop <- image_trim(neureffect.raw)
neureffect <- ggdraw() + draw_image(neureffect.crop)

png("../Figures/../Figures/dictVecKD_normalMag.png", width =720, height = 300, units = "px",res = 300)
plot_grid(dseffect, emceffect, neureffect, ncol = 3, nrow = 1, 
         scale = c(0.9, 0.9, 0.9), 
          labels = c("ds 1x", "emc 10x", "neur 2x"),
          label_size = 8,
          rel_heights = (c(0.9, 0.9, 0.9)) )

dev.off()


#######################################
Fdsbur.raw <- image_read("../Figures/dsSel_F_wingBlur.png")
Fdsbur.crop <- image_trim(Fdsbur.raw)
Fdsbur.crop.r <- image_flop(Fdsbur.crop)
Fdsbur <- ggdraw() + draw_image(Fdsbur.crop.r)


Mdsbur.raw <- image_read("../Figures/dsSel_M_wingBlur.png")
Mdsbur.crop <- image_trim(Mdsbur.raw)
Mdsbur.crop.r <- image_flop(Mdsbur.crop)
Mdsbur <- ggdraw() + draw_image(Mdsbur.crop.r)

png("../Figures/dsSelectionBlurs.png", height = 250, width = 500)
plot_grid(Fdsbur, Mdsbur, labels = c("A", "B"))
dev.off()
