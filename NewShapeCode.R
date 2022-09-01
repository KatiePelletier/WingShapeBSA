####Wild wings shape stuff for the paper and making final figures. 
library(tidyverse)
library(emmeans)

########################
source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R" )
source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R" )
###########################

##### Edited from the WRP_functions file to round better. 
panel.cor <- function(x, y) {
  r <- cor(x, y)
  par( usr =c(0, 1, 0, 1))
  Cor <- formatC(c(r, 0.123456789), format = "f", digits=2)[1]
  text(0.5, 0.5, paste(Cor), cex=1.5)
}


projFunction <- function(x, y) {
  scalarProj <- (x %*% y) / norm(y, type = "2")
  return(as.numeric(scalarProj))
}

######

#reading in the file with all wings superimposed. 
#this data has already been cleaned up and collapsed. 
#still not 100% convinced that these are all in the same eigenspace? This could just be from Will's code (although there was never a line to write this file either?). Talk to Ian about what he thinks. Getting the splines is actually a GIANT headache. 
wildwings <- read.csv("../Data/BSA_all_wings.csv")

#this is a wing. 
WingPlot(PrcCrds=colMeans(wildwings[,5:100] ) )

str(wildwings)

hist(wildwings$CSize)

plot( wildwings$CSize ~ wildwings$sex )


#Had to leave the QC checks in other code. It uses other sets of Will's wings to build M -> F and sim -> mel vectors for LDA. 
#These did look like we had reliable classifications. 

#Projection function 
proj_tails <- function( vector, wings=as.matrix( wildwings[,5:100] ), trunk_size=100, ID=wildwings$new_ID ) {
  proj_vec <- wings %*% vector
  n_wings <- nrow( wings )
  proj_tail_right <- ID[ order( proj_vec ) ][ 1:trunk_size ]
  proj_tail_left <- ID[ order( proj_vec ) ][ ((n_wings+1) - trunk_size):n_wings ]
  list( right=proj_tail_right, left=proj_tail_left, both=c( proj_tail_right, proj_tail_left ) )
}

#The ds projection has already been done in this data set but again to check myself 

selvec <- read.csv( "../Data/seldict_vectors.csv" )
str( selvec )

WingPlot(as.matrix(colMeans(wildwings[,5:100])))

#Houle wings are way bigger (specifically, 100x bigger) than ours. This plot looks funny. Mean center data?
WingEffect( meanshape=as.matrix( colMeans( wildwings[,5:100] ) ), effectplus=as.matrix( selvec[1,3:98]/100 ),
            effectminus=as.matrix( selvec[1,3:98]/100), winglabel=paste(selvec[1,1]),
            scale.factor=1, wingcol=c("black", "blue", "red") )
#The other vecs (maybe for a supp??)
#add later 


#making the actual vectors. 
ds_vec <- as.matrix(selvec[1,3:98])
emc_vec <- as.matrix(selvec[2,3:98])

##Also need Neur vec####
vec2 <- read.csv("../Data/newdict_vectors_uncorrelated_to_ds.csv")
neur <- t(vec2[5,2:97])

#The effect vector from that experiment. Again ds looks crazy. 
WingEffect( as.matrix(selvec[3,3:98])/100, ds_vec/100, ds_vec/100, scale.factor=1)
WingEffect( as.matrix(selvec[3,3:98])/100, emc_vec/100, emc_vec/100, scale.factor=2)

WingEffect( as.matrix(selvec[3,3:98])/100, neur/100, neur/100, scale.factor=1)

#magnitude of ds effect is 10x that of emc effect. 
norm(ds_vec, type = "2") #5.542832
norm(emc_vec, type = "2") #0.4469862
norm(neur, type = "2") #2.847019 why so giant? 


#####F14with females plot for the paper####
#I am very sorry for the shitty names in here. 
with(wildwings, table(pop, sex, block))

f14all <- (wildwings %>%
             filter(pop == "fvw") %>%
             filter(block == 2))

f14all.size <- lm( as.matrix(f14all[,5:100] ) ~ f14all$CSize + f14all$sex)$resid

f14all_all <- data.frame(f14all, f14all.size)

#The projection function is below. I am so sorry 

f14all_all$ds <- projFunction(as.matrix(f14all.size), t(ds_vec))
f14all_all$emc <- projFunction(as.matrix(f14all.size), t(emc_vec))
f14all_all$neur <- projFunction(as.matrix(f14all.size), as.matrix(neur))

f14all_all <- cbind(f14all_all, prcomp(f14all_all[,111:206])$x[,1:57])

png("../Figures/dsneur_wildwingsproj_bothSex_f14_PC1to3.png")
pairs( f14all_all[,207:212], lower.panel=panel.cor )
dev.off()

##################boys only################################
boysonly <- filter(wildwings, sex == "M")

#Now I want to model out CS to account for allometry 
boyssize <- lm( as.matrix(boysonly[,5:100] ) ~ boysonly$CSize)$resid

boys_all <- data.frame(boysonly, boyssize)
str(boys_all)
names(boys_all)


#projection 

boys_all$ds <- projFunction(as.matrix(boyssize), t(ds_vec))


#checking against Will's 
#close but not perfect... ok.... so something is up here.
#Will included the F in the model so data was a little different (did not subset boys)
cor(boys_all$ds_proj, boys_all$ds)

boys_all$emc <- projFunction(as.matrix(boyssize), t(emc_vec))

# boysonly$ct <- projFunction(as.matrix(boyssize), 
                              #t(( as.matrix( ct_vectors[1,2:97] ))

# boysonly$ptc <- projFunction(as.matrix(boyssize), 
#                             t( as.matrix(ptc_vectors[1,2:97]))

# boysonly$vg <- projFunction(as.matrix(boyssize), 
#                             t(as.matrix(newvecs[1,2:97] )))

# boysonly$yki <- projFunction(as.matrix(boyssize), 
#                             t(as.matrix(newvecs[2,2:97])))

# boysonly$hh <- projFunction(as.matrix(boyssize), 
#                             t(as.matrix(newvecs[3,2:97])))

# boysonly$vn <- 

boys_all$neur <- projFunction(as.matrix(boyssize), as.matrix(neur))

# boysonly$bx <- projFunction(as.matrix(boyssize), 
#                             t(as.matrix(newvecs[6,2:97])))

# boysonly$sd <- projFunction(as.matrix(boyssize), 
#                             t(sd_vector)



################some extra paper plots###########################

#what a dumb way to do this. Not sure why I worte it like this???
boys_all$PC1 <- prcomp(boyssize )$x[,1]
boys_all$PC2 <- prcomp( boyssize )$x[,2]
boys_all$PC3 <- prcomp( boyssize )$x[,3]
boys_all$PC4 <- prcomp( boyssize )$x[,4]
boys_all$PC5 <- prcomp( boyssize )$x[,5]


#need to make fvw12 -> fvw13 to make everything match.
boys_all$pop_yr <- gsub("fvw12", "fvw13", boys_all$pop_yr)

pc12 <- ggplot(boys_all, aes(x = PC1, y = PC2, col = pop_yr)) + 
  geom_point(alpha = 0.3) + 
  theme(legend.position="none") + 
  theme(text = element_text(size = 20))
  

pc34 <- ggplot(boys_all, aes(x = PC3, y = PC4, col = pop_yr)) + 
  geom_point(alpha = 0.3) + 
  labs(col = "Population") + 
  theme(legend.position="none") +
  theme(text = element_text(size = 20))

forkey <- ggplot(boys_all, aes(x = PC3, y = PC4, col = pop_yr)) + 
  geom_point() + 
  labs(col = "Population") + 
  theme(text = element_text(size = 20))

key <- get_legend(forkey + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") 
)

#additionally want to do this with size left in (not modeling out as above). This will help (kinda) to compare allometries 

#as a quick little check for me. 
#Pho is a little smaller than the rest but overall not crazy

png("../Figures/wildsize_density.png", width = 480, height = 480, units = "px")
ggplot(boys_all, aes(x = CSize, fill = pop_yr)) + 
  geom_density(alpha = 0.5) + 
  xlab("Centroid Size") + 
  labs(fill = "Population") + 
  theme(text = element_text(size = 20))
dev.off()

#and a quick model 

sizemod <- lm(CSize ~ pop_yr, data = boys_all)
#these numbers are a little tough because it is in CS and not in mm^2.Pretty much shows what you think. FVW collections are the same, PHO is diffrent. Possibly an effect of time of year/micro enviormental stuff. I don't think this is particularly interesting
summary(sizemod)  

#This is not at all a suprising result. 
anova(sizemod)

#fitted val
plot(emmeans(sizemod, "pop_yr"))

#pairs.because I would rather a computer adds for me. 
pairs(emmeans(sizemod, "pop_yr"))


#all I really need
withSize_pca <- prcomp(boys_all[,5:100])$x[,1:4]

size_pc12 <- ggplot(boys_all, aes(x = withSize_pca[,1], y = withSize_pca[,2], col = pop_yr)) + 
  geom_point(alpha = 0.3) + 
  theme(legend.position="none") + 
  xlab("PC1") + 
  ylab("PC2") + 
  theme(text = element_text(size = 20))

size_pc34 <- ggplot(boys_all, aes(x = withSize_pca[,3], y = withSize_pca[,4], col = pop_yr)) + 
  geom_point(alpha = 0.3) + 
  theme(legend.position="none") + 
  xlab("PC3") + 
  ylab("PC4") + 
  theme(text = element_text(size = 20))


#library(cowplot)

size_pan <- plot_grid(size_pc12, size_pc34, labels = c("A", "B"))
no_size_pan <- plot_grid(pc12, pc34, labels = c("C", "D"))


png("../Figures/wildPop_PCA.png", width = 700, height = 700, units = "px")
plot_grid(size_pan, no_size_pan, key, ncol = 1, rel_widths = c(1,1,0.05))
dev.off()
  
#Ian also recomended plotting the ds and neur projections

png("../Figures/wild_proj_scatterplot.png", width = 480, height = 480, units = "px")
ggplot(boys_all, aes(x = ds, y = neur, col = pop_yr)) + 
  geom_point(alpha = 0.4) + 
  labs(col = "Population") + 
  theme(text = element_text(size = 20)) + 
  xlab(expression(paste(italic("ds"), " shape vector"))) + ylab(expression(paste(italic("neur"), " shape vector")))

dev.off()

cor(boys_all$ds, boys_all$neur)


#the final plot Ian recomended was taking the line means and doing a PCA/PCoA and then projecting the data onto those PCs. Should look really similar to the lines themselves. 

#using the size corrected shape residuals. 
m_pho <- colMeans(boys_all[boys_all$pop_yr == "pho14", 116:211])
m_cmo <- colMeans(boys_all[boys_all$pop_yr == "cmo14", 116:211])
m_fvw13 <- colMeans(boys_all[boys_all$pop_yr == "fvw13", 116:211])
m_fvw14 <- colMeans(boys_all[boys_all$pop_yr == "fvw14", 116:211])

#binding them together. 
m_pops <- rbind(m_pho, m_cmo, m_fvw13, m_fvw14)
dim(m_pops)

# Euclidian Distances between shape means. Should match from geomoph estimates. 
dist(m_pops)
#PCA
m_pca <- prcomp(m_pops)
summary(m_pca)

#use these for projection
m_pca$rotation

pop_yr <-rownames(m_pops) 
all_means_pca <- data.frame(pop_yr, m_pops, m_pca$x)

#truly useless. but looks like you would expect from PDs. 


ggplot(all_means_pca, aes(x = PC1, y = PC2, col = pop_yr)) + 
  geom_point() 


means_PC1 <- as.numeric(m_pca$rotation[,1])
means_PC2 <- as.numeric(m_pca$rotation[,2])

#projecting wild wings onto these vectors. 
boys_all$mPC1proj <- projFunction(as.matrix(boyssize), means_PC1)
boys_all$mPC2proj <- projFunction(as.matrix(boyssize), means_PC2)

png("../Figures/wild_PCsproj_scatterplot.png", width = 480, height = 480, units = "px")
ggplot(boys_all, aes(x = mPC1proj, y = mPC2proj, col = pop_yr)) + 
  geom_point(alpha = 0.2) + 
  labs(col = "Population") + 
  theme(text = element_text(size = 20)) + 
  xlab("group means PC1 projection") + 
  ylab("group means PC2 projection")

dev.off()


#################################################################

#going to split populations to look at PCs and plot 

pho <- filter(boys_all, pop_yr == "pho14")
pho <- cbind(pho, prcomp(pho[,111:206])$x[,1:57])

png("../Figures/dsneur_wildwingsproj_pho_PC1to3.png")
pairs( pho[,207:212], lower.panel=panel.cor )
dev.off()

cmo <- filter(boys_all, pop_yr == "cmo14")
cmo <- cbind(cmo, prcomp(cmo[,111:206])$x[,1:57])

png("../Figures/dsneur_wildwingsproj_cmo_PC1to3.png")
pairs( cmo[,207:212], lower.panel=panel.cor )
dev.off()


cmo.cor.plot <- pairs( cmo[,207:212], lower.panel=panel.cor )

f13 <- filter(boys_all, pop_yr == "fvw13")
f13 <- cbind(f12, prcomp(f12[,111:206])$x[,1:57])

png("../Figures/dsneur_wildwingsproj_f13_PC1to3.png")
pairs( f13[,207:212], lower.panel=panel.cor )
dev.off()

f14 <- filter(boys_all, pop_yr == "fvw14")
f14 <- cbind(f14, prcomp(f14[,111:206])$x[,1:57])

png("../Figures/dsneur_wildwingsproj_f14_PC1to3.png")
pairs( f14[,207:212], lower.panel=panel.cor )
dev.off()


#plotting 
#pairs( boysonly[, c(111, 113:115)], lower.panel=panel.cor )
# colors <- c("red", "black", "grey", "blue")[unclass(boysonly$pop_yr)]
#pairs( boysonly[, 111:114 ], lower.panel=panel.cor, bg = boysonly$pop_yr)

#pairs( boysonly[, 111:114 ], lower.panel=panel.cor, col=c("red", "black", "grey", "blue") )

#install.packages("GGally")

# png("../Figures/dsneur_wildwingsproj_PCto5.png")
# pairs( boysonly[, c(111, 113:118)], lower.panel=panel.cor )
# dev.off()

#add sup figure with F included. 

 
 
#ggplot(boysonly, aes(x = PC1, y = PC2, col = pop_yr)) + geom_point()
#ggplot(boysonly, aes(x = PC1, y = PC3, col = pop_yr)) + geom_point()
#ggplot(boysonly, aes(x = PC2, y = PC3, col = pop_yr)) + geom_point()


#PCs<- prcomp( boyssize )
#summary(PCs)

#PC2<- prcomp(wildwings[,5:100])
#summary(PC2)

#all <- cbind(wildwings, PC2$x[,1:56])
#ggplot(all, aes(x = CSize, y = PC1, col = pop_yr)) + geom_point()
#ggplot(all, aes(x = PC2, y = PC3, col = pop_yr)) + geom_point()


#Now a quick look in geomorph 
#reading in both sets of data (with and without F)

library(geomorph)
cord <- as.matrix(wildwings[,5:100])
shape <- arrayspecs(cord, 48, 2)
gdf <- geomorph.data.frame(shape = shape,
                           CS = wildwings$CSize,
                           pop = wildwings$pop,
                           pop_year = wildwings$pop_yr,
                           sex = wildwings$sex,
                           ind = wildwings$ind,
                           block = wildwings$block)



cord <- as.matrix(boysonly[,5:100])
shape <- arrayspecs(cord, 48, 2)
boys <- geomorph.data.frame(shape = shape,
                            CS = boysonly$CSize,
                            pop = boysonly$pop,
                            pop_year = boysonly$pop_yr,
                            ind = boysonly$ind,
                            block = boysonly$block)

source("~/Dropbox/KatiePelletier/KP_geomorphwingfunctions.R")

png("../Figures/all_wild_boysPlot.png")
plotAllSpecimens(boys$shape, 
                 links = wing.links)
dev.off()

boysmod <- procD.lm(shape ~ CS + pop_year, data = boys)

#uniqueallomod <- procD.lm(shape ~ CS * pop_year, data = boys)

#the CS:pop doesn't seem to add much. Can really use either model. 
#anova(uniqueallomod)

#This Rsq value is the same as when females are included (there population also explains about 15% of variance)
#Redo with a type 2. 
anova(boysmod)


boys.pair <- pairwise(boysmod, groups = boys$pop_year)


#Looks like PHO is the diffrent one?  
summary(boys.pair, test.type = "dist")



pairtable <- summary(boys.pair, test.type = "dist")[["summary.table"]]


write.csv(pairtable, file = "../Tables/WildPopulation_pairwiseTest.csv", quote = FALSE, row.names = TRUE)


#I want to plot the mean shapes of the pools together. 
ds_pools <- read.csv("../Data/selectedshape_75tails_ds.csv")

cord <- as.matrix(ds_pools[,10:105])
shape <- arrayspecs(cord, 48, 2)
ds_pools <- geomorph.data.frame(shape = shape,
                            CS = ds_pools$CSize,
                            pop = ds_pools$pop,
                            pop_year = ds_pools$pop_yr,
                            ind = ds_pools$ind,
                            block = ds_pools$block, 
                            tail = ds_pools$Tail)


group <- factor(paste(ds_pools$pop_year, ds_pools$tail))
levels(group)

new.coords <- coords.subset(ds_pools$shape, group = group)
names(new.coords) # see the list levels
# group shape means

pool_mean <- lapply(new.coords, mshape)
wing_flip <- matrix(rep(c(1, -1), 75), nrow = 75, ncol = 2, byrow = TRUE)


png("../Figures/ds_wildpools_shapechange_cmo_mag2.png")
plotRefToTarget(pool_mean[["cmo14 Left"]], 
                pool_mean[["cmo14 Right"]], 
                links = wing.links, method = "points", mag = 2, 
                gridPars=wing.spec )
dev.off()

png("../Figures/ds_wildpools_shapechange_pho_mag3.png")
plotRefToTarget(pool_mean[["pho14 Left"]], 
                pool_mean[["pho14 Right"]], 
                links = wing.links, method = "points", mag = 3, 
                gridPars=wing.spec)
dev.off()

png("../Figures/ds_wildpools_shapechange_fvw14_mag3.png")
plotRefToTarget(pool_mean[["fvw14 Left"]], 
                pool_mean[["fvw14 Right"]], 
                links = wing.links, method = "points", mag = 3, 
                gridPars=wing.spec)
dev.off()


png("../Figures/ds_wildpools_shapechange_fvw12_mag3.png")
plotRefToTarget(pool_mean[["fvw12 Left"]], 
                pool_mean[["fvw12 Right"]], 
                links = wing.links, method = "points", mag = 2, 
                gridPars=wing.spec)
dev.off()


#and the neur pools 
neur_dat <- read.csv("../Data/selectedshape_75tails_neur.csv")

cord <- as.matrix(neur_dat[,9:104])
shape <- arrayspecs(cord, 48, 2)
neur_pools <- geomorph.data.frame(shape = shape,
                                CS = neur_dat$CSize,
                                pop = neur_dat$pop,
                                pop_year = neur_dat$pop_yr,
                                ind = neur_dat$ind,
                                block = neur_dat$block, 
                                tail = neur_dat$tail)


group <- factor(paste(neur_pools$pop_year, neur_pools$tail))
levels(group)

new.coords <- coords.subset(neur_pools$shape, group = group)
names(new.coords) # see the list levels
# group shape means

pool_mean <- lapply(new.coords, mshape)



png("../Figures/wildpools_shapechange_cmo_mag2.png")
plotRefToTarget(pool_mean[["cmo14 L"]], 
                pool_mean[["cmo14 R"]], 
                links = wing.links, method = "points", mag = 3, 
                gridPars=wing.spec )
dev.off()

png("../Figures/_wildpools_shapechange_pho_mag2.png")
plotRefToTarget(pool_mean[["pho14 L"]], 
                pool_mean[["pho14 R"]], 
                links = wing.links, method = "points", mag = 2, 
                gridPars=wing.spec)
dev.off()

png("../Figures/_wildpools_shapechange_fvw14_mag2.png")
plotRefToTarget(pool_mean[["fvw14 L"]], 
                pool_mean[["fvw14 R"]], 
                links = wing.links, method = "points", mag = 2, 
                gridPars=wing.spec)
dev.off()


png("../Figures/_wildpools_shapechange_fvw12_mag2.png")
plotRefToTarget(pool_mean[["fvw12 L"]], 
                pool_mean[["fvw12 R"]], 
                links = wing.links, method = "points", mag = 1, 
                gridPars=wing.spec)
dev.off()



