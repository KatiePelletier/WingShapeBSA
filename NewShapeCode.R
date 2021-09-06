####Wild wings shape stuff for the paper and making final figures. 
library(tidyverse)

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
norm(ds_vec) #2.058989
norm(emc_vec) #0.2057204
norm(neur) #21.12176 why so giant? 


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
#Will included the F in the model so data was a little diffrent (did not subset boys)
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

# boysonly$PC1 <- prcomp(boyssize )$x[,1]
# boysonly$PC2 <- prcomp( boyssize )$x[,2]
# boysonly$PC3 <- prcomp( boyssize )$x[,3]
# boysonly$PC4 <- prcomp( boyssize )$x[,4]
# boysonly$PC5 <- prcomp( boyssize )$x[,5]

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

f12 <- filter(boys_all, pop_yr == "fvw12")
f12 <- cbind(f12, prcomp(f12[,111:206])$x[,1:57])

png("../Figures/dsneur_wildwingsproj_f12_PC1to3.png")
pairs( f12[,207:212], lower.panel=panel.cor )
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

plotAllSpecimens(boys$shape, 
                 links = wing.links)


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


write.csv(pairtable, file = "../Tables/WildPopulation_pairwiseTest.csv", quote = FALSE, row.names = FALSE)


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


png("../Figures/ds_wildpools_shapechange_cmo.png")
plotRefToTarget(pool_mean[["cmo14 Left"]], 
                pool_mean[["cmo14 Right"]], 
                links = wing.links, method = "points", mag = 1, 
                gridPars=wing.spec )
dev.off()

png("../Figures/ds_wildpools_shapechange_pho.png")
plotRefToTarget(pool_mean[["pho14 Left"]], 
                pool_mean[["pho14 Right"]], 
                links = wing.links, method = "points", mag = 1, 
                gridPars=wing.spec)
dev.off()

png("../Figures/ds_wildpools_shapechange_fvw14.png")
plotRefToTarget(pool_mean[["fvw14 Left"]], 
                pool_mean[["fvw14 Right"]], 
                links = wing.links, method = "points", mag = 1, 
                gridPars=wing.spec)
dev.off()


png("../Figures/ds_wildpools_shapechange_fvw12.png")
plotRefToTarget(pool_mean[["fvw12 Left"]], 
                pool_mean[["fvw12 Right"]], 
                links = wing.links, method = "points", mag = 1, 
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



png("../Figures/_wildpools_shapechange_cmo.png")
plotRefToTarget(pool_mean[["cmo14 L"]], 
                pool_mean[["cmo14 R"]], 
                links = wing.links, method = "points", mag = 1, 
                gridPars=wing.spec )
dev.off()

png("../Figures/_wildpools_shapechange_pho.png")
plotRefToTarget(pool_mean[["pho14 L"]], 
                pool_mean[["pho14 R"]], 
                links = wing.links, method = "points", mag = 1, 
                gridPars=wing.spec)
dev.off()

png("../Figures/_wildpools_shapechange_fvw14.png")
plotRefToTarget(pool_mean[["fvw14 L"]], 
                pool_mean[["fvw14 R"]], 
                links = wing.links, method = "points", mag = 1, 
                gridPars=wing.spec)
dev.off()


png("../Figures/_wildpools_shapechange_fvw12.png")
plotRefToTarget(pool_mean[["fvw12 L"]], 
                pool_mean[["fvw12 R"]], 
                links = wing.links, method = "points", mag = 1, 
                gridPars=wing.spec)
dev.off()



