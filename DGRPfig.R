#Drawing the cor plots between selection vecs and DGRP shape to get the numbers and make the figure for the paper. 
#using the landmarks because the estimated G didn't have this relationship 
source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R" )
source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R" )

# For Ian only
#source( "~/Dropbox/Dworkin_lab/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R" )
#source( "~/Dropbox/Dworkin_lab/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R" )

#library(tidyverse)
library(lme4)
#from the 2019 paper.
dgrp <- read.csv("../Data/BothLabs_Wings_28Oct.csv",
                 header=TRUE, nrows=22923, stringsAsFactors = TRUE)

projFunction <- function(x, y) {
  scalarProj <- (x %*% y) / norm(y, type = "2")
  return(scalarProj)
}

##### Edited from the WRP_functions file to round better.
panel.cor <- function(x, y) {
  r <- cor(x, y)
  par( usr =c(0, 1, 0, 1))
  Cor <- formatC(c(r, 0.123456789), format = "f", digits=2)[1]
  text(0.5, 0.5, paste(Cor), cex=1.5)
}
######

dgrp <- read.csv("../Data/BothLabs_Wings_28Oct.csv")

str(dgrp)
#line 761 has some 'probably M' this is one of the selected lines for the art sel experiment.
with(dgrp, table(Line, Sex))

hist(dgrp$Csize)
#Looks like a wing
WingPlot(PrcCrds=colMeans(dgrp[1,10:105]))

dwo_wings <- dgrp[dgrp$Lab == "Dwo", ]
# dwo_parent <- droplevels(dwo_wings[dwo_wings$Line %in% parents$Line,  ])
# head(dwo_wings)
# head(dwo_parent)
# with(dwo_parent, table(Line, Sex))
# nrow(with(dwo_parent, table(Line, Sex))) #only 19/30

#This is just total shape, will check common allometry just in case
# This ignores sex totally.
lm_shape <- lm(as.matrix(dwo_wings[,10:105]) ~ 0 + Line,
               data = dwo_wings)

# With sex

estimates_shape <- coef(lm_shape)
estimates_shape[1:5, 1:5]

# looks like a bunch of wings..
WingPlot(estimates_shape[1,])
WingPlot(estimates_shape[2,], add = T)
WingPlot(estimates_shape[3,], add = T, wingcol = "red")

WingPlot(colMeans(estimates_shape))

dim(estimates_shape)
rownames(estimates_shape)




# With sex effects accounted for as an average sex effect on shape
lm_shape_sex <- lm(as.matrix(dwo_wings[,10:105]) ~ 0 + Line:Sex,
               data = dwo_wings,
               subset = dwo_wings$Sex != "probablyM")

crap <- coef(lm_shape_sex)
dim(crap)

#first DGRP
WingPlot(crap[1,], wingcol = "red") # F
WingPlot(crap[(1+153),], wingcol = "red", winglty = 2, add = T) # M

#second DGRP
WingPlot(crap[2,], wingcol = "blue", add = T) # F
WingPlot(crap[(2+153),], wingcol = "blue", winglty = 2, add = T) # M

WingPlot(crap[3,], wingcol = "purple", add = T) # F
WingPlot(crap[(3+153),], wingcol = "purple", winglty = 2, add = T) # M

# Create a sex averaged estimate (new_crap has the sex averaged estimates)

new_crap <- matrix(NA, nrow = 153, ncol = 96)

for (i in 1:nrow(new_crap)) {
  temp_dat <- colMeans(crap[c(i,i+153),])
  new_crap[i,] <- temp_dat
}


WingPlot(new_crap [1,], wingcol = "blue")
WingPlot(new_crap [2,],  wingcol = "purple", add = T)


#Now need to grab the PCs for projection.
lines <- rownames(estimates_shape)

# please note KP... that this way is using the line means for the PCA, so it is kind of a quick and dirty G matrix.
# Will Pitchers, may have used the PCA for the full set of observations.
# So one potential source of differences in the projections as relating to the PCs is whether you are (like below) using the total phenotypic VCV matrix or the approximation of the genetic VCV matrix.

line.PC <- data.frame(lines, estimates_shape, prcomp(estimates_shape)$x[,1:56])
#just looks cleaner.
rownames(line.PC) <- c()
head(line.PC)


########################Shape change vecs######################

selvec <- read.csv("../Data/seldict_vectors.csv")
str(selvec)
vec2 <- read.csv("../Data/newdict_vectors_uncorrelated_to_ds.csv")
str(vec2)
neur <- vec2[5,2:97]


#now we want to get the relationship between selection vectors for the paper.

#ds and emc
#0.6553712
cor(as.numeric(selvec[1,3:98]), as.numeric(selvec[2,3:98]))

#ds and neur
#0.03464236
cor(as.numeric(selvec[1,3:98]), as.numeric(neur))

#emc and neur
#0.3008779
cor(as.numeric(selvec[2,3:98]), as.numeric(neur))

mshape <- colMeans(estimates_shape)

#This is all messed up?
#This doesn't look like the same as what I have drawn before. I think it is the DGRP data set.
#Also looks fucked witht the dgrp means. or with a single line mean.
# Maybe not common superimposition?
WingEffect( meanshape=as.matrix(mshape) , effectplus=as.matrix( selvec[1,3:98]/100 ),
            effectminus=as.matrix( selvec[1,3:98]/100), winglabel=paste(selvec[1,1]),
            scale.factor = 0.3, wingcol=c("black", "blue", "red") )


WingEffect( meanshape = estimates_shape[1,] , effectplus=as.matrix( selvec[1,3:98]/100 ),
            effectminus=as.matrix( selvec[1,3:98]/100), winglabel=paste(selvec[1,1]),
            scale.factor = 0.25, wingcol=c("black", "blue", "red") , meanline = T)

WingEffect( meanshape = estimates_shape[1,] , effectplus=as.matrix( selvec[2,3:98]/100 ),
            effectminus=as.matrix( selvec[2,3:98]/100), winglabel=paste(selvec[2,1]),
            scale.factor=5, wingcol=c("black", "blue", "red") , meanline = T)

diff_crap <- estimates_shape[1,] - estimates_shape[2,]


WingEffect( meanshape=estimates_shape[6,] , effectplus=as.matrix(diff_crap ),
            effectminus=as.matrix(diff_crap),
            scale.factor=2, wingcol=c("black", "blue", "red") , meanline = T)

#Now projecting the data onto the shape change vecs.
#using model estimated means.
line.PC$ds <- as.matrix(line.PC[,2:97]) %*% t(selvec[1,3:98])
line.PC$ds1 <- projFunction(x = as.matrix(line.PC[,2:97]), y = t(as.matrix(selvec[1,3:98])))


line.PC$emc <- as.matrix(line.PC[,2:97]) %*% t(selvec[2,3:98])
line.PC$neur <- as.matrix(line.PC[,2:97]) %*% t(neur)


#png(file = "../Figures/dgrp_all_corplot_DworkinWings.png")
pairs(line.PC[, c(154:156, 98:100)], lower.panel = panel.cor)
#dev.off()


##############Going to try with both labs#############

#fitted a more complete model here and using both labs.
str(dgrp)
with(dgrp, table(Line, Sex))
#going to take out the probably M line
dgrp <- droplevels(dgrp[dgrp$Sex != "probablyM",])
with(dgrp, table(Lab, Sex))
#also want to remove lines with less than 10 observations/sex.
#this is only actually two lines. line_306 &  line_256

dgrp <- dgrp[dgrp$Line != "line_306",]
dgrp <- dgrp[dgrp$Line != "line_256",]
dgrp <- droplevels(dgrp)
with(dgrp, table(Line, Sex))
dim(with(dgrp, table(Line, Sex)))

lines <- rownames(with(dgrp, table(Line, Sex)))

lm_dgrp <- lm(as.matrix(dgrp[, 10:105]) ~ 0 + Line:Sex,
              data = dgrp)

crap <- coef(lm_dgrp)
dim(crap)
rownames(crap)

#looks fine
#first DGRP
WingPlot(crap[1,], wingcol = "red") # F
WingPlot(crap[(1+170),], wingcol = "red", winglty = 2, add = T) # M

#second DGRP
WingPlot(crap[2,], wingcol = "blue", add = T) # F
WingPlot(crap[(2+170),], wingcol = "blue", winglty = 2, add = T) # M

WingPlot(crap[3,], wingcol = "purple", add = T) # F
WingPlot(crap[(3+170),], wingcol = "purple", winglty = 2, add = T) # M

# Create a sex averaged estimate (new_crap has the sex averaged estimates)

new_crap <- matrix(NA, nrow = 170, ncol = 96)

for (i in 1:nrow(new_crap)) {
  temp_dat <- colMeans(crap[c(i,i+170),])
  new_crap[i,] <- temp_dat
}


WingPlot(new_crap [1,], wingcol = "blue")
WingPlot(new_crap [2,],  wingcol = "purple", add = T)


#Quick and dirty G matrix
line.PC <- data.frame(lines, new_crap, prcomp(new_crap)$x[,1:56])
#just looks cleaner.
rownames(line.PC) <- c()
head(line.PC)


#Now projecting the data onto the shape change vecs.
#using model estimated means.
line.PC$ds <- projFunction(x = as.matrix(line.PC[,2:97]),
                           y = t(as.matrix(selvec[1,3:98])))

line.PC$emc <- projFunction(x = as.matrix(line.PC[,2:97]),
                            y = t(as.matrix(selvec[2,3:98])))

line.PC$neur <- projFunction(x = as.matrix(line.PC[,2:97]),
                             y = t(as.matrix(neur)))


#png(file = "../Figures/dgrp_all_corplot_bothlabs.png")
pairs(line.PC[, c(154:156, 98:100)], lower.panel = panel.cor)
#dev.off()

###########one more time but accounting for allometry############
#mean centering CS
dgrp$Csize_c <- dgrp$Csize - mean(dgrp$Csize)


#common allometrey is a bad assumption .
lm_dgrp <- lm(as.matrix(dgrp[, 10:105]) ~ 0 + Csize + Line:Sex,
              data = dgrp)

allo.coef <- coef(lm_dgrp)
allo.coef
dim(allo.coef)
rownames(allo.coef)

#looks fine
#first line is now CS
#first DGRP
#Do wings just look weird with the allmoetric shape taken out?
WingPlot(allo.coef[2,], wingcol = "red") # F
WingPlot(allo.coef[(2+170),], wingcol = "red", winglty = 2, add = T) # M

#second DGRP
WingPlot(allo.coef[3,], wingcol = "blue", add = T) # F
WingPlot(allo.coef[(3+170),], wingcol = "blue", winglty = 2, add = T) # M

# Create a sex averaged estimate (new_crap has the sex averaged estimates)
allo.coef2 <- matrix(NA, nrow = 170, ncol = 96)

for (i in 1:nrow(allo.coef2)) {
  temp_dat <- colMeans(allo.coef[c(i + 1,i+171),])
  allo.coef2[i,] <- temp_dat
}


WingPlot(allo.coef2[1,], wingcol = "blue")
WingPlot(allo.coef2[2,],  wingcol = "purple", add = T)


#Quick and dirty G matrix
allo.PC <- data.frame(lines, allo.coef2, prcomp(allo.coef2)$x[,1:56])
#just looks cleaner.
rownames(allo.PC) <- c()
head(allo.PC)

#PCA from landmarks
lm.pc <- prcomp(dgrp[,10:105])$x[,1:56]
lm.pc.allo <- data.frame()

#Now projecting the data onto the shape change vecs.
#using model estimated means.
allo.PC$ds <- projFunction(x = as.matrix(allo.PC[,2:97]),
                           y = t(as.matrix(selvec[1,3:98])))

allo.PC$emc <- projFunction(x = as.matrix(allo.PC[,2:97]),
                            y = t(as.matrix(selvec[2,3:98])))

allo.PC$neur <- projFunction(x = as.matrix(allo.PC[,2:97]),
                             y = t(as.matrix(neur)))


#png(file = "../Figures/dgrp_all_corplot_bothlabs_noAllo.png")
pairs(allo.PC[, c(154:156, 98:100)], lower.panel = panel.cor)
#dev.off()

###################Using all the data for projection?##############

#First want to model out sex and allometry
#This isn't working.
#I have no idea what I am doing.This is not right. I have no idea what Ian was talking about
lm_landmarks <- lm(as.matrix(dgrp[,10:105]) ~  Csize*Sex ,
                   data = dgrp)

summary(lm_landmarks)

mod.resid <- lm_landmarks$residuals

#this works.
dim(dgrp)
dim(mod.resid)

#PCA from landmarks
lm.pc <- prcomp(mod.resid)$x[,1:56]
lm.pc.allo <- data.frame(dgrp, mod.resid, lm.pc)

#now projecting onto the landmarks.
lm.pc.allo$ds <- projFunction(x = as.matrix(mod.resid),
                           y = t(as.matrix(selvec[1,3:98])))

lm.pc.allo$emc <- projFunction(x = as.matrix(mod.resid),
                            y = t(as.matrix(selvec[2,3:98])))

lm.pc.allo$neur <- projFunction(x = as.matrix(mod.resid),
                             y = t(as.matrix(neur)))

#Now I want to get the mean for each line.
lm.project.means <- aggregate(lm.pc.allo[,204:262],
                       by=list( Line=lm.pc.allo$Line),
                       FUN=mean )

dim(lm.project.means)


#png(file = "../Figures/dgrp_all_corplot_bothlabs_modelFromlm.png")
pairs(lm.project.means[, c(58:60, 2:4)], lower.panel = panel.cor)
#dev.off()


#Wondering if it is something in the data creating the cor or if it is the vectors. I think it is the data?
ds <- unlist(selvec[1,3:98])
emc <- unlist(selvec[2,3:98])
neur2 <- unlist(neur)

cor(ds, emc)
cor(ds, neur2)
cor(emc, neur2)


#########I want to deal with the funky looking wing effect plot by using geomorph rather than messing with this.

cord <- as.matrix(dwo_wings[,10:105])
shape <- arrayspecs(cord, 48, 2)
dgrp <- geomorph.data.frame(shape = shape,
                            CS = dwo_wings$Csize,
                            line = dwo_wings$Line,
                            rep = dwo_wings$Rep,
                            sex = dwo_wings$Sex)

#This has Maria's stuff added.
source("~/Dropbox/KatiePelletier/KP_geomorphwingfunctions.R")

#need a ref.
m <- mshape(dgrp$shape)
plot(m, links = wing.links)

#there is no clear way to plot the vector change.
#for now I am going to try to just calculate this with subtraction?

#first need put the ds shape change vec into the same format as this one
geo_ds <- matrix(ds, nrow = 48, ncol = 2, byrow = TRUE)

ds_change <- matrix(NA, nrow = 48, ncol = 2)
ds_change <- m - geo_ds/100

png("../Figures/ds_KDeffectPlot_geomorph_1x.png")
plotRefToTarget(m, ds_change, method = "points",
                links = wing.links, mag = 1,
                gridPars=wing.spec)
dev.off()

geo_emc <- matrix(emc, nrow = 48, ncol = 2, byrow = TRUE)

emc_change <- matrix(NA, nrow = 48, ncol = 2)
emc_change <- m - geo_emc/100

png("../Figures/emc_KDeffectPlot_geomorph_2x.png")
plotRefToTarget(m, emc_change, method = "points",
                links = wing.links, mag = 2,
                gridPars=wing.spec)
dev.off()


geo_neur <- matrix(as.numeric(neur), nrow = 48, ncol = 2, byrow = TRUE)

neur_change <- matrix(NA, nrow = 48, ncol = 2)
neur_change <- m - geo_neur/100

png("../Figures/neur_KDeffectPlot_geomorph_2x.png")
plotRefToTarget(m, neur_change, method = "points",
                links = wing.links, mag = 2,
                gridPars=wing.spec)
dev.off()


