poolwings <- read.csv("../Data/selectedshape_75tails_ds.csv", stringsAsFactors = T)

str(poolwings)
source('~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R', chdir = TRUE)

selvec <- read.csv( "../Data/seldict_vectors.csv" )
str( selvec )

WingPlot(as.matrix(colMeans(poolwings[,10:105])))

#Houle wings are way bigger (specifically, 100x bigger) than ours. This plot looks funny. Mean center data?
WingEffect( meanshape=as.matrix( colMeans( poolwings[,10:105] ) ), effectplus=as.matrix( selvec[1,3:98]/100 ),
            effectminus=as.matrix( selvec[1,3:98]/100), winglabel=paste(selvec[1,1]),
            scale.factor=1, wingcol=c("black", "blue", "red") )

ds_vec <- as.matrix(selvec[1,3:98])

vec2 <- read.csv("../Data/newdict_vectors_uncorrelated_to_ds.csv")
neur_vec <- vec2[5,2:97]

norm_vec <- function(x) sqrt(sum(x^2))

#need to get mean shape of each pool 

cmo14.left <- colMeans(filter(poolwings %>%
                       filter(pop_yr == "cmo14") %>%
                       filter(Tail == "Left") %>%
                       dplyr::select(x1:y48)
                     ))

cmo14.right <- colMeans(filter(poolwings %>%
                                filter(pop_yr == "cmo14") %>%
                                filter(Tail == "Right") %>%
                                dplyr::select(x1:y48)
))

cmo14.PD <- norm_vec(as.matrix(cmo14.left - cmo14.right))
#0.03328883
cmo14.PD 

pho14.left <- colMeans(filter(poolwings %>%
                                filter(pop_yr == "pho14") %>%
                                filter(Tail == "Left") %>%
                                dplyr::select(x1:y48)
))

pho14.right <- colMeans(filter(poolwings %>%
                                 filter(pop_yr == "pho14") %>%
                                 filter(Tail == "Right") %>%
                                 dplyr::select(x1:y48)
))

pho14.PD <- norm_vec(as.matrix(pho14.left - pho14.right))
#0.03618959
pho14.PD 

fvw14.left <- colMeans(filter(poolwings %>%
                                filter(pop_yr == "fvw14") %>%
                                filter(Tail == "Left") %>%
                                dplyr::select(x1:y48)
))

fvw14.right <- colMeans(filter(poolwings %>%
                                 filter(pop_yr == "fvw14") %>%
                                 filter(Tail == "Right") %>%
                                 dplyr::select(x1:y48)
))

fvw14.PD <- norm_vec(as.matrix(fvw14.left - fvw14.right))
#0.04126227
fvw14.PD 

fvw12.left <- colMeans(filter(poolwings %>%
                                filter(pop_yr == "fvw12") %>%
                                filter(Tail == "Left") %>%
                                dplyr::select(x1:y48)
))

fvw12.right <- colMeans(filter(poolwings %>%
                                 filter(pop_yr == "fvw12") %>%
                                 filter(Tail == "Right") %>%
                                 dplyr::select(x1:y48)
))

fvw12.PD <- norm_vec(as.matrix(fvw12.left - fvw12.right))
#0.03966965
fvw12.PD 

#I want an avg of all for of these 
#0.03760258
mean(c(cmo14.PD, pho14.PD, fvw12.PD, fvw14.PD))

cor((fvw12.left - fvw12.right), ds_vec[1,]) #0.9187666
cor((fvw14.left - fvw14.right), ds_vec[1,]) #0.9029401
cor((pho14.left - pho14.right), ds_vec[1,]) #0.7881412
cor((cmo14.left - cmo14.right), ds_vec[1,]) #0.9349219

#plotting means vs eachother
#Actually looks pretty legit?
WingPlot(cmo14.left, wingcol="black")
WingPlot(cmo14.right, wingcol="red", add=T)

WingPlot(pho14.left, wingcol="black")
WingPlot(pho14.right, wingcol="red", add=T)

WingPlot(fvw14.left, wingcol="black")
WingPlot(fvw14.right, wingcol="red", add=T)

WingPlot(fvw12.left, wingcol="black")
WingPlot(fvw12.right, wingcol="red", add=T)


#now for the neur pools 
neurwings <- read.csv("../Data/selectedshape_75tails_neur.csv")

neur.cmo14.left <- colMeans(filter(neurwings %>%
                                filter(pop_yr == "cmo14") %>%
                                filter(tail == "L") %>%
                                dplyr::select(x1:y48)
))

neur.cmo14.right <- colMeans(filter(neurwings %>%
                                 filter(pop_yr == "cmo14") %>%
                                 filter(tail == "R") %>%
                                 dplyr::select(x1:y48)
))

neur.cmo14.PD <- norm_vec(as.matrix(neur.cmo14.left - neur.cmo14.right))
#0.02748431
neur.cmo14.PD 

neur.pho14.left <- colMeans(filter(neurwings %>%
                                filter(pop_yr == "pho14") %>%
                                filter(tail == "L") %>%
                                dplyr::select(x1:y48)
))

neur.pho14.right <- colMeans(filter(neurwings %>%
                                 filter(pop_yr == "pho14") %>%
                                 filter(tail == "R") %>%
                                 dplyr::select(x1:y48)
))

neur.pho14.PD <- norm_vec(as.matrix(neur.pho14.left - neur.pho14.right))
#0.02863765
neur.pho14.PD 

neur.fvw14.left <- colMeans(filter(neurwings %>%
                                filter(pop_yr == "fvw14") %>%
                                filter(tail == "L") %>%
                                dplyr::select(x1:y48)
))

neur.fvw14.right <- colMeans(filter(neurwings %>%
                                 filter(pop_yr == "fvw14") %>%
                                 filter(tail == "R") %>%
                                 dplyr::select(x1:y48)
))

neur.fvw14.PD <- norm_vec(as.matrix(neur.fvw14.left - neur.fvw14.right))
#0.04126227
neur.fvw14.PD 

neur.fvw12.left <- colMeans(filter(neurwings %>%
                                filter(pop_yr == "fvw12") %>%
                                filter(tail == "L") %>%
                                dplyr::select(x1:y48)
))

neur.fvw12.right <- colMeans(filter(neurwings %>%
                                 filter(pop_yr == "fvw12") %>%
                                 filter(tail == "R") %>%
                                 dplyr::select(x1:y48)
))

neur.fvw12.PD <- norm_vec(as.matrix(neur.fvw12.left - neur.fvw12.right))
#0.03816062
neur.fvw12.PD



#I want an avg of all for of these 
#0.03232756
mean(c(neur.cmo14.PD, neur.pho14.PD, neur.fvw12.PD, neur.fvw14.PD))

#plotting means vs eachother
#Actually looks pretty legit?
WingPlot(neur.cmo14.left, wingcol="black")
WingPlot(neur.cmo14.right, wingcol="red", add=T)

WingPlot(neur.pho14.left, wingcol="black")
WingPlot(neur.pho14.right, wingcol="red", add=T)

WingPlot(neur.fvw14.left, wingcol="black")
WingPlot(neur.fvw14.right, wingcol="red", add=T)

WingPlot(neur.fvw12.left, wingcol="black")
WingPlot(neur.fvw12.right, wingcol="red", add=T)


cor((neur.fvw12.left - neur.fvw12.right), t(as.matrix(neur_vec[1,]))) #0.7924834
cor(as.matrix(neur.fvw14.left - neur.fvw14.right), t(neur_vec[1,])) #0.8014295
cor((neur.pho14.left - neur.pho14.right), t(neur_vec[1,])) #0.5859181
cor((neur.cmo14.left - neur.cmo14.right), t(neur_vec[1,]))#0.8329623

#do this 

# library(geomorph)
# cord <- as.matrix(poolwings[,10:105])
# shape <- arrayspecs(cord, 48, 2)
# gdf <- geomorph.data.frame(shape = shape,
#                            CS = poolwings$CSize,
#                            pop = poolwings$pop,
#                            pop_year = poolwings$pop_yr,
#                            sex = poolwings$sex,
#                            ind = poolwings$ind,
#                            tail = poolwings$Tail)
# 
# 
# 
# 
# #I want to model shape change and then use pairwise to look at diffrences between pools to get an estimate of PD.
# 
# mod <- procD.lm(shape ~ CS*pop_year*sex + tail, data = gdf)
# 
# #why no third order interactions? Not allowed in geomorph?
# summary(mod)
# 
# anova(mod)
# 
# 
# pair <- pairwise(mod, groups = poolwings$Tail)
# 
# #0.03627186 
# #About the same as what I saw before. 
# #There must be lots of variance in the data that makes the means look diffrent like this. 
# summary(pair, test.type = "dist")
