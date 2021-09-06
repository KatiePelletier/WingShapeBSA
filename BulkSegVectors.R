rm( list= ls() )

source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R" )
source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R" )

#setwd( "~/Dropbox/IanShare/BulkSegregant/" )
wings <- read.csv( "BSA_wings.csv" )

names( wings )

wings$ID

# the ID string needs to be parsed to extract categories, thus:
wings$year <- factor( ifelse( grepl( 'fvw_[a-z]', as.character( wings$ID ) ), "2014", "2013" ))
wings$sex <- factor( ifelse( grepl( '_f_', as.character( wings$ID ), fix=T ), "F", "M" ))
wings$side <- factor( ifelse( grepl( '_r', as.character( wings$ID ), fix=T ), "R", "L" ))
wings$ind <- as.integer( gsub( '[fvwlrm_]', '', as.character( wings$ID ) ) )

# now we'll need a scaled centroid size vector, with plot as sanity-check
wings$CSize <- wings$CSize * 0.007
plot( wings$CSize ~ wings$sex )

# shape sanity-check  <- looks good
WingPlot( PrcCrds=colMeans( wings[,2:98] ) )

###### INSERTED LATER  <- grabbing a single wing in order to make a nice figure for Annat (25-Mar-2014) #####
# using image 'FVW_0009_L.tif'
WingPlot( PrcCrds=wings[ wings$year=="2013" & wings$sex=="M" & wings$side=="L" & wings$ind==9 ,2:98],
          winglwd=4, wingcol="blue", wingpch=16, wingcex=1.5 )

# take mean shape within individual (i.e. collapse by side)
wings1 <- aggregate( wings[, 2:98], wings[c("year", "sex", "ind")], mean )
# sanity-check <- still good
WingPlot( PrcCrds=colMeans( wings1[,4:99] ) )

# write out for ease of Eladio
write.csv( wings1, "BSA_wings_by_ind.csv", row.names=F, quote=F )

###### SELECTION VECTORS FROM ELADIO ######
    # quoted from email 24/02/2014 #
# I'm attaching here the relevant selection vectors, in configuration space (96D).
# The first two rows are the selection gradients. If you want to visualize them, just
# add them to the reference configuration (third column), which is the global reference
# used to align all of the experiment's data. The rest of the rows contain the
# eigenvectors we used to project these data onto the same shape (58D) space as the
# selection gradients. This is all in shape, not shape-size space.

selvec <- read.csv( "seldict_vectors.csv" )
str( selvec )

# key information contained in:
selvec[1:3,]

# sanity-check on EM's mean shape
WingPlot( (as.matrix( selvec[3,3:98] ) / 100) )
# our wings are in reversed orientation & EM's LM coords are 100x the size!

# ...thus I shall build an 100x x-inverted (horizontal-flipped) version of our data
h_flip_mat <- matrix( c( -100, 100), ncol=96, nrow=nrow( wings1 ), byrow=T )
newwings <- data.frame( wings1[,1:3], ( wings1[,4:99] * h_flip_mat ), wings1$CSize )

# ...and re-do the sanity check
png( "mean_shape_sanity_check.png" )
  WingPlot( (as.matrix( selvec[3,3:98] ) / 100), winglwd=2 )
  WingPlot( PrcCrds=colMeans( newwings[,4:99]/100 ), add=T, wingcol="red", winglwd=2 )
  text( -0.1,0.2, "Houle", col="blue" ) ; text( 0.1,0.2, "Dworkin", col="red" )
dev.off()

# visualise vectors for ID & EM sanity-check-by-proxy

# ds
png( "ds_vec_EM.png" )
WingEffect( meanshape=as.matrix( selvec[3,3:98] )/100, effectplus=as.matrix( selvec[1,3:98] )/100,
            effectminus=as.matrix( selvec[1,3:98] )/100, winglabel=paste(selvec[1,1]),
            scale.factor=1, wingcol=c("black", "blue", "red") )
dev.off()

png( "ds_vec_WP.png" )
WingEffect( meanshape=as.matrix( colMeans( newwings[,4:99] )/100 ), effectplus=as.matrix( selvec[1,3:98] )/100,
            effectminus=as.matrix( selvec[1,3:98] )/100, winglabel=paste(selvec[1,1]),
            scale.factor=1, wingcol=c("black", "blue", "red") )
dev.off()

# emc
png( "emc_vec_EM.png" )
WingEffect( meanshape=as.matrix( selvec[3,3:98] )/100, effectplus=as.matrix( selvec[2,3:98] )/100,
            effectminus=as.matrix( selvec[2,3:98] )/100, winglabel=paste(selvec[2,1]),
            scale.factor=5, wingcol=c("black", "blue", "red") )
dev.off()

png( "emc_vec_WP.png" )
WingEffect( meanshape=as.matrix( colMeans( newwings[,4:99] )/100 ), effectplus=as.matrix( selvec[2,3:98] )/100,
            effectminus=as.matrix( selvec[2,3:98] )/100, winglabel=paste(selvec[2,1]),
            scale.factor=5, wingcol=c("black", "blue", "red") )
dev.off()


# now to project our new data onto EM's selected vectors

ds_vec <- as.matrix( selvec[1,3:98] )
emc_vec <- as.matrix( selvec[2,3:98] )

cor( t(selvec[1,3:98]), t(selvec[2,3:98]) )   # selection vectors themselves correlated at 0.66
ang.vec( t(selvec[1,3:98]), t(selvec[2,3:98]) )

newwings$ds <- as.matrix( newwings[,4:99] ) %*% t( ds_vec )
newwings$emc <- as.matrix( newwings[,4:99] ) %*% t( emc_vec )

# I can't resist having a look at the distribution...
png( "sel_vecs_correlated.png" )
  plot( newwings$ds, newwings$emc, col=densCols( newwings$ds, newwings$emc ), xlab="LM's * ds", ylab="LM's * emc",
        main=paste( "correlation=", round( cor( newwings$ds, newwings$emc ), 2) ), pch=16)
  abline( lm( newwings$emc ~ newwings$ds ), col="red" )
dev.off()

# ...this might be a problem: projecting our new data onto these vectors increases the correlation to 0.91

# let's see how much overlap we get...
num <- 50   # number of selected individuals

ds_BSA_up <- newwings[ order( newwings$ds, decreasing=T ), c(1:3, 101) ][1:num,]
ds_BSA_dn <- newwings[ order( newwings$ds, decreasing=F ), c(1:3, 101) ][1:num,]

emc_BSA_up <- newwings[ order( newwings$emc, decreasing=T ), c(1:3, 101) ][1:num,]
emc_BSA_dn <- newwings[ order( newwings$emc, decreasing=F ), c(1:3, 101) ][1:num,]

# not as bad as I'd feared <- we get 50-70% overlap between the selected individuals...
length( intersect( paste(ds_BSA_up$sex, ds_BSA_up$year, ds_BSA_up$ind, sep="_"),
                   paste(emc_BSA_up$sex, emc_BSA_up$year, emc_BSA_up$ind, sep="_") )) / num

length( intersect( paste(ds_BSA_dn$sex, ds_BSA_dn$year, ds_BSA_dn$ind, sep="_"),
                   paste(emc_BSA_dn$sex, emc_BSA_dn$year, emc_BSA_dn$ind, sep="_") )) / num


# NB:  ask EM for raw data to see how far art. sel. has gotten

##### 1 more check <- I'm going to regress out size, sex etc. to see if they drive the correlation

shapesize <- lm( as.matrix( newwings[,4:99] ) ~ newwings$wings1.CSize * newwings$sex * newwings$year )$resid

newwings2 <- data.frame( newwings, shapesize )

newwings2$ds2 <- as.matrix( shapesize ) %*% t( ds_vec )
newwings2$emc2 <- as.matrix( shapesize ) %*% t( emc_vec )

png( "sel_vecs_correlated1.png" )
  plot( newwings2$ds2, newwings2$emc2, col=densCols( newwings2$ds2, newwings2$emc2 ), xlab="LM's * ds", ylab="LM's * emc",
        main=paste( "correlation=", round( cor( newwings2$ds2, newwings2$emc2 )[1], 2) ), pch=16)
  abline( lm( newwings2$emc2 ~ newwings2$ds2 ), col="red" )
dev.off()



##### Ian-check for wonkiness #####
# compare the pattern above with ptc and ct vectors... (I used the )

ct_vectors <- read.csv( "ct_3panel_vectors.csv" )
ptc_vectors <- read.csv( "ptc_3_panel_coords.csv" )

# ct vs ptc VC between mutant effects <- 0.47
ang.vec.abs( unlist( ct_vectors[1,2:97] ), unlist( ptc_vectors[1,2:97] ) )

# ct vs ptc VC between DGRP SNPs <- 0.22
ang.vec.abs( unlist( ct_vectors[2,2:97] ), unlist( ptc_vectors[2,2:97] ) )

# ct vs ptc VC between MaNC2 SNPs <- 0.18
ang.vec.abs( unlist( ct_vectors[3,2:97] ), unlist( ptc_vectors[3,2:97] ) )

newwings2$ct <- as.matrix( shapesize ) %*% t( as.matrix( ct_vectors[1,2:97] )*100 )
newwings2$ptc <- as.matrix( shapesize ) %*% t( as.matrix( ptc_vectors[1,2:97] )*100 )

png( "ct_vs_ptc.png" )
  plot( newwings2$ct, newwings2$ptc, col=densCols( newwings2$ct, newwings2$ptc ), xlab="LM's * ct", ylab="LM's * ptc",
        main=paste( "correlation=", round( cor( newwings2$ct, newwings2$ptc )[1], 2) ), pch=16)
  abline( lm( newwings2$ptc ~ newwings2$ct ), col="red" )
dev.off()

##### second set of (uncorrelated) vectors from EM 8/3/14 #####

newvecs <- read.csv( "newdict_vectors_uncorrelated_to_ds.csv" )
str(newvecs)

# sanity check
WingEffect( as.matrix(selvec[3,3:98])/100, as.matrix( newvecs[3,2:97] )/100,
            as.matrix( newvecs[3,2:97] )/100, scale.factor=2 )

sd_effect <- read.csv( "sd_effect.csv" )
sd_vec <- t( sd_effect[,2] )

# sanity check
WingEffect( as.matrix(selvec[3,3:98])/100, sd_vec, sd_vec, scale.factor=2 )


# project onto new vectors
newwings2$vg <- as.matrix( shapesize ) %*% t( as.matrix( newvecs[1,2:97] ) )
newwings2$yki <- as.matrix( shapesize ) %*% t( as.matrix( newvecs[2,2:97] ) )
newwings2$hh <- as.matrix( shapesize ) %*% t( as.matrix( newvecs[3,2:97] ) )
newwings2$vn <- as.matrix( shapesize ) %*% t( as.matrix( newvecs[4,2:97] ) )
newwings2$neur <- as.matrix( shapesize ) %*% t( as.matrix( newvecs[5,2:97] ) )
newwings2$bx <- as.matrix( shapesize ) %*% t( as.matrix( newvecs[6,2:97] ) )
newwings2$sd <- as.matrix( shapesize ) %*% t( sd_vec )

newwings2$PC1 <- prcomp( newwings2[,4:99] )$x[,1]
newwings2$PC2 <- prcomp( newwings2[,4:99] )$x[,2]

# express the correlations between projected vectors for interpretation
# 1st: a pairs plot
png( "relationship_among_projected_vectors.png", width=1000, height=1000 )
  pairs( newwings2[, 199:211 ], lower.panel=panel.cor )
dev.off()

#Katie: for drawing: rename some of these. 

# 2nd: a matrix of correlations
cormat <- cor( newwings2[, 199:209] )
write.csv( round(cormat, 3), "projected_vector_cor_matrix.csv" )

# 3rd: a heatmap
png( "proj_vector_cor_heatmap.png" )
  heatmap( cormat, symm=T, Rowv=NA, Colv=NA, revC=T, col=heat.colors( n=10 ) )
dev.off()

#I like the cor matrix better than this heat map... (Katie)
