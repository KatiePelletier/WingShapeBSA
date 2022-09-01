library(data.table)
library(tidyverse)
library(lme4)
library(car)
library(effects)
library(emmeans)
library(glmmTMB)
library(broom.mixed)

projFunction <- function(x, y) {
  scalarProj <- (x %*% y) / norm(y, type = "2")
  return(scalarProj)
}

#making the figure to look at the response to selection (phenotype)

dat <- fread("../Data/slid_lmdata0_SDVds.dat")
head(dat)

#I want to also add the selection data for getting rel h^2 
#Why are there three more flies in this file. FUCK ITTTTTTT
sel.dat <-fread("../Data/ds_allDat_selection.dat")

#This dropps those three extras. 
dat.all <- left_join(dat, sel.dat)

#Just trying to figure out what I have
wings <- separate(dat.all, "Tags", into = c("exp", "line", "rep", "gen"))

#only ds data .
levels(as.factor(wings$exp))

#ug. case problems.
levels(as.factor(wings$line))
wings$line <- ifelse(wings$line == "dn", "DN", wings$line)
levels(as.factor(wings$line))


levels(as.factor(wings$rep))
#Do we not have shape data for gen 8? DID WE ONLY DO 7 generations?? Have I been lied to THIS WHOLE TIME????
levels(as.factor(wings$gen))
wings$gen <- as.numeric(gsub("Gen", "", wings$gen))


#Now I think I want to calculate a shape score for each wing?

#reading in the ds vec 
selvec <- read.csv( "../Data/seldict_vectors.csv" )
head( selvec )
ds_vec <- as.matrix(selvec[1,3:98])


#first, I am going to model out alometry. 

shaperesid <- lm(as.matrix(wings[,16:111]) ~ log(CS), data = wings)$residuals

wings$ds <- projFunction(x = as.matrix(shaperesid) , y = t(ds_vec ))

head(wings$ds)

wings2 <- data.frame(wings, shaperesid)


#now I want to plot the response to selection. 
#this looks like poop. 
# ggplot(wings, aes(x = gen, y = ds, col = interaction(wings$line, wings$rep))) + 
#   geom_smooth(method = "lm")
# 
# 
# ggplot(wings, aes(x = gen, y = ds, col = interaction(wings$line, wings$rep))) + 
#   geom_line()

#insted I will just plot the mean (this should probably be modeled?)
#this confounds sex in there too. 

means <- (wings %>%
            group_by(interaction(wings$line, wings$rep, wings$gen, wings$Sex))  %>%
            summarise(ds = mean(ds)) %>%
            separate("interaction(wings$line, wings$rep, wings$gen, wings$Sex)", into = c("line", "rep", "gen", "Sex")) )

means
means$gen <- as.numeric(means$gen)

male_means <- filter(means, Sex == "M")
female_means <- filter(means, Sex == "F")

male.line.avg <- (male_means %>%
                    unite(line.gen, line, gen) %>%
                    mutate(line.gen2 = as.factor(line.gen)) %>%
                    group_by(line.gen2) %>%
                    summarise(ds = mean(ds)) %>%
                    separate(line.gen2, into = c("line", "gen"))
                      )

female.line.avg <- (female_means %>%
                    unite(line.gen, line, gen) %>%
                    mutate(line.gen2 = as.factor(line.gen)) %>%
                    group_by(line.gen2) %>%
                    summarise(ds = mean(ds)) %>%
                    separate(line.gen2, into = c("line", "gen"))
)



#Need to add the shape stuff. 
#A problem here is that I don't have generation 0, or the starting population. In the plot I have they just mark generation 0 as having a score of 0. 
m.response <- ggplot(male_means, aes(x = gen, y = ds, col = rep, shape = line)) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  theme_classic() + 
  ylab("ds shape score") +
  xlab("Generation") + 
  theme(legend.position = "none") +
  scale_x_continuous(breaks = 1:7)

ggplot(female_means, aes(x = gen, y = ds, col = rep, shape = line)) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  theme_classic() + 
  ylab("ds shape score") +
  xlab("Generation") + 
  theme(legend.position = "none") +
  scale_x_continuous(breaks = 1:7)

# ggplot(means, aes(x = gen, y = ds, col = line, shape = rep)) + 
#   geom_point() + 
#   geom_line(alpha = 0.5) + 
#   theme_classic() + 
#   ylab("ds shape score") +
#   theme(legend.position = "none")



#I also want to plot this with real data points and model estimates to see if this looks like trash or not. 

wings$CS <- as.numeric(wings$CS)
wings$Sex <- as.factor(wings$Sex)
wings$line <- as.factor(wings$line)
wings$gen <- as.numeric(wings$gen)

wings$gen0 <- wings$gen - 1 
wings$rep <- as.factor(wings$rep)
wings$CS_0 <- wings$CS - mean(wings$CS)


# ds.sel.shape_mod <- lmer(ds ~ CS + Sex + line*gen0 + (1|line:rep), 
#                          data = wings)
# summary(ds.sel.shape_mod)
# 
# 
# 
# ds.sel.shape_mod3 <- lmer(ds ~ CS + Sex + line*gen0 + (1 + gen0||line:rep) , 
#                          data = wings)
# 
# summary(ds.sel.shape_mod3)



ds.sel.shape_mod2 <- lmer(ds ~ (CS + Sex + line + gen0)^3 + (1 + gen0||line:rep) , 
                          data = wings)
summary(ds.sel.shape_mod2)
Anova(ds.sel.shape_mod2)

# crap <- emmeans(ds.sel.shape_mod, ~ line|Sex)
# crap
# plot(crap)

crap2 <- emmeans(ds.sel.shape_mod2, ~gen0|line)
summary(crap2)

#plot(allEffects(ds.sel.shape_mod2))

allEffects(ds.sel.shape_mod2)
str(wings)

#sex does matter 
Anova(ds.sel.shape_mod2)

#line:gen effect is what I really care about 

#allEffects(ds.sel.shape_mod)

linegen_effect <- predictorEffect("gen0", by = c("line"), ds.sel.shape_mod)
linegen_effect2 <- data.frame(linegen_effect)

linegen_effect2$gen <- linegen_effect2$gen0 + 1

linegen_effect2$ds <- linegen_effect2$fit

#I think I actually want to keep lines seperate rather than smooth over them here. 
ggplot(female_means, aes(x = gen, y = ds, shape = line, col = rep)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.2) + 
  scale_colour_manual(values=c('black', 'grey46', 'grey57')) +
  geom_line(data = linegen_effect2, col = "red") + 
  theme_classic() + 
  ylab(expression(paste(italic("ds"), " shape score")))+
  xlab("Generation") + 
  theme(legend.position = "none")  + 
  scale_x_continuous(breaks = 1:7)  


#####an aside to look at size in response to selection##################
#is there a size change assosiated with selection 

ds.sel.size.mod <- lmer(CS ~ (Sex * line * gen0) - Sex:line + 
                          (1 + gen0|line:rep) , 
                        data = wings, subset = wings$line != "CR")

summary(ds.sel.size.mod)
#                   Estimate Std. Error t value
#lineUP:gen0      -0.0097745  0.0108937  -0.897

#No line:gen (what we really care about). 
#                   Chisq Df Pr(>Chisq)
#line:gen0         1.2319  1    0.26704  
Anova(ds.sel.size.mod)

plot(predictorEffect("gen0", by = "line", ds.sel.size.mod))

size_lineGen_effect <- predictorEffect("gen0", by = "line", ds.sel.size.mod)
size_lineGen_effect2 <- data.frame(size_lineGen_effect)

size_lineGen_effect2$gen <- size_lineGen_effect2$gen0 + 1

size_lineGen_effect2$size <- size_lineGen_effect2$fit

#I want to make a nice (eveything on one plot) figure for the supp for this. 
size_lineGen_effect2$selection <- size_lineGen_effect2$line

png("../Figures/ds_size_dsSelection_relationship.png", width =2000, height = 2000, units = "px",res = 300)
ggplot(size_lineGen_effect2, aes(x = gen0, y = size, linetype = Sex, col = selection)) + 
  geom_line() +  
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = line), alpha = 0.3, color = NA) + 
  theme_classic() +   
  ylab("Centroid Size") +
  xlab("Generation") 

dev.off()


#I will use Will's plotting functions here. 
source('~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R', chdir = TRUE)

source('~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R', chdir = TRUE)

####################################
#plots for the paper
#I want to do a wing blur comparing the variation between up and down wings. 

crapwingsf <- (wings %>% 
                filter(line == "DN" | line == "UP") %>%
                filter(Sex == "F") %>%
                 filter(gen0 == 6))

crapwingsf$line <- droplevels(crapwingsf$line)
wingflip <- rep(c(-1, 1), 48)


png("../Figures/dsSel_F_wingBlur.png")
WingBlur3(as.matrix(crapwingsf[,16:111]), grouping_var = crapwingsf$line, groups = T)
dev.off()

crapwingsm <- (wings %>% 
                 filter(line == "DN" | line == "UP") %>%
                 filter(Sex == "M") %>%
                 filter(gen0 == 6))

crapwingsm$line <- droplevels(crapwingsm$line)

png("../Figures/dsSel_M_wingBlur.png")
WingBlur3(as.matrix(crapwingsm[,16:111]), grouping_var = crapwingsf$line, groups = T)
dev.off()
  
#Old shape change stuff that I changed to do in geomorph anyways. 

# ups <- filter(wings, line == "UP", gen == 7)
# up_mean <- colMeans(ups[,16:111])
# 
# down <- filter(wings, line == "DN", gen == 7)
# down_mean <- colMeans(down[,16:111])
# 
# C <- filter(wings, line == "CR", gen == 7)
# C_mean <- colMeans(C[,16:111])
# 
# Gen1_all <- filter(wings, gen == 1)
# all_start_mean <- colMeans(Gen1_all[,16:111])
# 
# Gen1_control <- filter(wings, line == "CR", gen == 1)
# C_start_mean <- colMeans(Gen1_control[,16:111])
# 
# Gen1_down <- filter(wings, line == "DN", gen == 1)
# DN_start_mean <- colMeans(Gen1_down[,16:111])
# 
# Gen1_up <- filter(wings, line == "UP", gen == 1)
# UP_start_mean <- colMeans(Gen1_up[,16:111])
# 
# 
# #Now I want to calculate PD between gen 1 an d gen 2 as well as find the cor between the 1 -> 7 shape change vec and ds 
# 
# #flip this so the sign is +
# 
# PD(C_start_mean - C_mean) #0.005146645
# cor((C_mean - C_start_mean), ds_vec[1,]) #0.-260007
# 
# c_diff <- C_start_mean - C_mean
# 
# c_reflect <- rep(c(-1, 1), 48) * C_start_mean
# 
# png("../Figures/ds_control_wing_gen1v7selection.png",width = 1000, height = 1000, res = 300, units = "px", bg = "transparent")
# WingPlot(C_start_mean, wingcol="black")
# WingPlot(C_mean, wingcol="red", add=T)
# dev.off()
# 
# PD(UP_start_mean - up_mean) #0.04200999
# cor((UP_mean - up_start_mean), ds_vec[1,]) # -0.9006606
# 
# png("../Figures/ds_up_wing_gen1v7selection.png", width = 1000, height = 1000, units = "px", res = 300,bg = "transparent")
# WingPlot(UP_start_mean, wingcol="black")
# WingPlot(up_mean, wingcol="red", add=T)
# dev.off()
# 
# #or using wing effect 
# up_diff <- UP_start_mean - up_mean
# 
# up_reflect <- rep(c(-1, 1), 48) * UP_start_mean
# 
# 
# png("../Figures/ds_up_wing_selEffect.png",width =1000, height = 1000, units = "px", res = 300,bg = "transparent")
# 
# WingEffect(up_reflect, up_diff, up_diff,
#                    wingcol=c("black", "black", "red"),
#                    scale.factor = 0.5,
#                    scale.display = FALSE,
#                    wingframe = FALSE,
#                    winglwd=c(1, 1, 1))
# dev.off()
# 
# 
# PD(DN_start_mean - down_mean) # 0.02206913
# cor((DN_mean - down_start_mean), ds_vec[1,]) # 0.8080928
# 
# down_diff <- DN_start_mean - down_mean
# 
# down_reflect <- rep(c(-1, 1), 48) * DN_start_mean
# 
# png("../Figures/ds_down_wing_gen1v7selection.png",width =1000, height = 1000, units = "px", res = 300,bg = "transparent")
# WingPlot(DN_start_mean, wingcol="black")
# WingPlot(down_mean, wingcol="red", add=T)
# dev.off()

#Plotting these weird ass shape changes in geomorph because I am over dealing with solving the strechy problem. 

source("~/Dropbox/KatiePelletier/KP_geomorphwingfunctions.R")
library(geomorph)

cord <- as.matrix(wings[,16:111])
shape <- arrayspecs(cord, 48, 2)
wing.dat <- geomorph.data.frame(shape = shape,
                                CS = wings$CS,
                                line = wings$line,
                                rep = wings$rep,
                                gen = wings$gen,
                                sex = wings$Sex)


group <- factor(paste(wing.dat$gen, wing.dat$line, wing.dat$sex))
levels(group)

new.coords <- coords.subset(wing.dat$shape, group = group)
names(new.coords) # see the list levels
# group shape means

group_means <- lapply(new.coords, mshape)



png("../Figures/ds_selection_shapechange_UP_F1v7_2x.png")
plotRefToTarget(group_means[["1 UP F"]], 
                group_means[["7 UP F"]], 
                links = wing.links, method = "points", mag = 2, 
                gridPars=wing.spec )
dev.off()


png("../Figures/ds_selection_shapechange_UP_M1v7_2x.png")
plotRefToTarget(group_means[["1 UP M"]], 
                group_means[["7 UP M"]], 
                links = wing.links, method = "points", mag = 2, 
                gridPars=wing.spec )
dev.off()


png("../Figures/ds_selection_shapechange_DN_F1v7_2x.png")
plotRefToTarget(group_means[["1 DN F"]], 
                group_means[["7 DN F"]], 
                links = wing.links, method = "points", mag = 2, 
                gridPars=wing.spec )
dev.off()

png("../Figures/ds_selection_shapechange_DN_M1v7_2x.png")
plotRefToTarget(group_means[["1 DN M"]], 
                group_means[["7 DN M"]], 
                links = wing.links, method = "points", mag = 2, 
                gridPars=wing.spec )
dev.off()

png("../Figures/ds_selection_shapechange_CN_F1v7_2x.png")
plotRefToTarget(group_means[["1 CR F"]], 
                group_means[["7 CR F"]], 
                links = wing.links, method = "points", mag = 1, 
                gridPars=wing.spec )
dev.off()

png("../Figures/ds_selection_shapechange_CN_M1v7_2x.png")
plotRefToTarget(group_means[["1 CR M"]], 
                group_means[["7 CR M"]], 
                links = wing.links, method = "points", mag = 1, 
                gridPars=wing.spec )
dev.off()


###########vec cor and mean shape change#####################


wings.shape.means <- (wings %>%
                        group_by(interaction(line, gen0, Sex, rep)) %>%
                        summarise_at(names(.)[16:111],
                                     .funs = c(mean="mean")) %>%
                        separate("interaction(line, gen0, Sex, rep)", 
                                 into = c("line", "gen0", "Sex", "rep"))
)

wings.shape.means

wings.shape.means[1,]
#My code is way faster. But it is a tiny bit diffrent. why
colMeans(wings[wings$line == "CR" & wings$gen0 == 0 & wings$Sex == "F", wings$rep == "A", 16:111])


C.gen1.f <- colMeans(wings.shape.means[wings.shape.means$line == "CR" & 
                                         wings.shape.means$gen0 == 0 &
                                         wings.shape.means$Sex == "F", 5:100])
C.gen7.f <- colMeans(wings.shape.means[wings.shape.means$line == "CR" & 
                                         wings.shape.means$gen0 == 6 & 
                                         wings.shape.means$Sex == "F", 5:100])

deltaC.f <- C.gen7.f - C.gen1.f 

#-0.3001403
cor(deltaC.f, as.numeric(selvec[1,3:98]))

PD(deltaC.f) # 0.005956335

C.gen1.m <- colMeans(wings.shape.means[wings.shape.means$line == "CR" & 
                                         wings.shape.means$gen0 == 0 & 
                                         wings.shape.means$Sex == "M", 5:100])

C.gen7.m <- colMeans(wings.shape.means[wings.shape.means$line == "CR" & 
                                         wings.shape.means$gen0 == 6 & 
                                         wings.shape.means$Sex == "M", 5:100])

deltaC.m <- C.gen7.m - C.gen1.m
# -0.1728336
cor(deltaC.m, as.numeric(selvec[1,3:98]))
PD(deltaC.m) # 0.005141316


dn.gen1.f <- colMeans(wings.shape.means[wings.shape.means$line == "DN" & 
                                           wings.shape.means$gen0 == 0 & 
                                           wings.shape.means$Sex == "F", 5:100])

dn.gen7.f <-  colMeans(wings.shape.means[wings.shape.means$line == "DN" & 
                                           wings.shape.means$gen0 == 6 & 
                                           wings.shape.means$Sex == "F", 5:100])

delta.dn.f <- dn.gen7.f - dn.gen1.f 
# -0.8237332
cor(delta.dn.f, as.numeric(selvec[1,3:98]))

PD(delta.dn.f) #0.02189198

dn.gen1.m <- colMeans(wings.shape.means[wings.shape.means$line == "DN" & 
                                          wings.shape.means$gen0 == 0 & 
                                          wings.shape.means$Sex == "M", 5:100])
dn.gen7.m <- colMeans(wings.shape.means[wings.shape.means$line == "DN" & 
                                          wings.shape.means$gen0 == 6 & 
                                          wings.shape.means$Sex == "M", 5:100])

delta.dn.m <- dn.gen7.m - dn.gen1.m
# -0.7738024
cor(delta.dn.m, as.numeric(selvec[1,3:98]))

PD(delta.dn.m ) #0.02278963

up.gen1.f <- colMeans(wings.shape.means[wings.shape.means$line == "UP" & 
                                          wings.shape.means$gen0 == 0 & 
                                          wings.shape.means$Sex == "F", 5:100])
up.gen7.f <- colMeans(wings.shape.means[wings.shape.means$line == "UP" & 
                                          wings.shape.means$gen0 == 6 & 
                                          wings.shape.means$Sex == "F", 5:100])

delta.up.f <- up.gen7.f - up.gen1.f
# 0.8962785
cor(delta.up.f, as.numeric(selvec[1,3:98]))

PD(delta.up.f) #0.03947258


up.gen1.m <- colMeans(wings.shape.means[wings.shape.means$line == "UP" & 
                                          wings.shape.means$gen0 == 0 & 
                                          wings.shape.means$Sex == "M", 5:100])
up.gen7.m <- colMeans(wings.shape.means[wings.shape.means$line == "UP" & 
                                          wings.shape.means$gen0 == 6 & 
                                          wings.shape.means$Sex == "M", 5:100])

delta.up.m <- up.gen7.m - up.gen1.m
# 0.8949903
cor(delta.up.m, as.numeric(selvec[1,3:98]))

PD(delta.up.m) #0.04460296

##########Gettin the rel  h^2 from selection differentals##############

#Starting with the ups. 
str(wings)
#isSelected is the T/F for which indiv are breeders for the next generation. 

wings.up <- filter(wings, line == "UP")
str(wings.up)

#getting means of each generation (total)

wings.up.gen.mean <- (wings.up %>%
                        group_by(interaction(gen0, Sex, rep)) %>%
                        summarize(gen_mean = mean(ds)) %>%
                        separate("interaction(gen0, Sex, rep)", 
                                 into = c("gen0", "Sex", "rep"))
)

wings.up.gen.mean

wings.up.sel.mean <- (wings.up %>%
                        filter(isSelected == 1) %>%
                        group_by(interaction(gen0, Sex, rep)) %>%
                        summarize(sel_mean = mean(ds)) %>%
                        separate("interaction(gen0, Sex, rep)", 
                                 into = c("gen0", "Sex", "rep"))
)

wings.up.sel.mean

wings.up.means <- left_join(wings.up.gen.mean, wings.up.sel.mean)

#I don't know if I will need this later, but now I have this. 
wings.up.means$gen0.sex.rep <- paste(wings.up.means$gen0, wings.up.means$Sex, wings.up.means$rep, sep = ".")
wings.up.means

wings.up.means$S <- wings.up.means$sel_mean - wings.up.means$gen_mean
wings.up.means

#getting the diffrences between generations within line. 
#I can't find a good way to get this so I'm just going to do it by hand and figure out the pretty way later (or not)

#Going to fill a df so that I can plot and do the regression after. 

Rc.F.A.up <- cumsum(diff(wings.up.means[wings.up.means$Sex == "F" & 
                                          wings.up.means$rep == "A",]$gen_mean, 
                         lag = 1))

Sc.F.A.up  <- cumsum(wings.up.means[wings.up.means$Sex == "F" & wings.up.means$rep == "A",]$S)[1:6]

relh.df.up <- data.frame(Sc.F.A.up , Rc.F.A.up )

#filling in the rest. 
relh.df.up$Rc.F.B.up  <- cumsum(diff(wings.up.means[wings.up.means$Sex == "F" & 
                                                      wings.up.means$rep == "B",]$gen_mean, lag = 1))
relh.df.up$Rc.F.C.up  <- cumsum(diff(wings.up.means[wings.up.means$Sex == "F" & 
                                                      wings.up.means$rep == "C",]$gen_mean,lag = 1))
relh.df.up$Rc.M.A.up  <- cumsum(diff(wings.up.means[wings.up.means$Sex == "M" & 
                                                      wings.up.means$rep == "A",]$gen_mean,lag = 1))
relh.df.up$Rc.M.B.up  <- cumsum(diff(wings.up.means[wings.up.means$Sex == "M" & 
                                                      wings.up.means$rep == "B",]$gen_mean, lag = 1))
relh.df.up$Rc.M.C.up <- cumsum(diff(wings.up.means[wings.up.means$Sex == "M" & 
                                                     wings.up.means$rep == "C",]$gen_mean, lag = 1))


relh.df.up$Sc.F.B.up <-cumsum(wings.up.means[wings.up.means$Sex == "F" & wings.up.means$rep == "B",]$S)[1:6]
relh.df.up$Sc.F.C.up <-cumsum(wings.up.means[wings.up.means$Sex == "F" & wings.up.means$rep == "C",]$S)[1:6]
relh.df.up$Sc.M.A.up <-cumsum(wings.up.means[wings.up.means$Sex == "M" & wings.up.means$rep == "A",]$S)[1:6]
relh.df.up$Sc.M.B.up <-cumsum(wings.up.means[wings.up.means$Sex == "M" & wings.up.means$rep == "B",]$S)[1:6]
relh.df.up$Sc.M.C.up <-cumsum(wings.up.means[wings.up.means$Sex == "M" & wings.up.means$rep == "C",]$S)[1:6]

#adding gen 
relh.df.up$gen <- seq(1:6)

#Now to make this long. 

relh.long.up <- pivot_longer(relh.df.up, cols = c("Sc.F.A.up":"Sc.M.C.up"),names_to = "tag", values_to = "score")

relh.long.up <- separate(relh.long.up, tag, into = c("stat", "sex", "rep", "treat"))

#this is dumb but now I want to fix the two stats I think 
relh.final.up <- pivot_wider(relh.long.up, names_from = "stat", values_from = "score")

relh.final.up

#This seems ok. Everything is in the right direction. 
ggplot(relh.final.up, aes(x = Sc, y = Rc, col = rep, shape = sex )) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  theme_classic() + 
  ylab("cumulative resopnse") +
  xlab("cumulative selection diff")

#Now to do the regression. 

up.mod <- lmer(Rc ~ Sc*sex + (1 |rep) , data = relh.final.up)
#up.mod_v2 <- glmmTMB(Rc ~ Sc*sex + (1 |rep) , data = relh.final.up)

summary(up.mod)
#summary(up.mod_v2)

#The slope is the realized haritability
#Defult is 95% conf intrival

#fixed  Sc            0.377    CI:  0.254      0.501  
tidy(up.mod_v2, effects="fixed", conf.int = T)

#Now for the down. 
wings.dn <- filter(wings, line == "DN")
str(wings.dn)

#getting means of each generation (total)

wings.dn.gen.mean <- (wings.dn %>%
                        group_by(interaction(gen0, Sex, rep)) %>%
                        summarize(gen_mean = mean(ds)) %>%
                        separate("interaction(gen0, Sex, rep)", 
                                 into = c("gen0", "Sex", "rep"))
)

wings.dn.gen.mean

wings.dn.sel.mean <- (wings.dn %>%
                        filter(isSelected == 1) %>%
                        group_by(interaction(gen0, Sex, rep)) %>%
                        summarize(sel_mean = mean(ds)) %>%
                        separate("interaction(gen0, Sex, rep)", 
                                 into = c("gen0", "Sex", "rep"))
)

wings.dn.sel.mean


wings.dn.means <- left_join(wings.dn.gen.mean, wings.dn.sel.mean)

#I don't know if I will need this later, but now I have this. 
wings.dn.means$gen0.sex.rep <- paste(wings.dn.means$gen0, wings.dn.means$Sex, wings.dn.means$rep, sep = ".")
wings.dn.means

wings.dn.means$S <- wings.dn.means$sel_mean - wings.dn.means$gen_mean
wings.dn.means

#getting the diffrences between generations within line. 
#I can't find a good way to get this so I'm just going to do it by hand and figure out the pretty way later (or not)

#Going to fill a df so that I can plot and do the regression after. 

Rc.F.A.dn <- cumsum(diff(wings.dn.means[wings.dn.means$Sex == "F" & 
                                          wings.dn.means$rep == "A",]$gen_mean, lag = 1))

Sc.F.A.dn  <- cumsum(wings.dn.means[wings.dn.means$Sex == "F" & wings.dn.means$rep == "A",]$S)[1:6]

relh.df.dn <- data.frame(Sc.F.A.dn , Rc.F.A.dn )

#filling in the rest. 
relh.df.dn$Rc.F.B.dn  <- cumsum(diff(wings.dn.means[wings.dn.means$Sex == "F" & 
                                                      wings.dn.means$rep == "B",]$gen_mean, lag = 1))
relh.df.dn$Rc.F.C.dn  <- cumsum(diff(wings.dn.means[wings.dn.means$Sex == "F" & 
                                                      wings.dn.means$rep == "C",]$gen_mean,lag = 1))
relh.df.dn$Rc.M.A.dn  <- cumsum(diff(wings.dn.means[wings.dn.means$Sex == "M" & 
                                                      wings.dn.means$rep == "A",]$gen_mean,lag = 1))
relh.df.dn$Rc.M.B.dn  <- cumsum(diff(wings.dn.means[wings.dn.means$Sex == "M" & 
                                                      wings.dn.means$rep == "B",]$gen_mean, lag = 1))
relh.df.dn$Rc.M.C.dn <- cumsum(diff(wings.dn.means[wings.dn.means$Sex == "M" & 
                                                     wings.dn.means$rep == "C",]$gen_mean, lag = 1))


relh.df.dn$Sc.F.B.dn <-cumsum(wings.dn.means[wings.dn.means$Sex == "F" & wings.dn.means$rep == "B",]$S)[1:6]
relh.df.dn$Sc.F.C.dn <-cumsum(wings.dn.means[wings.dn.means$Sex == "F" & wings.dn.means$rep == "C",]$S)[1:6]
relh.df.dn$Sc.M.A.dn <-cumsum(wings.dn.means[wings.dn.means$Sex == "M" & wings.dn.means$rep == "A",]$S)[1:6]
relh.df.dn$Sc.M.B.dn <-cumsum(wings.dn.means[wings.dn.means$Sex == "M" & wings.dn.means$rep == "B",]$S)[1:6]
relh.df.dn$Sc.M.C.dn <-cumsum(wings.dn.means[wings.dn.means$Sex == "M" & wings.dn.means$rep == "C",]$S)[1:6]

#adding gen 
relh.df.dn$gen <- seq(1:6)

relh.df.dn

#Now to make this long. 

relh.long.dn <- pivot_longer(relh.df.dn, cols = c("Sc.F.A.dn":"Sc.M.C.dn"),names_to = "tag", values_to = "score")

relh.long.dn <- separate(relh.long.dn, tag, into = c("stat", "sex", "rep", "treat"))

#this is dumb but now I want to fix the two stats I think 
relh.final.dn <- pivot_wider(relh.long.dn, names_from = "stat", values_from = "score")

relh.final.dn

#This seems ok. Everything is in the right direction. 
ggplot(relh.final.dn, aes(x = Sc, y = Rc, col = rep, shape = sex )) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  theme_classic() + 
  ylab("cumulative resopnse") +
  xlab("cumulative selection diff")

#Now to do the regression. 
#This fit is sigular
dn.mod <- lmer(Rc ~ Sc * sex + (1 |rep) , data = relh.final.dn)

#dn.mod_v2 <- glmmTMB(Rc ~ Sc*sex + (1 |rep) , data = relh.final.dn)

summary(dn.mod)

#The slope is the realized haritability but I am unsure exactly how to get that. 

#Defult is 95% conf intrival
#2 fixed  Sc           0.283  CI: 0.214   0.352    

#fixed  Sc           0.394     CI:    0.292    0.496 
tidy(dn.mod, effects="fixed", conf.int = T)

#control liniages 
wings.cn <- filter(wings, line == "CR")

#getting means of each generation (total)

wings.cn.gen.mean <- (wings.cn %>%
                        group_by(interaction(gen0, Sex, rep)) %>%
                        summarize(gen_mean = mean(ds)) %>%
                        separate("interaction(gen0, Sex, rep)", 
                                 into = c("gen0", "Sex", "rep"))
)

wings.cn.gen.mean

wings.cn.sel.mean <- (wings.cn %>%
                        filter(isSelected == 1) %>%
                        group_by(interaction(gen0, Sex, rep)) %>%
                        summarize(sel_mean = mean(ds)) %>%
                        separate("interaction(gen0, Sex, rep)", 
                                 into = c("gen0", "Sex", "rep"))
)

wings.cn.sel.mean


wings.cn.means <- left_join(wings.cn.gen.mean, wings.cn.sel.mean)

#I don't know if I will need this later, but now I have this. 
wings.cn.means$gen0.sex.rep <- paste(wings.cn.means$gen0, wings.cn.means$Sex, wings.cn.means$rep, sep = ".")
wings.cn.means

wings.cn.means$S <- wings.cn.means$sel_mean - wings.cn.means$gen_mean
wings.cn.means

#getting the diffrences between generations within line. 
#I can't find a good way to get this so I'm just going to do it by hand and figure out the pretty way later (or not)

#Going to fill a df so that I can plot and do the regression after. 

Rc.F.A.cn <- cumsum(diff(wings.cn.means[wings.cn.means$Sex == "F" & 
                                          wings.cn.means$rep == "A",]$gen_mean, 
                         lag = 1))

Sc.F.A.cn  <- cumsum(wings.cn.means[wings.cn.means$Sex == "F" & wings.cn.means$rep == "A",]$S)[1:6]

relh.df.cn <- data.frame(Sc.F.A.cn , Rc.F.A.cn )

#filling in the rest. 
relh.df.cn$Rc.F.B.cn  <- cumsum(diff(wings.cn.means[wings.cn.means$Sex == "F" & 
                                                      wings.cn.means$rep == "B",]$gen_mean, lag = 1))
relh.df.cn$Rc.F.C.cn  <- cumsum(diff(wings.cn.means[wings.cn.means$Sex == "F" & 
                                                      wings.cn.means$rep == "C",]$gen_mean,lag = 1))
relh.df.cn$Rc.M.A.cn  <- cumsum(diff(wings.cn.means[wings.cn.means$Sex == "M" & 
                                                      wings.cn.means$rep == "A",]$gen_mean,lag = 1))
relh.df.cn$Rc.M.B.cn  <- cumsum(diff(wings.cn.means[wings.cn.means$Sex == "M" & 
                                                      wings.cn.means$rep == "B",]$gen_mean, lag = 1))
relh.df.cn$Rc.M.C.cn <- cumsum(diff(wings.cn.means[wings.cn.means$Sex == "M" & 
                                                     wings.cn.means$rep == "C",]$gen_mean, lag = 1))


relh.df.cn$Sc.F.B.cn <-cumsum(wings.cn.means[wings.cn.means$Sex == "F" & wings.cn.means$rep == "B",]$S)[1:6]
relh.df.cn$Sc.F.C.cn <-cumsum(wings.cn.means[wings.cn.means$Sex == "F" & wings.cn.means$rep == "C",]$S)[1:6]
relh.df.cn$Sc.M.A.cn <-cumsum(wings.cn.means[wings.cn.means$Sex == "M" & wings.cn.means$rep == "A",]$S)[1:6]
relh.df.cn$Sc.M.B.cn <-cumsum(wings.cn.means[wings.cn.means$Sex == "M" & wings.cn.means$rep == "B",]$S)[1:6]
relh.df.cn$Sc.M.C.cn <-cumsum(wings.cn.means[wings.cn.means$Sex == "M" & wings.cn.means$rep == "C",]$S)[1:6]

#adding gen 
relh.df.cn$gen <- seq(1:6)

#Now to make this long. 

relh.long.cn <- pivot_longer(relh.df.cn, cols = c("Sc.F.A.cn":"Sc.M.C.cn"),names_to = "tag", values_to = "score")

relh.long.cn <- separate(relh.long.cn, tag, into = c("stat", "sex", "rep", "treat"))

#this is dumb but now I want to fix the two stats I think 
relh.final.cn <- pivot_wider(relh.long.cn, names_from = "stat", values_from = "score")

relh.final.cn

#This seems ok. Everything is in the right direction. 
ggplot(relh.final.cn, aes(x = Sc, y = Rc, col = rep, shape = sex )) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  theme_classic() + 
  ylab("cumulative resopnse") +
  xlab("cumulative selection diff")

#Now to do the regression. 
#This gives a warning about everything being on a diffrent scale. Because it is. 

cn.mod <- lmer(Rc ~ Sc * sex + (1 |rep) , data = relh.final.cn)
#cn.mod_v2 <- glmmTMB(Rc ~ Sc * sex + (1 |rep) , data = relh.final.cn)

summary(cn.mod)

#The slope is the realized haritability but I am unsure exactly how to get that. 
#Defult is 95% conf intrival

#fixed  cond      Sc      1.66    CI:   6.95e-1  2.63 
#tidy(cn.mod, effects="fixed", conf.int = T)

#plotting all together

all.relh <- rbind( relh.final.up, relh.final.dn)

#this looks dumb because they are all on diffrent scales. I will not be using this. 
#remove controls from this. 
ggplot(all.relh[all.relh$sex == "F",], aes(x = abs(Sc), y = Rc, col = treat, shape = rep )) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  theme_classic() + 
  ylab("cumulative resopnse") +
  xlab("cumulative selection diff")

ggplot(all.relh[all.relh$sex == "M",], aes(x = abs(Sc), y = Rc, col = treat, shape = rep )) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  theme_classic() + 
  ylab("cumulative resopnse") +
  xlab("cumulative selection diff")

#fixing because I dont want to play with ggplot anymore 
all.relh$selection <- all.relh$treat


png("../Figures/ds_h2_plot.png", width =2000, height = 2000, units = "px",res = 300)

ggplot(all.relh, aes(x = abs(Sc), y = Rc, col = selection, shape = rep, linetype = sex )) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  theme_classic() + 
  ylab("Cumulative Resopnse") +
  xlab("Cumulative Selection Differential") 

dev.off()

