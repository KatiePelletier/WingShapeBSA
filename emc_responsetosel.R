library(data.table)
library(tidyverse)
library(lme4)
library(car)
library(effects)
library(emmeans)
library(glmmTMB)
library(broom.mixed)
library(boot)
library(geomorph)


#I will use Will's plotting functions here. 
source('~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R', chdir = TRUE)

source('~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R', chdir = TRUE)

projFunction <- function(x, y) {
  scalarProj <- (x %*% y) / norm(y, type = "2")
  return(scalarProj)
}

#making the figure to look at the response to selection (phenotype) for emc
#This is transformed from the landmark data from the Ian share dropbox, used the rotation matrix that is with the selection vectors. 
#I think this is mean centered and but I don't know what happened with size. 

wings <- fread("../Data/emc_selection_all_lmTransform.csv")
#head(dat)

#I am guessing this ds label is a mistake (see notes in cleanup script) but something to note if stats are wonky.
levels(as.factor(wings$sel))

#all the same. good 
levels(as.factor(wings$treat))

#case problems. 
levels(as.factor(wings$Sex))
wings$Sex <- ifelse(wings$Sex == "f", "F", wings$Sex)
levels(as.factor(wings$Sex))

levels(as.factor(wings$gen))
wings$gen <- as.numeric(gsub("Gen", "", wings$gen))
#just to check
str(wings$gen)

#reading in the ds vec 
selvec <- read.csv( "../Data/seldict_vectors.csv" )
head( selvec )
emc_vec <- as.matrix(selvec[2,3:98])
norm(emc_vec, type = "2")
#projection 
wings$emc <- projFunction(x = as.matrix(wings[,22:117]) , y = t( emc_vec ))


#Just a quick look at this
means <- (wings %>%
            group_by(interaction(wings$treat, wings$rep, wings$gen, wings$Sex))  %>%
            summarise(emc = mean(emc)) %>%
            separate("interaction(wings$treat, wings$rep, wings$gen, wings$Sex)", into = c("treat", "rep", "gen", "Sex")) )

means
means$gen <- as.numeric(means$gen)


male_means <- filter(means, Sex == "M")
female_means <- filter(means, Sex == "F")

male.treat.avg <- (male_means %>%
                    unite(treat.gen, treat, gen) %>%
                    mutate(treat.gen2 = as.factor(treat.gen)) %>%
                    group_by(treat.gen2) %>%
                    summarise(emc = mean(emc)) %>%
                    separate(treat.gen2, into = c("treat", "gen"))
)




#quick plots. 
m.response <- ggplot(male_means, aes(x = gen, y = emc, col = rep, shape = treat)) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  theme_classic() + 
  ylab("emc shape score") +
  xlab("Generation") + 
  theme(legend.position = "none") +
  scale_x_continuous(breaks = 1:7)

ggplot(female_means, aes(x = gen, y = emc, col = rep, shape = treat)) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  theme_classic() + 
  ylab("emc shape score") +
  xlab("Generation") + 
  theme(legend.position = "none") +
  scale_x_continuous(breaks = 1:7)

#now some better tests
#missing????
#wings$CS <- as.numeric(wings$CS)
wings$Sex <- as.factor(wings$Sex)
wings$treat <- as.factor(wings$treat)
wings$gen <- as.numeric(wings$gen)

wings$gen0 <- wings$gen - 1 
wings$rep <- as.factor(wings$rep)
#wings$CS_0 <- wings$CS - mean(wings$CS)


#Model for plotting. 
emc.sel.shape_mod <- lmer(emc ~  Sex + treat*gen0 + (1|treat:rep),
                         data = wings)
summary(emc.sel.shape_mod)

# emc.sel.shape_mod3 <- lmer(emc ~ Sex + treat*gen0 + (1 + gen0||treat:rep) , 
#                           data = wings)
# 
# summary(emc.sel.shape_mod3)
emc.sel.shape_mod2 <- lmer(emc ~ (Sex + treat + gen0)^3 + (1 + gen0||treat:rep) , 
                          data = wings)
summary(emc.sel.shape_mod2)
Anova(emc.sel.shape_mod2)

# crap <- emmeans(emc.sel.shape_mod, ~ treat|Sex)
# crap
# plot(crap)

crap2 <- emmeans(emc.sel.shape_mod2, ~gen0|treat)
summary(crap2)

plot(allEffects(emc.sel.shape_mod2))

allEffects(emc.sel.shape_mod2)
str(wings)

#sex does matter 
Anova(emc.sel.shape_mod)

#line:gen effect is what I really care about 

allEffects(emc.sel.shape_mod)

linegen_effect <- predictorEffect("gen0", by = c("treat", "Sex"), 
                                  emc.sel.shape_mod2)

linegen_effect2 <- data.frame(linegen_effect)

linegen_effect2$gen <- linegen_effect2$gen0 + 1

linegen_effect2$emc <- linegen_effect2$fit

ggplot(female_means, aes(x = gen, y = emc, shape = treat, col = rep)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.2) + 
  scale_colour_manual(values=c('black', 'grey46', 'grey57')) +
  geom_line(data = linegen_effect2[linegen_effect2$Sex =="F",], col = "red") + 
  theme_classic() + 
  ylab(expression(paste(italic("emc"), " shape score")))+
  xlab("Generation") + 
  theme(legend.position = "none")  + 
  scale_x_continuous(breaks = 1:7)  

#I think I actually want to keep lines seperate rather than smooth over them here. 
plot(linegen_effect)

#I actually plot this in the figure making  code so I can combine easily with the genomic plot. 


############################Realized h^2################################
#to do this I want to calculate the mean shape score in each generation for both the total population and selected individuals
#Then I want to find the diffrence between selected individuals and total population (S)
#I can then calculate the diffrence between mean in generation 2 - mean in gen 1 (R)
#From this I can get the cummulative S in each generation and plot this against cumulative R. 
#From this regression I can extimate the h^2 

#Starting with the ups. 
str(wings)
#isSelected is the T/F for which indiv are breeders for the next generation. 

wings.up <- filter(wings, treat == "UP")
str(wings.up)

#getting means of each generation (total)

wings.up.gen.mean <- (wings.up %>%
                     group_by(interaction(gen0, Sex, rep)) %>%
                     summarize(gen_mean = mean(emc)) %>%
                       separate("interaction(gen0, Sex, rep)", 
                                into = c("gen0", "Sex", "rep"))
                     )

wings.up.gen.mean

wings.up.sel.mean <- (wings.up %>%
                        filter(isSelected == 1) %>%
                        group_by(interaction(gen0, Sex, rep)) %>%
                        summarize(sel_mean = mean(emc)) %>%
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
up.mod_v2 <- glmmTMB(Rc ~ Sc*sex + (1 |rep) , data = relh.final.up)


summary(up.mod)
summary(up.mod_v2)

#The slope is the realized haritability

#Well that is a handy line of code! I finally understand why they used this all the time in QMEE *from broom.mixed package*
#Defult is 95% conf intrival

#fixed  Sc           0.381    CI:  0.290      0.471  
tidy(up.mod_v2, effects="fixed", conf.int = T)

#This is for nice plotting if I want to do that. 
# allEffects(up.mod)
# relh_effect <- predictorEffect("Sc",  up.mod)
# plot(relh_effect)
# relh_effect2 <- data.frame(relh_effect)


#Now for the down. 
wings.dn <- filter(wings, treat == "DN")
str(wings.dn)

#getting means of each generation (total)

wings.dn.gen.mean <- (wings.dn %>%
                        group_by(interaction(gen0, Sex, rep)) %>%
                        summarize(gen_mean = mean(emc)) %>%
                        separate("interaction(gen0, Sex, rep)", 
                                 into = c("gen0", "Sex", "rep"))
)

wings.dn.gen.mean

wings.dn.sel.mean <- (wings.dn %>%
                        filter(isSelected == 1) %>%
                        group_by(interaction(gen0, Sex, rep)) %>%
                        summarize(sel_mean = mean(emc)) %>%
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
dn.mod <- lmer(Rc ~ Sc  + sex + (1 |rep) , data = relh.final.dn)

dn.mod_v2 <- glmmTMB(Rc ~ Sc*sex + (1 |rep) , data = relh.final.dn)

summary(dn.mod_v2)

#The slope is the realized haritability but I am unsure exactly how to get that. 

#Defult is 95% conf intrival
#2 fixed  Sc           0.283  CI: 0.214   0.352    


#fixed  Sc           0.381    CI:  0.290      0.471  
tidy(dn.mod, effects="fixed", conf.int = T)

#control liniages 
wings.cn <- filter(wings, treat == "CR")

#getting means of each generation (total)

wings.cn.gen.mean <- (wings.cn %>%
                        group_by(interaction(gen0, Sex, rep)) %>%
                        summarize(gen_mean = mean(emc)) %>%
                        separate("interaction(gen0, Sex, rep)", 
                                 into = c("gen0", "Sex", "rep"))
)

wings.cn.gen.mean

wings.cn.sel.mean <- (wings.cn %>%
                        filter(isSelected == 1) %>%
                        group_by(interaction(gen0, Sex, rep)) %>%
                        summarize(sel_mean = mean(emc)) %>%
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
cn.mod_v2 <- glmmTMB(Rc ~ Sc * sex + (1 |rep) , data = relh.final.cn)

#can I even mean centre this data? Is that even valid? Does it even matter?
relh.final.cn$Rc_c <- relh.final.cn$Rc - mean(relh.final.cn$Rc)
relh.final.cn$Sc_c <- relh.final.cn$Sc - mean(relh.final.cn$Sc)

#this did not help. Fuck it. 
cn.mod_v3 <- lmer(Rc_c ~ Sc_c * sex + (1 |rep) , data = relh.final.cn)

summary(cn.mod)
summary(cn.mod_v2)
summary(cn.mod_v3)

#The slope is the realized haritability but I am unsure exactly how to get that. 

#Defult is 95% conf intrival

#fixed  cond      Sc      1.66    CI:   6.95e-1  2.63 
tidy(cn.mod_v2, effects="fixed", conf.int = T)

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

all.relh$treatment <- all.relh$treat

png("../Figures/emc_h2_plot.png", width =2000, height = 2000, units = "px",res = 300)
ggplot(all.relh, aes(x = abs(Sc), y = Rc, col = treat, shape = rep, linetype = sex )) + 
  geom_point() + 
  geom_line(alpha = 0.5) + 
  theme_classic() + 
  ylab("Cumulative Resopnse") +
  xlab("Cumulative Selection Differential")
dev.off()

###################vector cor between real change and selection vec###########
#for this I want to find the mean for each rep:treatment at gen 1 and 7 then take the mean of means. 
#Then I will get the diffrence between gen 1 and 7 (can scale to be 1 gen by dividing by 6? or was this a diffrent disscusion)
#I can cor this vector with the selection vector (emc and ds in this case) 
#for some CI I can boot strap this by sampling with replacement within replicate and generation (strata in the boot function should help here) However, ID says this may not behave well because this tends to do weird things when near a boundry. 


wings.shape.means <- (wings %>%
                        group_by(interaction(treat, gen0, Sex, rep)) %>%
                        summarise_at(names(.)[22:117],
                                     .funs = c(mean="mean")) %>%
                        separate("interaction(treat, gen0, Sex, rep)", 
                                 into = c("treat", "gen0", "Sex", "rep"))
)

wings.shape.means

wings.shape.means[1,]
#My code is way faster. But it is a tiny bit diffrent. why
colMeans(wings[wings$treat == "CR" & wings$gen0 == 0 & wings$Sex == "F", wings$rep == "A", 22:117])


C.gen1.f <- colMeans(wings[wings$treat == "CR" & wings$gen0 == 0 & wings$Sex == "F", 22:117])
C.gen7.f <- colMeans(wings[wings$treat == "CR" & wings$gen0 == 6 & wings$Sex == "F", 22:117])

deltaC.f <- C.gen7.f - C.gen1.f 
#-0.1337788
cor(deltaC.f, as.numeric(emc_vec))
#0.04906536
cor(deltaC.f, as.numeric(selvec[1,3:98]))

PD(deltaC.f) #0.009667443

C.gen1.m <- colMeans(wings[wings$treat == "CR" & wings$gen0 == 0 & wings$Sex == "M", 22:117])
C.gen7.m <- colMeans(wings[wings$treat == "CR" & wings$gen0 == 6 & wings$Sex == "M", 22:117])

deltaC.m <- C.gen7.m - C.gen1.m
#-0.08776019
cor(deltaC.m, as.numeric(emc_vec))

#0.1087453
cor(deltaC.m, as.numeric(selvec[1,3:98]))

PD(deltaC.m) #0.009667443


PD(deltaC.m) # 0.008376809

dn.gen1.f <- colMeans(wings[wings$treat == "DN" & wings$gen0 == 0 & wings$Sex == "F", 22:117])
dn.gen7.f <- colMeans(wings[wings$treat == "DN" & wings$gen0 == 6 & wings$Sex == "F", 22:117])

delta.dn.f <- dn.gen7.f - dn.gen1.f 
# -0.686719
cor(delta.dn.f, as.numeric(emc_vec))
#-0.6362707
cor(delta.dn.f, as.numeric(selvec[1,3:98]))
PD(delta.dn.f) # 0.02093135

dn.gen1.m <- colMeans(wings[wings$treat == "DN" & wings$gen0 == 0 & wings$Sex == "M", 22:117])
dn.gen7.m <- colMeans(wings[wings$treat == "DN" & wings$gen0 == 6 & wings$Sex == "M", 22:117])

delta.dn.m <- dn.gen7.m - dn.gen1.m
# -0.7528295
cor(delta.dn.m, as.numeric(emc_vec))
# -0.6485443
cor(delta.dn.m, as.numeric(selvec[1,3:98]))

PD(delta.dn.m) # 0.02029452

up.gen1.f <- colMeans(wings[wings$treat == "UP" & wings$gen0 == 0 & wings$Sex == "F", 22:117])
up.gen7.f <- colMeans(wings[wings$treat == "UP" & wings$gen0 == 6 & wings$Sex == "F", 22:117])

delta.up.f <- up.gen7.f - up.gen1.f
# 0.7530443
cor(delta.up.f, as.numeric(emc_vec))
# 0.7769413
cor(delta.up.f, as.numeric(selvec[1,3:98]))

PD(delta.up.f) #0.0425766

up.gen1.m <- colMeans(wings[wings$treat == "UP" & wings$gen0 == 0 & wings$Sex == "M", 22:117])
up.gen7.m <- colMeans(wings[wings$treat == "UP" & wings$gen0 == 6 & wings$Sex == "M", 22:117])

delta.up.m <- up.gen7.m - up.gen1.m
# 0.6911699
cor(delta.up.m, as.numeric(emc_vec))
# 0.7036082
cor(delta.up.m, as.numeric(selvec[1,3:98]))

PD(delta.up.m) #0.03968325

####################blurs for paper 

#I want to do a wing blur comparing the variation between up and down wings. 

crapwingsf <- (wings %>% 
                 filter(treat == "DN" | treat == "UP") %>%
                 filter(Sex == "F") %>%
                 filter(gen0 == 6))

crapwingsf$treat <- droplevels(crapwingsf$treat)
#needs a mean shape added back for plotting. effect size is made smaller for plotting reasons

#this looks so shitty and I don't know why. Works with the geomorph function?
meanshape <- as.numeric(selvec[3,3:98])/10

crapwingsf.plot <- as.matrix(crapwingsf[,22:117]) + rep(meanshape, each = nrow(crapwingsf))

WingBlur3(crapwingsf.plot)

# png("../Figures/emcSel_F_wingBlur.png")
# WingBlur3(crapwingsf.plot, grouping_var = crapwingsf$treat, groups = T)
# dev.off()
# 
# crapwingsm <- (wings %>% 
#                  filter(treat == "DN" | treat == "UP") %>%
#                  filter(Sex == "M") %>%
#                  filter(gen0 == 6))
# 
# crapwingsm$treat <- droplevels(crapwingsm$treat)
# 
# png("../Figures/emcSel_M_wingBlur.png")
# WingBlur3(as.matrix(crapwingsm[,16:111]), grouping_var = crapwingsf$line, groups = T)
# dev.off()


#############Shape change figures for the paper##############

source("~/Dropbox/KatiePelletier/KP_geomorphwingfunctions.R")

#this magnitude is WAY too big and you can't see anything.
meanshape <- as.numeric(selvec[3,3:98])/10

#I need to add a mean shape back to these landmarks so they will plot
wings.withmeans <- as.matrix(wings[,22:117]) + rep(meanshape, each = nrow(wings))

cord <- as.matrix(wings.withmeans)
shape <- arrayspecs(cord, 48, 2)
wing.dat <- geomorph.data.frame(shape = shape,
                                treat = wings$treat,
                                rep = wings$rep,
                                gen = wings$gen,
                                sex = wings$Sex)

group <- factor(paste(wing.dat$gen, wing.dat$treat, wing.dat$sex))
levels(group)

new.coords <- coords.subset(wing.dat$shape, group = group)
names(new.coords) # see the list levels
# group shape means

group_means <- lapply(new.coords, mshape)

png("../Figures/emc_selection_shapechange_UP_F1v7_15xmag.png")
plotRefToTarget(group_means[["1 UP F"]],
                group_means[["7 UP F"]],
                links = wing.links, method = "points", mag = 15,
                gridPars=wing.spec )
dev.off()

png("../Figures/emc_selection_shapechange_DN_F1v7_15xmag.png")
plotRefToTarget(group_means[["1 DN F"]],
                group_means[["7 DN F"]],
                links = wing.links, method = "points", mag = 15,
                gridPars=wing.spec )
dev.off()

png("../Figures/emc_selection_shapechange_CR_F1v7_15xmag.png")
plotRefToTarget(group_means[["1 CR F"]],
                group_means[["7 CR F"]],
                links = wing.links, method = "points", mag = 15,
                gridPars=wing.spec )
dev.off()

png("../Figures/emc_selection_shapechange_UP_M1v7_15xmag.png")
plotRefToTarget(group_means[["1 UP M"]],
                group_means[["7 UP M"]],
                links = wing.links, method = "points", mag = 10,
                gridPars=wing.spec )
dev.off()

png("../Figures/emc_selection_shapechange_DN_M1v7_15xmag.png")
plotRefToTarget(group_means[["1 DN M"]],
                group_means[["7 DN M"]],
                links = wing.links, method = "points", mag = 15,
                gridPars=wing.spec )
dev.off()

png("../Figures/emc_selection_shapechange_CR_M1v7_5xmag.png")
plotRefToTarget(group_means[["1 CR M"]],
                group_means[["7 CR M"]],
                links = wing.links, method = "points", mag = 5,
                gridPars=wing.spec )
dev.off()



##########################################################################
#Now to take a shot at the bootstrapping. 
#reverse order as above, because controls are boring!! 

#just working with one sex at a time. 
wings.f <- filter(wings, Sex == "F")

#I need to create a treat:rep:gen variable 
# names(wings.f)
# 
# wings.f$group <- paste(wings.f$rep, wings.f$gen0, 
#                        sep = ".")
# 
# wings.f$group <- as.factor(wings.f$group)
# levels(wings.f$group)
# 
# ds.cor <- function(x){
#   gen1mean <- colMeans(x[x$gen0 == 0,22:117])
#   gen7mean <- colMeans(x[x$gen0 == 6,22:117])
#   diffmean <- gen7mean - gen1mean
#   cor(diffmean, as.numeric(selvec[1,3:98]))
# }
# 
# emc.cor <- function(x){
#   gen1mean <- colMeans(x[x$gen0 == 0,22:117])
#   gen7mean <- colMeans(x[x$gen0 == 6,22:117])
#   diffmean <- gen7mean - gen1mean
#   cor(x, as.numeric(emc_vec))
# }
# 
# f.up.boot <- boot(wings.f[wings.f$treat == "UP"], emc.cor, R = 10,
#                   strata = wings.f$group)

