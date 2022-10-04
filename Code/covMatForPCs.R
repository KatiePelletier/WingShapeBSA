#####doing the PCA cors in a new script because I am tired of looking at things. 

cmo.crap <- read.csv("../Data/cmo_all_first3PC.csv")

wild.crap <- read.csv("../Data/wild_all_first3PC.csv")

dgrp.crap <- read.csv("../Data/dgrp_all_first3PC.csv")

#that fucking rawnames line that I forgot to cut out. fuck it 
all.crap <- data.frame(dgrp.crap[,2:4], wild.crap[,2:4], cmo.crap[,2:4])


cor(all.crap)

cor.mat.pc <- cor(all.crap)

#write out the cor mat to make into a table. 

write.csv(cor.mat.pc, "../Tables/shape_covPCmat.csv")
