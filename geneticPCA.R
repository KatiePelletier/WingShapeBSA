#trying to run SMARTPCA for the paper. Using the smartsnp package. 

library(smartsnp)
library(tidyverse)

crap <- fread("../Data/wildpops_controlPool.snp")

# pathToGenoFile = system.file("extdata", "dataSNP", package = "smartsnp")
# 
# my_groups <- as.factor(c(rep("A", 50), rep("B", 50))); cols = c("red", "blue")
# 
# groups <- c(1,2,3,4)
# 
# mvaR <- smart_mva(snp_data = crap, sample_group = groups)
# mvaR$pca$pca.eigenvalues # extract PCA eigenvalues
# mvaR$pca$pca.snp_loadings # extract principal coefficients (SNP loadings)
# mvaR$pca$pca.sample_coordinates # extract PCA principal components (sample position in PCA space)
# 
# plot(mvaR$pca$pca.sample_coordinates[,c("PC1","PC2")], cex = 2,
#      pch = 19, col = cols[my_groups], main = "genotype smartpca")
# legend("topleft", legend = levels(my_groups), cex = 1,
#        pch = 19, col = cols, text.col = cols)

#I think I want to scale these so every SNP 'votes' equally? I think that is what I am doing there? Ask Ian?
geneticpca <- prcomp(t(crap), center = TRUE, scale = TRUE)
shit <- data.frame(geneticpca[["x"]])
shit$pop <- c("CMO", "FVW13", "FVW14", "PHO")

png("morecrap_forID.png")
ggplot(shit, aes(x = PC1, y = PC2, col = pop)) + 
  geom_point()+ 
  labs(col = "Population")
dev.off()

ggplot(shit, aes(x = PC3, y = PC4, col = pop)) + 
  geom_point()
