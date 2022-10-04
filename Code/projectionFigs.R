#I can't make these look the way I want using base so I am going to convert to ggplot.
#plus, because they are all ggobjects, they will go directly into cowplot. 

library(ggplot2)
library(cowplot)
library(latex2exp)

# two raw vectors
x <-  seq(0, 1.5, 0.1)
y <-  seq(0, 1.5, 0.1)

dat1 <- data.frame(x  = seq(0, 1.5, 0.1),  y =  seq(0, 1.5, 0.1))

x_vector <- c(1, 1)
length_x <- sqrt(t(x_vector) %*% x_vector)
x_vector_scaled <- x_vector/length_x

x1dat <- data.frame(x  = c(0, 1/length_x),  y =  c(0, 1/length_x))


y_vector <- c(0.5, 0.25)
length_y <- sqrt(t(y_vector) %*% y_vector)
y_vector_scaled <- y_vector/length_y

dsdat <- data.frame(x  = c(0, 0.5/length_y),  y =  c(0, 0.25/length_y))


#Now I want to make the plot 

gmaxDs <- ggplot(x1dat, aes(x, y)) + 
  geom_line(arrow = arrow(), col = "blue") + 
  annotate("text", x = x1dat[2,1], y = x1dat[2,2] + 0.03, label = expression(paste(italic(g))[paste(italic(max))]), col = "blue", size = 8) + 
  geom_line(aes(dsdat$x, dsdat$y), arrow = arrow(), col = "red") +
  annotate("text", x = dsdat[2,1], y = dsdat[2,2] + 0.03, label = expression(paste(italic("ds"))), col = "red", size = 8) +
  annotate("text", x = 0.24, y = 0.18, label = expression(theta), col = "black", size = 10) + 
  annotate("text", x = 0.28, y = 0.7, label =  TeX(r"($ \textit{r} = \frac{\textbf{ds} \cdot \textbf{g_{max}}}{||\textbf{ds}|| \times ||\textbf{g_{max}}||}$)"), size = 8) + 
  xlab(expression(x[1])) +
  ylab(expression(x[2])) +
  theme_classic() +
  xlim(0, 1) +
  ylim(0, 1) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size=20))


gmaxDs

#Now onto the projection/shape score figure 

# x <-  seq(0, 1.5, 0.1)
# y <-  seq(0, 1.5, 0.1)
# 
# newds <- data.frame(x = seq(0, 1.5, 0.1), y = seq(0, 1.5, 0.1))

ds_new <- c(1.2, 0.6)
length_ds_new <- sqrt(t(y_vector) %*% y_vector)
ds_new_scaled <- y_vector/length_y

dsnewdat <- data.frame(x  = c(0, ds_new[1]),  y =  c(0, ds_new[2]))

#making the two vecs
y1_vector <- c(0.8, 1)
length_y1 <- sqrt(t(y1_vector) %*% y1_vector)
y1dat <- data.frame(x  = c(0, 0.8),  y =  c(0, 1))

proj_veckp = as.vector( (y1_vector %*% ds_new) /  norm(ds_new, type = "2")^2 ) * ds_new
projdf <- data.frame(x = c(y1dat[2,1], proj_veckp[1]), y = c(y1dat[2,2],  proj_veckp[2]))


y2_vector <- c(0.22, 0.85)
length_y2 <- sqrt(t(y2_vector) %*% y2_vector)
y2dat <- data.frame(x  = c(0, 0.22),  y =  c(0, 0.85))
     
proj_vec2 = as.vector( (y2_vector %*% ds_new) /  norm(ds_new, type = "2")^2 ) * ds_new   
projdf2 <- data.frame(x = c(y2dat[2,1], proj_vec2[1]), y = c(y2dat[2,2], proj_vec2[2]))


shapescore <- ggplot(dsnewdat, aes(x, y)) + 
  geom_line(arrow = arrow(), col = "red", size = 1.5) + 
  annotate("text", x = dsnewdat[2,1], y = dsnewdat[2,2] + 0.05, label = expression(paste(italic("ds"))), col = "red", size = 8) +
  geom_line(data = y1dat, aes(x, y), arrow = arrow(), col = "black", size = 1.5) + 
  annotate("text", x = y1dat[2,1], y = y1dat[2,2] + 0.05, label = expression(y[1]), col = "black", size = 10) +
  geom_line(data = projdf, aes(x, y), arrow = arrow(), col = "black", linetype = "dashed") + 
   geom_line(data = y2dat, aes(x, y), arrow = arrow(), col = "darkgray", size = 1.5) +
   annotate("text", x = y2dat[2,1], y = y2dat[2,2] + 0.05, label = expression(y[2]), col = "darkgrey", size = 10) +
  geom_line(data = projdf2, aes(x, y), arrow = arrow(), col = "darkgray", linetype = "dashed") + 
  annotate("text", x = 0.7, y = 0.1, label =  TeX(r"($ || proj_{\textbf{ds}} \textbf{y} || = \frac{\textbf{y} \cdot \textbf{ds}}{ ||\textbf{ds}||}$)"), size = 8) +
  xlab(expression(x[1])) +
  ylab(expression(x[2])) +
  theme_classic() +
  xlim(0, 1.25) +
  ylim(0, 1.25) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size=20))

png("../Figures/shapeScorefig.png")
shapescore
dev.off()

#####
#Now for the same figure with gmax 

# y1_vector <- c(0.8, 1)
# length_y1 <- sqrt(t(y1_vector) %*% y1_vector)
# 
# y2_vector <- c(0.22, 0.85)
# length_y2 <- sqrt(t(y2_vector) %*% y2_vector)

gmaxvec <- c(1.1, 0.9)
length_gmaxvec <- sqrt(t(gmaxvec) %*% gmaxvec)
gmaxvec_scaled <- gmaxvec /length_gmaxvec

gmaxdat <- data.frame(x = c(0,gmaxvec[1]), y = c(0, gmaxvec[2]))

proj_vec3 = as.vector( (y1_vector %*% gmaxvec) /  norm(gmaxvec, type = "2")^2 ) * gmaxvec  
projdf3 <- data.frame(x = c(y1dat[2,1], proj_vec3[1]), y = c(y1dat[2,2], proj_vec3[2]))

proj_vec4 = as.vector( (y2_vector %*% gmaxvec) /  norm(gmaxvec, type = "2")^2 ) * gmaxvec  
projdf4 <- data.frame(x = c(y2dat[2,1], proj_vec4[1]), y = c(y2dat[2,2], proj_vec4[2]))

gmaxproj <- ggplot(gmaxdat, aes(x, y)) + 
  geom_line(arrow = arrow(), col = "blue", size = 1.5) + 
  annotate("text", x = gmaxdat[2,1], y = gmaxdat[2,2] + 0.05, label = expression(paste(italic(g))[paste(italic(max))]) , col = "blue", size = 8) +
  geom_line(data = y1dat, aes(x, y), arrow = arrow(), col = "black", size = 1.5) + 
  annotate("text", x = y1dat[2,1], y = y1dat[2,2] + 0.05, label = expression(y[1]), col = "black", size = 10) +
  geom_line(data = projdf3, aes(x, y), arrow = arrow(), col = "black", linetype = "dashed") + 
  geom_line(data = y2dat, aes(x, y), arrow = arrow(), col = "darkgray", size = 1.5) +
  annotate("text", x = y2dat[2,1], y = y2dat[2,2] + 0.05, label = expression(y[2]), col = "darkgrey", size = 10) +
  geom_line(data = projdf4, aes(x, y), arrow = arrow(), col = "darkgray", linetype = "dashed") + 
  xlab(expression(x[1])) +
  ylab(expression(x[2])) +
  theme_classic() +
  xlim(0, 1.25) +
  ylim(0, 1.25) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size=20))

gmaxproj



###############
#and the final pannel with the points to show what that corralation looks like 

#first I need a bunch of points, each with a PC1 score and ds shape score. 
#I should look at the actual mean and sd of the data and fix this 

#Pulled these from the wild wing estimates 
sig <- matrix(c(0.1, 0.1,0.001,0.2),2,2)

points <- data.frame(mvrnorm(n =100, c(2, 5), sig))
  
  # data.frame(gmax = rnorm(100, 4.7007e-20, 0.009032069), ds = rnorm(100, -5.16757e-19, 0.008005567))

scatterplot <- ggplot(points, aes(x = points[,1], y = points[,2])) + 
  geom_smooth(method = "lm", se = FALSE, colour = "grey") +
  geom_point() + 
  theme_classic() + 
  ylab(expression(paste(italic("ds"), " shape score"))) + 
  xlab(expression(paste(italic(g))[paste(italic(max))]))  + 
        theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size=20))

scatterplot


#Making the full figure 
library(cowplot)

projectionplot <- plot_grid(shapescore, gmaxproj, gmaxDs, scatterplot,
          labels = c("A", "B", "C", "D"))


png("../Figures/projectionFigures.png", width =2000, height = 2000, units = "px",res = 200)
projectionplot
dev.off()


