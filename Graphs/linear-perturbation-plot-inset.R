dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/linear-perturbation-plot.csv")
dat$n <- factor(dat$n, levels=c("351","526","701"),labels=c("n=350","n=525","n=700"))

library(ggplot2)
library(scales) 
library(grid) 
library(gtable)

cbPalette <- c("#e41a1c","#377eb8","#4daf4a")

p1 <- ggplot(data=dat, aes(x=x, y= Htilde,group=n,shape=n,color=n))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(ylim=c(0.575,0.59),xlim=c(0,0.015),expand=FALSE) +
  scale_x_continuous(breaks=seq(0,0.015, by=0.003))+
  scale_y_continuous(breaks=seq(0.575,0.59, by=0.005))+
  ylab(expression(tilde(H)*xi^-s )) + xlab(expression(xi)) +
  theme_custom() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"),
  legend.justification=c(0,1), legend.position=c(0,1), legend.text = element_text(size=8,family="Palatino")
  ,legend.key=element_rect(colour = NA))

  
p2 <-ggplot(data=dat, aes(x=x, y= Htilde,group=n,shape=n,color=n))  +
    geom_line(linetype="solid", size=0.3) +
  geom_point(size=1) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(ylim=c(0.5755,0.5777),xlim=c(0,0.00014),expand=FALSE) +
  theme_custom() +
  theme(axis.text  = element_text(size=6),legend.position="none",
  axis.title.x=element_blank(),axis.title.y=element_blank())

g1 = plot_custom(p1)
g2 = plot_custom(p2) 
    
p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() +
     coord_cartesian(xlim=c(0,500),ylim=c(0,220),expand=FALSE) +
     theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())

pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/linear-perturbation-plot.pdf", useDingbats = FALSE, width = 5, height = 2.2)
p + annotation_custom(grob = g1, xmin = 0, xmax = 500, ymin = 0, ymax = 220) +
annotation_custom(grob = g2, xmin = 310, xmax = 478, ymin = 37, ymax = 140)
 dev.off()