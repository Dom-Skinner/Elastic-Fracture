pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-x.pdf", useDingbats = FALSE, width = 4.5, height = 2)
dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-x.csv")
dat$K <- factor(dat$K, levels=c("one","two","three","five","interp"),labels=c("K=0.46","K=0.39","K=0.28","K=0.21","Extrapolated"))
#After some reflection, there doesn't appear to be any easy way to extract the K values from the data file, so I've just hard coded them in here...

library(ggplot2)
library(scales)     # Need the scales package
library(grid) 
library(gtable)
cbPalette <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#000000")

p <- ggplot(data=dat, aes(x=x, y= hprime,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.2) +
  geom_point(size=1) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(ylim=c(0.352,0.375),xlim=c(0,0.015),expand=FALSE) +
  ylab(expression("H'"*xi^{1/3} )) + xlab(expression(xi)) +
  theme_custom()+  
  scale_x_continuous(breaks=seq(0,0.014, by=0.002))+
  guides(col=guide_legend(ncol=3))+
  theme(axis.title = element_text(size=8),
  legend.justification=c(1,1), legend.position=c(0.75,1.06), legend.text = element_text(size=8)
  ,legend.key=element_rect(colour = NA))
   
   
grid.draw(plot_custom(p))

  
dev.off()