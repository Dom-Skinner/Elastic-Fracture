pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-x.pdf", useDingbats = FALSE, width = 5.2, height = 2.5)
dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-x.csv")
dat$K <- factor(dat$K, levels=c("one","two","three","four","five","interp"),labels=c("K=0.28","K=0.25","K=0.21","K=0.17","K=0.12","Interpolated"))
#After some reflection, there doesn't appear to be any easy way to extract the K values from the data file, so I've just hard coded them in here...

library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data=dat, aes(x=x, y= hprime,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.5) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,8,4), name  = NULL) +
  coord_cartesian(ylim=c(0.351,0.365),xlim=c(0,0.001),expand=FALSE) +
  ylab(expression(Hprime *xi^{1/3} )) + xlab(expression(xi)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=10,family="Palatino"), axis.text  = element_text(size=6),
  legend.justification=c(1,1), legend.position=c(1,1),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
  legend.text = element_text(size=10,family="Palatino"),legend.direction = "horizontal")

  
  
dev.off()