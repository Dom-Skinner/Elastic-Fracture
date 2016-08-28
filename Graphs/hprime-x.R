pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-x.pdf", useDingbats = FALSE, width = 5.2, height = 2.5)
dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-x.csv")
dat$K <- factor(dat$K, levels=c("one","two","three","five","interp"),labels=c("K=0.46","K=0.39","K=0.28","K=0.21","Interpolated"))
#After some reflection, there doesn't appear to be any easy way to extract the K values from the data file, so I've just hard coded them in here...

library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data=dat, aes(x=x, y= hprime,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.2) +
  geom_point(size=0.8) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(ylim=c(0.352,0.375),xlim=c(0,0.015),expand=FALSE) +
  ylab(expression(Hprime *xi^{1/3} )) + xlab(expression(xi)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
  legend.justification=c(1,1), legend.position=c(1,1),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
  legend.text = element_text(size=8,family="Palatino"),legend.direction = "horizontal")+
  annotate("rect",xmin = 0.0117, xmax = 0.0145, ymin = 0.353, ymax = 0.355,color="#808080",fill="#F0F0F0", size=0.25)+
  annotate("text",label = "n=465, xend=819", x = 0.0131, y = 0.354 , size = 2.5,family="Palatino")
  
  
dev.off()