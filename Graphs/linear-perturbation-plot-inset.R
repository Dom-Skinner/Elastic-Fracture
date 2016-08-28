pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/linear-perturbation-plot.pdf", useDingbats = FALSE, width = 5.2, height = 2.5)
dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/linear-perturbation-plot.csv")
dat$n <- factor(dat$n, levels=c("350","524","815"),labels=c("n=350, xend=873","n=524, xend=846","n=815, xend=846"))

library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data=dat, aes(x=x, y= Htilde,group=n,shape=n,color=n))  +
  geom_line(linetype="solid", size=0.2) +
  geom_point(size=0.8) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(ylim=c(0.575,0.59),xlim=c(0,0.015),expand=FALSE) +
  ylab(expression(Htilde*xi^-s )) + xlab(expression(xi)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
  legend.justification=c(0,1), legend.position=c(0,1),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
  legend.text = element_text(size=8,family="Palatino"),legend.direction = "vertical")
 dev.off()
 
 pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/linear-perturbation-plot-inset.pdf", 
 useDingbats = FALSE, width = 2, height = 1.5)

 ggplot(data=dat, aes(x=x, y= Htilde,group=n,shape=n,color=n))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(ylim=c(0.575,0.577),xlim=c(0,0.00045),expand=FALSE) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0.1,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),legend.position="none",
  axis.title.x=element_blank(),axis.title.y=element_blank())
  
dev.off()