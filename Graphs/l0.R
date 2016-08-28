pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/l0.pdf", useDingbats = FALSE, width = 4.3, height = 2.3)
dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/l0.csv")
dat$type <- factor(dat$type, levels=c("Kval","linear","quadratic"),
  labels=c("Numerical values","Linear extrapolation","Quadratic extrapolation"))


library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data=dat, aes(x=K, y= lambda,group=type,shape=type,color=type,linetype=type))  +
  geom_line( size=0.4) +
  geom_point(size=1.5) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,26,26), name  = NULL) +
  # Note that shape size 26 does not actually exist, this is just a hack
  scale_linetype_manual(values=c(0,1,2), name  = NULL) +
  coord_cartesian(ylim=c(0.0569,0.0595),xlim=c(0,0.3),expand=FALSE) +
  ylab(expression(lambda)) + xlab(expression(K[I]^u)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
  legend.justification=c(1,1), legend.position=c(1,1),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
  legend.text = element_text(size=8,family="Palatino")) +
  annotate("rect",xmin = 0.015, xmax = 0.125, ymin = 0.05708, ymax = 0.0577,color="#808080",fill="#F0F0F0", size=0.25)+
  annotate("text",label = c("Estimated percentage error", "due to extrapolation = 0.002%","n=524, xend=846"), 
  x = c(0.066,0.07,0.07), y = c(0.0576,0.0574,0.0572)
  , size = 2.5,family="Palatino")
  
  
dev.off()