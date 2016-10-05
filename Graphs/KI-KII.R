
dat1 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-KII-1.csv")
dat1$lambda <- factor(dat1$lambda, levels=c("0.02","0.04","0.06","0.08","0.095"),
  labels=c("lambda=0.02","lambda=0.04","lambda=0.06","lambda=0.08","lambda=0.095"))


library(ggplot2)
library(scales)    
library(grid) 
library(gtable)

cbPalette1 <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")


p1 <- ggplot(data=dat1, aes(x=L, y= KI,group=lambda,shape=lambda,color=lambda))  +
  geom_line(linetype="solid", size=0.2) +
  geom_point(size=1.5) +
  scale_colour_manual(values=cbPalette1,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,8), name  = NULL) +
  coord_cartesian(ylim=c(0,3),xlim=c(0,3.4),expand=FALSE) +
  theme_custom() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),axis.title = element_text(size=8,family="Palatino"), legend.position="none")


p2 <- ggplot(data=dat1, aes(x=L, y= KII,group=lambda,shape=lambda,color=lambda))  +
  geom_line(linetype="solid", size=0.2) +
  geom_point(size=1.5) +
  scale_colour_manual(values=cbPalette1,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,8), name  = NULL) +
  coord_cartesian(ylim=c(-0.08,1.55),xlim=c(0,3.4),expand=FALSE) +
  geom_hline(yintercept=0,linetype="dashed") +
  theme_custom() +
  theme(axis.title = element_text(size=8,family="Palatino"),
  legend.justification=c(1,1), legend.position=c(1,1),legend.key=element_rect(colour = NA),
  legend.text = element_text(size=8,family="Palatino"))
 
 
dat2 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-KII-2.csv")
dat2$L <- factor(dat2$L, levels=c("0.01","0.4","0.8","1.6","3.2"),labels=c("L=0.01","L = 0.4","L = 0.8","L = 1.6","L = 3.2"))

p3 <- ggplot(data=dat2, aes(x=lambda, y= KI,group=L,shape=L,color=L))  +
  geom_line(linetype="solid", size=0.2) +
  geom_point(size=1) +
  scale_colour_manual(values=cbPalette1,name  = NULL) +
  scale_shape_manual(values = c(24,22,23,21,25), name  = NULL) +
  coord_cartesian(ylim=c(0,3),xlim=c(0,0.105),expand=FALSE) +
  
  theme_custom() +
  guides(col=guide_legend(ncol=2))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),
        legend.justification=c(1,1), legend.position=c(1.1,1.1),legend.key = element_rect(color=NA),
        legend.text = element_text(size=8,family="Palatino"))
 
p4 <- ggplot(data=dat2, aes(x=lambda, y= KII,group=L,shape=L,color=L))  +
  geom_line(linetype="solid", size=0.2) +
  geom_point(size=1) +
  scale_colour_manual(values=cbPalette1,name  = NULL) +
  scale_shape_manual(values = c(24,22,23,21,25), name  = NULL) +
  xlab(expression(lambda))+
  geom_hline(yintercept=0,linetype="dashed") +
  coord_cartesian(ylim=c(-0.08,1.55),xlim=c(0,0.105),expand=FALSE) +
  theme_custom() +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),legend.position="none")
 
g1 = plot_custom(p1)
g2 = plot_custom(p2)
g3 = plot_custom(p3)
g4 = plot_custom(p4)

p <- qplot(-20:-19,1:2) +
     coord_cartesian(xlim=c(0,600),ylim=c(0,400),expand=FALSE) +
     theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())

pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-KII.pdf", useDingbats = FALSE, width = 5.5, height = 3.5)
p + annotation_custom(grob = g1, xmin = 0   , xmax = 315, ymin = 215, ymax = 400) +
    annotation_custom(grob = g2, xmin = -9.5, xmax = 315, ymin = 0, ymax = 215) +
    annotation_custom(grob = g3, xmin = 315 , xmax = 600, ymin = 215, ymax = 400) + 
    annotation_custom(grob = g4, xmin = 315 , xmax = 600, ymin = 0, ymax = 215) 

  
dev.off()