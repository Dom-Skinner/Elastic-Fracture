pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-KII.pdf", useDingbats = FALSE, width = 3, height = 2)

dat1 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-KII-1.csv")
dat1$lambda <- factor(dat1$lambda, levels=c("0.02","0.04","0.06","0.08"),
  labels=c("lambda=0.02","lambda=0.04","lambda=0.06","lambda=0.08"))


library(ggplot2)
cbPalette1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
cbPalette2 <- c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(data=dat1, aes(x=L, y= KI,group=lambda,shape=lambda,color=lambda))  +
  geom_line(linetype="solid", size=0.2) +
  geom_point(size=1.5) +
  scale_colour_manual(values=cbPalette1,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18), name  = NULL) +
  coord_cartesian(ylim=c(0.8,3.2),xlim=c(0.8,3.4),expand=FALSE) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
  legend.justification=c(0,1), legend.position=c(0,1),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
  legend.text = element_text(size=8,family="Palatino"),legend.direction = "vertical")
 

dat2 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-KII-2.csv")
dat2$lambda <- factor(dat2$lambda, levels=c("0.02","0.04","0.06","0.08"),
  labels=c("lambda=0.02","lambda=0.04","lambda=0.06","lambda=0.08"))


ggplot(data=dat2, aes(x=L, y= KII,group=lambda,shape=lambda,color=lambda))  +
  geom_line(linetype="solid", size=0.2) +
  geom_point(size=1.5) +
  scale_colour_manual(values=cbPalette1,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18), name  = NULL) +
  coord_cartesian(ylim=c(-0.08,0.81),xlim=c(0.8,3.4),expand=FALSE) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
  legend.justification=c(1,1), legend.position=c(1,1),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
  legend.text = element_text(size=8,family="Palatino"),legend.direction = "vertical")
 

 
 
dat3 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-KII-3.csv")
dat3$L <- factor(dat3$L, levels=c("0.8","1.2","1.6","2"),labels=c("L = 0.8","L = 1.2","L = 1.6","L = 2"))

ggplot(data=dat3, aes(x=lambda, y= KI,group=L,shape=L,color=L))  +
  geom_line(linetype="solid", size=0.2) +
  geom_point(size=1) +
  scale_colour_manual(values=cbPalette2,name  = NULL) +
  scale_shape_manual(values = c(4,5,6,8), name  = NULL) +
  coord_cartesian(ylim=c(0.8,3.2),xlim=c(0,0.085),expand=FALSE) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
  legend.justification=c(0,1), legend.position=c(0,1),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
  legend.text = element_text(size=8,family="Palatino"),legend.direction = "vertical")
 
dat4 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-KII-4.csv")
dat4$L <- factor(dat4$L, levels=c("0.8","1.2","1.6","2"),labels=c("L = 0.8","L = 1.2","L = 1.6","L = 2"))


ggplot(data=dat4, aes(x=lambda, y= KII,group=L,shape=L,color=L))  +
  geom_line(linetype="solid", size=0.2) +
  geom_point(size=1) +
  scale_colour_manual(values=cbPalette2,name  = NULL) +
  scale_shape_manual(values = c(4,5,6,8), name  = NULL) +
  coord_cartesian(ylim=c(-0.08,0.81),xlim=c(0,0.085),expand=FALSE) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
  legend.justification=c(0,1), legend.position=c(0,1),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
  legend.text = element_text(size=8,family="Palatino"),legend.direction = "vertical")
 

  
dev.off()