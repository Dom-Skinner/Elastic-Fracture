pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-p-x-full.pdf", useDingbats = FALSE, width = 6, height = 2)
dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-p-x-full-1.csv")
dat$K <- factor(dat$K, labels=c("K=0.74","K=1.82"))

library(ggplot2)
  library(scales)     # Need the scales package

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p1 <- ggplot(data=dat, aes(x=x, y= hprime,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.2) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(xlim=c(0,2),ylim=c(0.8,2.5),expand=FALSE) +
  ylab(expression(Hprime )) + xlab(expression(xi)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),legend.position="none")
  
  
p2 <- ggplot(data=dat, aes(x=x, y= gprime,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.2) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(xlim=c(0,2),ylim=c(0.4,1.3),expand=FALSE) +
  ylab(expression(Gprime )) + xlab(expression(xi)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
  legend.justification=c(1,1), legend.position=c(1,1),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
  legend.text = element_text(size=8,family="Palatino"))
  

  
dat2 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-p-x-full-2.csv")
dat2$K <- factor(dat2$K, labels=c("K=0.74","K=1.82"))  


p3 <- ggplot(data=dat2, aes(x=z, y= pressure,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.2) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(xlim=c(0,2),ylim=c(-3,0),expand=FALSE) +
  ylab(expression(Pi )) + xlab(expression(xi)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
  legend.position="none")
  
  


  
p4 <- ggplot(data=dat, aes(x=x, y= hprime,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(xlim=c(0.0005,6),ylim=c(0.8,7),expand=FALSE) +
  ylab(expression(log(Hprime) )) + xlab(expression(log(xi))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),legend.position="none")



p5 <- ggplot(data=dat, aes(x=x, y= gprime,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(xlim=c(0.016,6),ylim=c(0.45,1.7),expand=FALSE) +
  ylab(expression(log(Gprime) )) + xlab(expression(log(xi))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
  legend.position="none")
    
multiplot(p1,p2,p3, cols=3)
dev.off()
pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-p-x-full-2.pdf", useDingbats = FALSE, width = 5, height = 2)

multiplot(p4,p5, cols=2)

dev.off()