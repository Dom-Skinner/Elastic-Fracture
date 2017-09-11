dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-p-x-full-1.csv")
dat$K <- factor(dat$K, labels=c("K=0.74","K=1.82"))

library(ggplot2)
library(scales)     # Need the scales package
library(grid) 
library(gtable)

cbPalette <- c("#e41a1c","#377eb8")

p1 <- ggplot(data=dat, aes(x=x, y= 12*hprime,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.2) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16), name  = NULL) +
  coord_cartesian(xlim=c(-0.005,2.005),ylim=c(10,12*2.5),expand=FALSE) +
  ylab(expression("H'" )) + xlab(expression(xi)) +
  theme_custom() +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(legend.position="none")
  
  
p2 <- ggplot(data=dat, aes(x=x, y= 12*gprime,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.2) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16), name  = NULL) +
  coord_cartesian(xlim=c(-0.005,2.005),ylim=c(5,12*1.3),expand=FALSE) +
  ylab(expression("G'" )) + xlab(expression(xi)) +
  theme_custom() +
  theme(plot.margin = unit(c(0.2,0.49,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(legend.justification=c(1,1), legend.position=c(0.9,0.9),
  legend.text = element_text(size=8),legend.key=element_rect(colour = NA))
  

  
dat2 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-p-x-full-2.csv")
dat2$K <- factor(dat2$K, labels=c("K=0.74","K=1.82"))  


p3 <- ggplot(data=dat2, aes(x=z, y= (3/pi)*pressure,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.2) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16), name  = NULL) +
  coord_cartesian(xlim=c(-0.005,2.005),ylim=c(-(3/pi),(3/pi)*0.0025),expand=FALSE) +
  ylab(expression(P )) + xlab(expression(xi)) +
  theme_custom()+
  theme(plot.margin = unit(c(0.2,0.44,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(legend.position="none")
  
  


  
p4 <- ggplot(data=dat, aes(x=x, y= 12*hprime,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="trlb",size=0.2,short = unit(0.035, "cm"), mid = unit(0.06, "cm"), long = unit(0.1, "cm")) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16), name  = NULL) +
  coord_cartesian(xlim=c(0.00000096,40.4),ylim=c(12*0.8,12*100),expand=FALSE) +
  annotate("segment", x=10^-6, xend=10^-3, y=12*10^(3-0.617), yend=12*10^(1.5-0.617))+
  annotate("segment", x=10^-6, xend=10^-4, y=12*10^(3-1.004), yend=12*10^(2-1.004))+
  annotate("segment", x=1, xend=40, y=12, yend=12*40)+
  ylab(expression("H'" )) + xlab(expression(xi)) +
  theme_custom() + theme(axis.ticks=element_blank()) +
  theme(plot.margin = unit(c(0.11,0.41,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(legend.position="none")



p5 <- ggplot(data=dat, aes(x=x, y= 12*gprime,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="trlb",size=0.2,short = unit(0.035, "cm"), mid = unit(0.06, "cm"), long = unit(0.1, "cm")) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(xlim=c(0.002,10.3),ylim=c(12*0.45,12*4),expand=FALSE) +
  ylab(expression("G'" )) + xlab(expression(xi)) +
  theme_custom() + theme(axis.ticks=element_blank()) +
  annotate("segment", x=10^-3, xend=10^-2, y=12*10^(1.5-0.7034), yend=12*10^(1-0.7034))+
  annotate("segment", x=10^-3, xend=10^-2, y=12*10^(1.5-0.7317), yend=12*10^(1-0.7317))+
  annotate("segment", x=0.1, xend=40, y=12*0.5, yend=12*0.5)+
  theme(plot.margin = unit(c(0.11,0.4,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(legend.position="none")

    
    
p6 <- ggplot(data=dat2, aes(x=z, y= (3/pi)*pressure,group=K,shape=K,color=K))  +
  geom_line(linetype="solid", size=0.3) +
  geom_point(size=1.2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides="tb",size=0.2,short = unit(0.035, "cm"), mid = unit(0.06, "cm"), long = unit(0.1, "cm")) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,16,17,18,4), name  = NULL) +
  coord_cartesian(xlim=c(0.0000003,10.1),ylim=c(-(3/pi)*14,(3/pi)*0.035),expand=FALSE) +
  ylab(expression(P )) + xlab(expression(xi)) +
  theme_custom() + theme(axis.ticks.x=element_blank()) +
  annotate("segment", x=10^-8, xend=10^-2, y=-(3/pi)*0.75, yend=(3/pi)*(-0.75 + 6*0.0987))+
  annotate("segment", x=10^-8, xend=10^-4, y=-(3/pi)*19, yend=(3/pi)*(-19 + 4*3.289))+
  theme(plot.margin = unit(c(0.11,0.43,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(legend.position="none")
  
  
  
#multiplot(p4,p5,p6, cols=3)
p <- qplot(-20:-19,1:2) +
     coord_cartesian(xlim=c(0,600),ylim=c(0,400),expand=FALSE) +
     theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
#p = qplot(1:10, 1:10) + theme_bw()

g1 = plot_custom(p1)
g2 = plot_custom(p2)
g3 = plot_custom(p3)
g4 = ggplotGrob(p4)
g5 = ggplotGrob(p5)
g6 = ggplotGrob(p6)

pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/hprime-p-x-full.pdf", useDingbats = FALSE, width = 6, height = 3.8)
p + annotation_custom(grob = g1, xmin = 0, xmax = 200, ymin = 200, ymax = 400) +
annotation_custom(grob = g2, xmin = 200,  xmax = 400,   ymin = 200, ymax = 400) +
annotation_custom(grob = g3, xmin = 390,  xmax = 600,   ymin = 200, ymax = 400) + 
annotation_custom(grob = g4, xmin = -8.5, xmax = 196.5, ymin = 0, ymax = 200) +
annotation_custom(grob = g5, xmin = 186.5,xmax = 396.5, ymin = 0, ymax = 200) +
annotation_custom(grob = g6, xmin = 400,  xmax = 600,   ymin = 0, ymax = 200) 


dev.off()