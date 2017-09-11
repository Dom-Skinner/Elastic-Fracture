
dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-0.csv")

library(ggplot2)
library(scales)    
library(grid) 
library(gtable)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p1 <- ggplot(data=dat, aes(x=12*(3/pi)*l0, y= KII2))  +
      geom_line(linetype="solid", size=0.3) +
      geom_point(size=1,shape =ifelse(dat$KII0 > 0,15,0)) +
      scale_colour_manual(values=cbPalette1,name  = NULL) +
      coord_cartesian(ylim=c(0,2),xlim=c(0.6,1.201),expand=FALSE) +
      ylab(expression(k[II]^2))+ xlab(expression(V[0]))+
      theme_custom() +
      #theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
      theme(axis.title = element_text(size=8,family="Palatino"))

dat$bool <- dat$KII0 >=  0

p2 <- ggplot(data=dat, aes(x=L, y= KII0,group=bool,shape=bool))  + 
      geom_line( size=c(.0,.0,.0,0.0,0.3,0.3,0.3,0.3,0.3,.3,.3,.3,.3,.3,0.3,0.3,0.3,0.3,0.3)) +
      geom_point(size=1,shape = c(15,15,15,15,15,15,15,15,15,15,15,15,15,15,0,0,0,0,0)) +
      # Apologies for the lazyness...
      scale_colour_manual(values=cbPalette1,name  = NULL) +
      coord_cartesian(ylim=c(-0.25,1.5),xlim=c(0,4),expand=FALSE) +
     ylab(expression(k[II]))+
      theme_custom() +
      theme(axis.title = element_text(size=8,family="Palatino"))
      
      
p3 <- ggplot(data=dat, aes(x=L, y= 12*(3/pi)*l0,group=bool,shape=bool))  +
 geom_line( size=c(.0,.0,.0,0.0,0.3,0.3,0.3,0.3,0.3,.3,.3,.3,.3,.3,0.3,0.3,0.3,0.3,0.3)) +
      geom_point(size=1,shape = c(15,15,15,15,15,15,15,15,15,15,15,15,15,15,0,0,0,0,0)) +
      scale_colour_manual(values=cbPalette1,name  = NULL) +
      coord_cartesian(ylim=c(12*(3/pi)*0.055,12*(3/pi)*0.105),xlim=c(0,4),expand=FALSE) +
      ylab(expression(V[0]))+ 
      theme_custom()+
      #theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
      theme(axis.title = element_text(size=8,family="Palatino"))
   

g1 = plot_custom(p1)
g2 = plot_custom(p2) 
g3 = plot_custom(p3) 
    
p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() +
     coord_cartesian(xlim=c(0,540),ylim=c(0,180),expand=FALSE) +
     theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())

pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-0.pdf", useDingbats = FALSE, width = 5.4, height = 1.8)
p + annotation_custom(grob = g1, xmin = -10, xmax = 180, ymin = -20, ymax = 180) +
annotation_custom(grob = g2, xmin = 175, xmax = 365, ymin = -17, ymax = 180) +
annotation_custom(grob = g3, xmin = 360, xmax = 550, ymin = -17, ymax = 180)
 dev.off()
