# Remember to run multiplot before using this script (should be in /Graphs)
pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-0.pdf", useDingbats = FALSE, width = 6, height = 2)

dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/KI-0.csv")

library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p1 <- ggplot(data=dat, aes(x=l0, y= KII2))  +
      geom_line(linetype="solid", size=0.2) +
      geom_point(size=0.75) +
      scale_colour_manual(values=cbPalette1,name  = NULL) +
     # scale_shape_manual(values = c(15,16,17,18), name  = NULL) +
      coord_cartesian(ylim=c(0,2),xlim=c(0.057,0.101),expand=FALSE) +
      ylab(expression(K[II]^2))+ xlab(expression(lambda[0]))+
      theme_bw() +
      theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
      theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6))


p2 <- ggplot(data=dat, aes(x=L, y= KII0))  +
      geom_line(linetype="solid", size=0.2) +
      geom_point(size=1) +
      scale_colour_manual(values=cbPalette1,name  = NULL) +
     # scale_shape_manual(values = c(15,16,17,18), name  = NULL) +
     # coord_cartesian(ylim=c(0.1,0.4),xlim=c(0.8,3.4),expand=FALSE) +
     ylab(expression(K[II]))+
      theme_bw() +
      theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
      theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6))
p3 <- ggplot(data=dat, aes(x=L, y= l0))  +
      geom_line(linetype="solid", size=0.2) +
      geom_point(size=1) +
      scale_colour_manual(values=cbPalette1,name  = NULL) +
     # scale_shape_manual(values = c(15,16,17,18), name  = NULL) +
     # coord_cartesian(ylim=c(0.1,0.4),xlim=c(0.8,3.4),expand=FALSE) +
     ylab(expression(lambda[0]))+ 
      theme_bw() +
      theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
      theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6))
   
 multiplot(p1,p2,p3, cols=3)
 
 dev.off()
