pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/catagory.pdf", useDingbats = FALSE, width = 3.2, height = 2)

dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/catagory.csv")

library(ggplot2)
library(scales)    
library(grid) 
library(gtable)

arrow_data=data.frame(x=c(0.2,0.6,1,1.4,2.25,2.25), y=c(1.58,1.58,1.58,1.58,1.465,1.415), vx=c(0,0,0,0,-0.4,-0.4), vy=c(-0.12,-0.12,-0.12,-0.12,0,0))


p<-ggplot(data=dat, aes(x=KI, y= KII))  +
      geom_line(linetype="solid", size=0.5) +
      geom_point(size=1.5,shape=15) +
      geom_hline(aes(yintercept=dat$KII[length(dat$KII)]))+
      annotate("segment", x=1.95, xend=3, y=1.505, yend=1.505)+
      geom_linerange(aes(x=dat$KI[1], y=NULL, ymin=dat$KII[1], ymax=4))+
      geom_segment(data=arrow_data, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), arrow=arrow(length=unit(1.2, "mm"),type = "closed"), size=0.5) +
      coord_cartesian(ylim=c(1.3,1.6),xlim=c(0,2.3),expand=FALSE) +
      ylab(expression(kappa[II]))+ xlab(expression(kappa[I]))+
      theme_custom() +
      theme(axis.title = element_text(size=8,family="Palatino"))
    
grid.draw(plot_custom(p))
 dev.off()