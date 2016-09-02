pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/catagory.pdf", useDingbats = FALSE, width = 4, height = 2)

dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/catagory.csv")

library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")



ggplot(data=dat, aes(x=KI, y= KII))  +
      geom_line(linetype="solid", size=0.5) +
      geom_point(size=1.5) +
      geom_hline(aes(yintercept=dat$KII[length(dat$KII)]))+
      #geom_vline(aes(xintercept=dat$KI[1]))+
      geom_linerange(aes(x=dat$KI[1], y=NULL, ymin=dat$KII[1], ymax=4))+
      scale_colour_manual(values=cbPalette1,name  = NULL) +
     # scale_shape_manual(values = c(15,16,17,18), name  = NULL) +
      coord_cartesian(ylim=c(1.3,1.6),xlim=c(0,2.3),expand=FALSE) +
      ylab(expression(kappa[II]))+ xlab(expression(kappa[I]))+
      theme_bw() +
      theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
      theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6))
    

 dev.off()