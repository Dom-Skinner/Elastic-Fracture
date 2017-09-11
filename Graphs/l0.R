pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/l0.pdf", useDingbats = FALSE, width = 3, height = 2)
dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/l0.csv")
dat$type <- factor(dat$type, levels=c("Kval","linear","quadratic"),
  labels=c("Numerical values","Linear extrapolation","Quadratic extrapolation"))


library(ggplot2)
library(grid) 
library(gtable)

cbPalette <- c("#000000","#e41a1c","#377eb8","#984ea3","#4daf4a")

p <-ggplot(data=dat, aes(x=K, y= 12*(3/pi)*lambda,group=type,shape=type,color=type,linetype=type))  +
  geom_line( size=0.4) +
  geom_point(size=0.9) +
  scale_colour_manual(values=cbPalette,name  = NULL) +
  scale_shape_manual(values = c(15,26,26), name  = NULL) +
  # Note that shape size 26 does not actually exist, this is just a hack
  scale_linetype_manual(values=c(0,1,2), name  = NULL) +
  coord_cartesian(ylim=c(0.65,0.68),xlim=c(0,0.275),expand=FALSE) +
  ylab(expression(V)) + xlab(expression(k[I]^u)) +
 # scale_y_continuous(breaks= seq(0.57,0.59,by=0.005), 
  #                   labels = insert_minor( seq(0.57, 0.59, by=0.005), 0 ))+ 
  scale_x_continuous(breaks= seq(0,0.25,by=0.05), 
                     labels = insert_minor( seq(0, 0.25, by=0.05), 0 ))+ 
  theme_custom()+
  #theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=8,family="Palatino"), legend.justification=c(1,1), legend.position="none")
  
grid.draw(plot_custom(p))

dev.off()