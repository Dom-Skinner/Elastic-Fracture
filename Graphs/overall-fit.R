#pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/overall-fit.pdf", useDingbats = FALSE, width = 5.6, height = 1.8)

dat1 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/overall-fit-1.csv")

library(ggplot2)
library(scales)    
library(grid) 
library(gtable)

cbPalette1 <- c("#000000","#000000")
cbPalette2 <- c("#000000","#000000")
cbPalette3 <- c("#000000","#000000")


KI        <-dat1$KI
KII        <-dat1$KII
lambda    <- dat1$lambda
estimated <- dat1$estimated
estimated2 <- dat1$estimated2

df <- data.frame(KI, lambda, estimated)
p1 <- ggplot(df, aes(x=KI, y = value, color = variable)) + 
    geom_line(aes(y = estimated, col = "Best fit"), size=0.3) + 
    geom_point(aes(y = lambda, col = "lambda"), size=1.3,shape=15) +
    scale_colour_manual(values=cbPalette1,name  = NULL) +
    coord_cartesian(ylim=c(0,0.061),xlim=c(0.3,2),expand=FALSE) +
    theme_custom() + 
    ylab(expression(lambda))+xlab(expression(kappa[I]))+
    theme(plot.margin = unit(c(0.1,0.2,0,0), "cm")) + # ("top", "right", "bottom", "left")
   theme(axis.title = element_text(size=8,family="Palatino"), legend.position="none")
   
df2 <- data.frame(KII, lambda, estimated2)
p2 <- ggplot(df2, aes(x=KII, y = value, color = variable)) + 
    geom_line(aes(y = estimated2, col = "Best fit"), size=0.3) + 
    geom_point(aes(y = lambda, col = "lambda"), size=1.3,shape=15) +
    scale_colour_manual(values=cbPalette2,name  = NULL) +
    coord_cartesian(ylim=c(0,0.061),xlim=c(1.38,1.51),expand=FALSE) +
    theme_custom() + 
    ylab(expression(lambda))+xlab(expression(kappa[II]))+
    theme(plot.margin = unit(c(0.1,0.2,0,0), "cm")) + # ("top", "right", "bottom", "left")
   theme(axis.title = element_text(size=8,family="Palatino"), legend.position="none")
   

   
dat2 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/overall-fit-2.csv")
KII              <-dat2$KII
lam          <- dat2$lambda
lambdaEst <- dat2$lambdaEst
L                <- dat2$L
LEst      <- dat2$LEst


df <- data.frame(KII, lam, lambdaEst)
p3 <- ggplot(df, aes(x=KII, y = value, color = variable)) + 
    geom_line(aes(y = lambdaEst, col = "Best fit"), size=0.3) + 
    geom_point(aes(y = lam, col = "lambda"), size=1.3,shape = c(15)) +
    scale_colour_manual(values=cbPalette3,name  = NULL) +
    coord_cartesian(ylim=c(0.058,0.105),xlim=c(0,1.6),expand=FALSE) +
    theme_custom() + 
    ylab(expression(lambda))+xlab(expression(kappa[II]))+
    theme(plot.margin = unit(c(0.1,0.2,0,0), "cm")) + # ("top", "right", "bottom", "left")
   theme(axis.title = element_text(size=8,family="Palatino"), legend.position="none")
   
#df2 <- data.frame(KII, L, LEst)
#p4 <- ggplot(df2, aes(x=KII, y = value, color = variable)) + 
#    geom_line(aes(y = LEst, col = "Best fit"), size=0.3) + 
#    geom_point(aes(y = L, col = "L"), size=1, shape=17) +
#    scale_colour_manual(values=cbPalette4,name  = NULL) +
#    coord_cartesian(ylim=c(0,3.5),xlim=c(-0.05,1.53),expand=FALSE) +
#    theme_bw() + 
#    ylab(expression(L))+xlab(expression(kappa[II]))+
#    theme(plot.margin = unit(c(0.1,0.2,0,0), "cm")) + # ("top", "right", "bottom", "left")
#   theme(axis.title = element_text(size=8,family="Palatino"), legend.position="none")

g1 = plot_custom(p1)
g2 = plot_custom(p2)
g3 = plot_custom(p3)

p <- qplot(-20:-19,1:2) +
     coord_cartesian(xlim=c(0,570),ylim=c(0,180),expand=FALSE) +
     theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())

pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/overall-fit.pdf", useDingbats = FALSE, width = 5.5, height = 1.8)
p + annotation_custom(grob = g1, xmin = 0   , xmax = 190, ymin = 0, ymax = 180) +
    annotation_custom(grob = g2, xmin = 190, xmax = 380, ymin = 0, ymax = 180) +
    annotation_custom(grob = g3, xmin = 380 , xmax = 570, ymin = 0, ymax = 180) 

dev.off()

