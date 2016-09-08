pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/overall-fit.pdf", useDingbats = FALSE, width = 5.4, height = 3.3)

dat1 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/overall-fit-1.csv")

library(ggplot2)
cbPalette1 <- c("#999999", "#E69F00")
cbPalette2 <- c("#56B4E9", "#009E73")
cbPalette3 <- c("#F0E442", "#0072B2")
cbPalette4 <- c("#D55E00", "#CC79A7")


KI        <-dat1$KI
KII        <-dat1$KII
lambda    <- dat1$lambda
estimated <- dat1$estimated
estimated2 <- dat1$estimated2

df <- data.frame(KI, lambda, estimated)
p1 <- ggplot(df, aes(x=KI, y = value, color = variable)) + 
    geom_line(aes(y = estimated, col = "Best fit"), size=0.3) + 
    geom_point(aes(y = lambda, col = "lambda"), size=1) +
    scale_colour_manual(values=cbPalette1,name  = NULL) +
    coord_cartesian(ylim=c(0,0.061),xlim=c(0.3,2),expand=FALSE) +
    theme_bw() + 
    ylab(expression(lambda))+xlab(expression(kappa[I]))+
    theme(plot.margin = unit(c(0.1,0.2,0,0), "cm")) + # ("top", "right", "bottom", "left")
   theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
   legend.justification=c(0,0), legend.position=c(0,0),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
   legend.text = element_text(size=8,family="Palatino"),legend.direction = "vertical")
   
df2 <- data.frame(KII, lambda, estimated2)
p2 <- ggplot(df2, aes(x=KII, y = value, color = variable)) + 
    geom_line(aes(y = estimated2, col = "Best fit"), size=0.3) + 
    geom_point(aes(y = lambda, col = "lambda"), size=1.5,shape=18) +
    scale_colour_manual(values=cbPalette2,name  = NULL) +
    coord_cartesian(ylim=c(0,0.061),xlim=c(1.38,1.51),expand=FALSE) +
    theme_bw() + 
    ylab(expression(lambda))+xlab(expression(kappa[II]))+
    theme(plot.margin = unit(c(0.1,0.2,0,0), "cm")) + # ("top", "right", "bottom", "left")
   theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
   legend.justification=c(0,0), legend.position=c(0,0),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
   legend.text = element_text(size=8,family="Palatino"),legend.direction = "vertical")   

   
dat2 <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/overall-fit-2.csv")
KII              <-dat2$KII
lam          <- dat2$lambda
lambdaEst <- dat2$lambdaEst
L                <- dat2$L
LEst      <- dat2$LEst


df <- data.frame(KII, lam, lambdaEst)
p3 <- ggplot(df, aes(x=KII, y = value, color = variable)) + 
    geom_line(aes(y = lambdaEst, col = "Best fit"), size=0.3) + 
    geom_point(aes(y = lam, col = "lambda"), size=1,shape = c(15)) +
    scale_colour_manual(values=cbPalette3,name  = NULL) +
    coord_cartesian(ylim=c(0.058,0.105),xlim=c(0,1.6),expand=FALSE) +
    theme_bw() + 
    ylab(expression(lambda))+xlab(expression(kappa[II]))+
    theme(plot.margin = unit(c(0.1,0.2,0,0), "cm")) + # ("top", "right", "bottom", "left")
   theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
   legend.justification=c(0,0), legend.position=c(0,0),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
   legend.text = element_text(size=8,family="Palatino"),legend.direction = "vertical")
   
df2 <- data.frame(KII, L, LEst)
p4 <- ggplot(df2, aes(x=KII, y = value, color = variable)) + 
    geom_line(aes(y = LEst, col = "Best fit"), size=0.3) + 
    geom_point(aes(y = L, col = "L"), size=1, shape=17) +
    scale_colour_manual(values=cbPalette4,name  = NULL) +
    coord_cartesian(ylim=c(0,3.5),xlim=c(-0.05,1.53),expand=FALSE) +
    theme_bw() + 
    ylab(expression(L))+xlab(expression(kappa[II]))+
    theme(plot.margin = unit(c(0.1,0.2,0,0), "cm")) + # ("top", "right", "bottom", "left")
   theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
   legend.justification=c(1,1), legend.position=c(1,1),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
   legend.text = element_text(size=8,family="Palatino"),legend.direction = "vertical")   



  multiplot(p1,p2,p3,p4, cols=2)
  dev.off()

