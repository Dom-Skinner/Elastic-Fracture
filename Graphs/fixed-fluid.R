pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/fixed-fluid.pdf", useDingbats = FALSE, width = 4.5, height = 2)

dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/fixed-fluid.csv")

library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

lambda <- dat$lambda
L08    <- dat$L08
L10    <- dat$L10
L22    <- dat$L22
KII    <- dat$KII0

df <- data.frame(lambda, L08,L10,L22,KII)
ggplot(df, aes(x=lambda, y = value, color = variable)) + 
    geom_line(aes(y = KII, col = "Calculated values"), size=0.5) + 
    geom_point(aes(y = L08, col = "L=0.8"), size=0.8) +
    geom_point(aes(y = L10, col = "L=1.0"), size=0.8, shape=17) +
    geom_point(aes(y = L22, col = "L=2.2"), size=0.8, shape = 15) +
    scale_colour_manual(values=cbPalette,name  = NULL) +
#    coord_cartesian(ylim=c(0,0.061),xlim=c(0.3,2)) +
    theme_bw() + 
    ylab(expression(K[II]))+xlab(expression(lambda))+
    theme(plot.margin = unit(c(0.1,0.2,0,0), "cm")) + # ("top", "right", "bottom", "left")
   theme(axis.title = element_text(size=8,family="Palatino"), axis.text  = element_text(size=6),
   legend.justification=c(0,0), legend.position=c(0,0),legend.background = element_rect(color="#808080",fill="#F0F0F0", size=0.25),
   legend.text = element_text(size=8,family="Palatino"),legend.direction = "vertical")
   
   dev.off()