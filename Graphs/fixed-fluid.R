pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/fixed-fluid.pdf", useDingbats = FALSE, width = 3.5, height = 1.5)

dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/fixed-fluid.csv")


library(ggplot2)
library(scales)    
library(grid) 
library(gtable)

cbPalette <- c("#000000","#e41a1c","#377eb8","#984ea3","#ff7f00")

lambda <- dat$lambda
L08    <- dat$L08
L10    <- dat$L10
L22    <- dat$L22
KII    <- dat$KII0

df <- data.frame(lambda, L08,L10,L22,KII)
q <- ggplot(df, aes(x=12*(3/pi)*lambda, y = value, color = variable)) + 
    geom_line(aes(y = KII, col = "Calculated values"), size=0.5) + 
    geom_point(aes(y = L08, col = "L=0.8"), size=0.8) +
    geom_point(aes(y = L10, col = "L=1.0"), size=0.8, shape=17) +
    geom_point(aes(y = L22, col = "L=2.2"), size=0.8, shape = 15) +
    scale_colour_manual(values=cbPalette,name  = NULL) +
    scale_x_continuous(breaks = seq(1.02,1.16, by=0.02))+
    coord_cartesian(ylim=c(0,0.8),xlim=c(1.03,1.16),expand=FALSE) +
    theme_custom() + 
    ylab(expression(K[II]))+xlab(expression(V))+
    theme(plot.margin = unit(c(0.1,0.2,0,0), "cm")) + # ("top", "right", "bottom", "left")
    guides(col=guide_legend(ncol=2))+
    theme(axis.title = element_text(size=8,family="Palatino"),
    legend.justification=c(0,0), legend.position=c(-0.02,-0.07),legend.key = element_rect(color=NA),
    legend.text = element_text(size=8,family="Palatino"))
   
   grid.draw(plot_custom(q))

   dev.off()