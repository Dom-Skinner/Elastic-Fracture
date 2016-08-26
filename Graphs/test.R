pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/K-lambda.pdf", useDingbats = FALSE, width = 4, height = 2)
dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/K-lambda.csv")

library(ggplot2)

ggplot(data=dat, aes(x=lambda, y=K)) + 
  geom_line(colour="black", linetype="solid", size=0.3) + 
  geom_point(colour="black", size=0.8, shape=21, fill="black") +
  coord_cartesian(ylim=c(0,2),xlim=c(0,0.06),expand=FALSE) +
  ylab(expression(kappa)) + xlab(expression(lambda)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.3,0,0), "cm")) + # ("top", "right", "bottom", "left")
  theme(axis.title = element_text(size=10,family="Palatino"), axis.text  = element_text(size=6))

dev.off()
