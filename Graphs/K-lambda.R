pdf(file = "./Documents/Summer-Project/Elastic-Fracture/Graphs/K-lambda.pdf", useDingbats = FALSE, width = 3, height = 2)
dat <- read.csv("./Documents/Summer-Project/Elastic-Fracture/Graphs/K-lambda.csv")

library(ggplot2)
library(grid) 
library(gtable)

insert_minor <- function(major_labs, n_minor) {labs <- 
                              c( sapply( major_labs, function(x) c(x, rep("", n_minor) ) ) )
                              labs[1:(length(labs)-n_minor)]}

p <- ggplot(data=dat, aes(x=12*(3/pi)*lambda, y=K)) + 
  geom_line(colour="#000000", linetype="solid", size=0.25) + 
  geom_point(colour="#000000", size=1.2, shape=15, fill="#000000") +
  coord_cartesian(ylim=c(0,2),xlim=c(0,0.7),expand=FALSE) +
  ylab(expression(k)) + xlab(expression(V)) +
  scale_y_continuous(breaks= seq(0,2,by=0.5), 
                     labels = insert_minor( seq(0, 2, by=0.5), 0 ))+ 
  scale_x_continuous(breaks= seq(0,0.7,by=0.1), 
                     labels = insert_minor( seq(0, 0.7, by=0.1), 0 ))+ 
  theme_custom() +
  theme(plot.margin = unit(c(0.2,0.3,0,0), "cm"))  # ("top", "right", "bottom", "left")
  

grid.draw(plot_custom(p))

dev.off()
