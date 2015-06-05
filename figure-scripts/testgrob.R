library(gridExtra)

plots <- lapply(1:5, function(x) qplot(1:10, rnorm(10), main=paste("plot", x)))

grid.arrange(arrangeGrob(plots[[1]], plots[[2]], widths=c(3,2), nrow=1),
             arrangeGrob(plots[[3]], plots[[4]], plots[[5]], nrow=1))

