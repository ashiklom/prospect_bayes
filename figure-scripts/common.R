#' Common variables for plots
library(ggplot2)
library(gridExtra)

succ.colors <- scale_color_manual(values = c("green3", "blue", "red"))

png.plot <- function(fname, h=4, w=4, ...){
    png(fname, height=h, width=w, res=300, units="in", ...)
}
