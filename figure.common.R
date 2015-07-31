#' Common variables for plots
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(PEcAnRTM)
library(reshape2)
load("../data/FFT.processed.RData")

succ.colors <- scale_color_manual(values = c("green3", "blue", "red"))

png.plot <- function(fname, h=4, w=4, ...){
    png(fname, height=h, width=w, res=300, units="in", ...)
}

get_legend <- function(g){
  tmp <- ggplot_gtable(ggplot_build(g))
  leg <- which(sapply(tmp$grobs, '[[', 'name') == 'guide-box')
  leg2 <- tmp$grobs[[leg]]
  return(leg2)
}
