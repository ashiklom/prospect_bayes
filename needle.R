source("preprocess.fft.R")
library(ggplot2)
library(data.table)

nplot <- ggplot(fft.c) + 
    aes(x=Cm.mu, y=LMA_g_DW_cm2, color=Age) +
    geom_point()
plot(nplot)
