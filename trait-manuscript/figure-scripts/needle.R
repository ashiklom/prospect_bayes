if(grepl("figure-scripts", getwd())) setwd("..")
source("preprocess.fft.R")
source("figure-scripts/common.R")
library(ggplot2)
library(gridExtra)

#' Compare new vs everything else
fft.c[Age == "N", new.needle := "New"]
fft.c[Age != "N", new.needle := "Old"]

#' Summary plot
ns.N <- ggplot(fft.c) + 
    aes(y=N.mu, x=new.needle) +
    geom_violin() +
    xlab("Needle age") +
    ylab("N") +
    theme_bw() +
    theme(text = element_text(size=11),
          axis.text = element_text(size=rel(0.6)),
          axis.title = element_text(size=rel(0.9)))
ns.Cab <- ns.N + aes(y=Cab.mu) + ylab("Cab")
ns.Cw <- ns.N + aes(y=Cw.mu) + ylab("Cw")
ns.Cm <- ns.N + aes(y=Cm.mu) + ylab("Cm")
png.plot("manuscript/figures/needle-summary.png", w=4, h=4)
grid.arrange(ns.N, ns.Cab, ns.Cw, ns.Cm, nrow=2)
dev.off()

#' LMA comparison
nplot <- ggplot(fft.c) + 
    aes(x=Cm.mu, y=LMA_g_DW_cm2, color=new.needle) +
    geom_point()
plot(nplot)
