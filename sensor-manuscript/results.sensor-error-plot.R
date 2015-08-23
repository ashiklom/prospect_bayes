# Big plot of sensor error
library(PEcAnRTM)
library(ggplot2)
library(gridExtra)
library(grid)
load("../data/simulation.samp.dat.RData")

simulation.dat[, sensor := factor(sensor, levels=sensor.list)]

gen.theme <- theme_bw() + 
    theme(axis.text = element_text(size=5),
          strip.text = element_text(size=8))
no.x <- theme(axis.title.x = element_blank())
gen.plot <- ggplot(simulation.dat) + 
    facet_grid(.~sensor) + geom_point(size=1) +
    geom_abline(linetype="dashed", color="red") +
    gen.theme
N.plot <- gen.plot + aes(x=N, y=N.mu) + ylab("N") + no.x
Cab.plot <- gen.plot + aes(x=Cab, y=Cab.mu) + ylab("Cab") + no.x
Car.plot <- gen.plot + aes(x=Car, y=Car.mu) + ylab("Car") + no.x
Cw.plot <- gen.plot + aes(x=Cw, y=Cw.mu) + ylab("Cw") + no.x
Cm.plot <- gen.plot + aes(x=Cm, y=Cm.mu) + ylab("Cm") + no.x

pdf("manuscript/figures/sensor-error.pdf", height=7, width=7)
grid.arrange(N.plot, Cab.plot, Car.plot, Cw.plot, Cm.plot, ncol=1)
dev.off()
