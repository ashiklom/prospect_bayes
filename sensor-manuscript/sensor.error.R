# Big plot of sensor error
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(PEcAnRTM)
load("../data/simulation.samp.dat.RData")

simulation.dat[, sensor := factor(sensor.sensor, levels=sensor.list)]

gen.theme <- theme_bw() + 
    theme(axis.text = element_text(size=5),
          strip.text = element_text(size=8))
no.x <- theme(axis.title.x = element_blank())
gen.plot <- ggplot(simulation.dat) + 
    facet_grid(.~sensor) + geom_point(size=1) +
    geom_abline(linetype="dashed", color="red") +
    gen.theme
N.plot <- gen.plot + aes(x=true.param.N.true, y=N.mu) + ylab("N") + no.x
Cab.plot <- gen.plot + aes(x=true.param.Cab.true, y=Cab.mu) + ylab("Cab") + no.x
Car.plot <- gen.plot + aes(x=true.param.Car.true, y=Car.mu) + ylab("Car") + no.x +
    ylim(0,max(simulation.dat[,true.param.Car.true]))
Cw.plot <- gen.plot + aes(x=true.param.Cw.true, y=Cw.mu) + ylab("Cw") + no.x
Cm.plot <- gen.plot + aes(x=true.param.Cm.true, y=Cm.mu) + ylab("Cm") + no.x

png("manuscript/figures/sensor-error.png", height=7, width=7,
    units="in", res=300)
grid.arrange(N.plot, Cab.plot, Car.plot, Cw.plot, Cm.plot, ncol=1)
dev.off()
