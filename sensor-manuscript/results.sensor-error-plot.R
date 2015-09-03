#' ---
#' title: Sensor error plot
#' author: Alexey Shiklomanov
#' output_format: html_document
#' ---

#' # Introduction
#' This script generates a large matrix of scatterplots that show the ability 
#' of each sensor to accurately retrieve parameters from the simulations. 

#' # Setup
#' First, we load dependencies and data. The PEcAnRTM package is used only for 
#' pre-defined sensor information. The remaining packages are used to generate 
#' the plots. For details on how the data file is generated, see the 
#' `load.sim.R` file. 

library(PEcAnRTM)
library(ggplot2)
library(gridExtra)
library(grid)
load("../data/simulation.samp.dat.RData")

#' To facilitate plot organization, we first change the sensor names to their 
#' pretty versions (`sensor.proper`) and then convert the representation of the 
#' sensor from character to factor with a pre-defined order based on spectral 
#' resolution.
simulation.dat[, sensor := sensor.proper[sensor]]
simulation.dat[, sensor := factor(sensor, levels=sensor.proper)]

#' We found that the inversion of Carotenoids had a number of substantial 
#' outliers that caused problems for the interpretation of this plot. To 
#' alleviate this problem, we extracted the outliers and their x values into a 
#' separate column to be displayed as annotated arrows pointing outside the 
#' plot extent and recoded the original values as `NA`. We arbitrarily define 
#' outliers as Car estimates that are more than double the true estimate.
simulation.dat[Car.mu/Car > 2, 
               c("Car.mu", "Car", "Car.out.x", "Car.out.y") := list(NA, NA, Car, Car.mu)]

#' # Plot specifications
#' The following lines provide graphical preferences for the plot. The `no.x` 
#' theme is an additional set of parameters that removes x axis labels. The 
#' `gen.plot` object is a generalized `ggplot` object describing the data's 
#' presentation, but without any data present (subsequent lines add this data 
#' in turn for each parameter). For each parameter (rows) and sensor (columns, 
#' in order of decreasing spectral resolution), this script draws one 
#' scatterplot of true parameter values (x) vs mean inversion estimates (y) 
#' with a dashed red one-to-one line. Specifically, each parameter (row) is 
#' created as a separate plot, and these plots are then stacked on top of each 
#' other using the `grid.arrange` function.

gen.theme <- theme_bw() + 
    theme(axis.text = element_text(size=5),
          strip.text = element_text(size=5))
no.x <- theme(axis.title.x = element_blank())
gen.plot <- ggplot(simulation.dat) + 
    facet_grid(.~sensor) + geom_point(size=1) +
    geom_abline(linetype="dashed", color="red") +
    gen.theme

N.plot <- gen.plot + aes(x=N, y=N.mu) + ylab("N") + no.x
Cab.plot <- gen.plot + aes(x=Cab, y=Cab.mu) + ylab("Cab") + no.x
y1 <- 22.5
y2 <- 23
Car.plot <- gen.plot + aes(x=Car, y=Car.mu) + ylab("Car") + no.x + ylim(0,y2) +
    geom_segment(aes(x=Car.out.x, y=y1, xend=Car.out.x, yend=y2), size=0.5, color="purple",
                 arrow=arrow(length=unit(0.02, "in"), type="closed"))
Cw.plot <- gen.plot + aes(x=Cw, y=Cw.mu) + ylab("Cw") + no.x
Cm.plot <- gen.plot + aes(x=Cm, y=Cm.mu) + ylab("Cm") + no.x

pdf("manuscript/figures/sensor-error.pdf", height=7, width=7)
grid.arrange(N.plot, Cab.plot, Car.plot, Cw.plot, Cm.plot, ncol=1)
dev.off()
