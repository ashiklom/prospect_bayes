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
y1 <- 21
y2 <- 23
Car.plot <- gen.plot + aes(x=Car, y=Car.mu) + ylab("Car") + no.x + ylim(0,y2) +
    geom_segment(aes(x=Car.out.x, y=y1, xend=Car.out.x, yend=y2), size=0.5, color="purple",
                 arrow=arrow(length=unit(0.02, "in"), type="closed"))
Cw.plot <- gen.plot + aes(x=Cw, y=Cw.mu) + ylab("Cw") + no.x
Cm.plot <- gen.plot + aes(x=Cm, y=Cm.mu) + ylab("Cm") + no.x

pdf("manuscript/figures/sensor-bias.pdf", height=7, width=7)
grid.arrange(N.plot, Cab.plot, Car.plot, Cw.plot, Cm.plot, ncol=1)
dev.off()

#' Here, we explore a few alternative ways of plotting the results of the 
#' sensor simulation experiment. The figure from this section that ended up in 
#' the supplementary information shows the relative width of the 95% confidence 
#' interval of the inversion parameters relative to the mean as a function of 
#' the true parameter value.  To facilitate generation and customization of 
#' this figure, we take advantage of the `sprintf` function for string 
#' processing combined with the `aes_string` feature from ggplot to create a 
#' common plotting function that we then `lapply` over the list of parameters 
#' and use `do.call` to arrange the resulting list of plots.

gp2 <- ggplot(simulation.dat) + facet_grid(.~sensor) + geom_point(size=1) + gen.theme
plt.error <- function(param, string)
    gp2 + aes_string(x=param, y=sprintf(string, param)) + no.x

# Plot of relative error vs.
#cv.str <- "%1$s.mu/%1$s - 1"
#cv.list <- c(lapply(params.prospect5, plt.error, string=cv.str), ncol=1)
#do.call(grid.arrange, cv.list)

riqr.str <- "(%1$s.q975 - %1$s.q25)/%1$s.mu"
riqr.list <- c(lapply(params.prospect5, plt.error, string=riqr.str), ncol=1)
pdf("manuscript/figures/sensor-riqr.pdf", height=7, width=7)
do.call(grid.arrange, riqr.list)
dev.off()

