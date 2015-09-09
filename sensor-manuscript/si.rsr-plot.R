#' ---
#' title: Spectral response curve plot
#' author: Alexey Shiklomanov
#' output_format: html_document
#' ---

#' # Introduction
#' This script is used to generate a figure showing the band widths and 
#' locations of the sensors used in this manuscript.

#' # Setup

#' First, we load dependencies. The `PEcAnRTM` package defines the sensor 
#' response curves. Furthermore, we have to load the `sensor.rsr` data object, 
#' which contains the RSR data.

library(PEcAnRTM)
data(sensor.rsr)

#' We then remove two of the hyperspectral sensors, given how similar they are 
#' and how difficult their plots are to interpret.

sensors.remove <- c("aviris.ng", "hyperion")
sensors.to.plot <- sensor.rsr[!names(sensor.rsr) %in% sensors.remove]

#' # Plot
#' This plot is generated using R base graphics, specifically using the 
#' `matplot` function, which makes plotting data organized in matrices very 
#' simple. First, we set up the plot layout (`mfrow`; 4 rows, 2 columns) and 
#' margins (`mar`; bottom, left, top, right). Then, we loop over the 
#' `sensors.to.plot` list. For each plot, we convert the index (1-2101) to a 
#' wavelength (400-2500) by adding 399, extract the proper name, and draw the 
#' line plot.

par(mfrow = c(4,2), mar=c(2.5, 2.5, 2, 1))
for(s.name in names(sensors.to.plot)){
    s <- sensors.to.plot[[s.name]]
    wl <- s[,"index"] + 399
    s.proper <- sensor.proper[s.name]
    matplot(wl, s[,-1], type='l', lty=1, xlim=c(400,2500), 
            xlab="", ylab="", main=s.proper)
}
