# Sensor covariance plots
library(ggplot2)
library(PEcAnRTM)
library(gridExtra)
library(grid)
library(reshape2)

datnum <- 100
parnames <- c("N", "Cab", "Car", "Cw", "Cm")
for(s in sensor.list){
    f.path <- sprintf("../data/some-simulations/%s.%d.RData", s, datnum)
    load(f.path)
    f.name <- gsub(".*/(.*)[.]RData", "\\1", f.path)
    f.list <- get(f.name)
    rm(list=f.name)
    sensor <- list(sensor = f.list$sensor)
    samples <- f.list$samples
    true.param <- as.list(f.list$true.param)
    names(true.param) <- sprintf("%s.true", parnames[1:5])

# Burn in and thin
    ngibbs <- nrow(samples)
    burnin <- floor(ngibbs/2)
    thin <- (ngibbs - burnin) / 5000
    bt <- seq(burnin, ngibbs, by=thin)
    samples.sub <- samples[bt,]
    samples.pairs <- samples.sub[,1:5]
    colnames(samples.pairs) <- parnames[1:5]

# Pairs plot
    plotname <- sprintf("manuscript/figures/pairs-%s.png", s)
    png(plotname, height=6, width=6, units="in", res=300)
    pairs(samples.pairs, pch=".", main=s)
    dev.off()
}


