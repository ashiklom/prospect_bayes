if(grepl("figure-scripts", getwd())) setwd("..")

library(data.table)
library(PEcAnRTM)
library(coda)
library(ggplot2)
library(GGally)
library(gridExtra)
library(gtable)
library(data.table)
library(reshape2)
library(MASS)

## Set paths
path.db <- file.path("~", "Documents", "Dropbox")
path.fft <- file.path(path.db, "FFT_spectra", "FFT_spectra_unique.csv")

## Read FFT spectra
fft.spec <- fread(path.fft, header=TRUE, drop=c(10,11,14:21,22:71))

## Sample inversion
obs <- fft.spec[Species == "ACRU", 12:2112, with=FALSE]
obs <- t(as.matrix(obs))[,1]
load("data/p4.defaults.RData")
#ss <- invert.fast("prospect_4", obs, inits, NULL, pmu, psd, plog, minp, ngibbs=20000)
#save(ss, file="data/testinv.RData")
load("data/testinv.RData")
par(mfrow=c(2,2))
for(i in 1:4) plot(ss[,i], type='l')

bt <- seq(10000, 20000)
samps <- ss[bt,]
sdat <- data.table(samps)
sdat[, Cw := Cw * 1000][, Cm := Cm * 1000]

#' Density plots
densfunc <- function(x,y,...){
    require(RColorBrewer)
    z <- kde2d(x, y, n=50)
    nlev <- 6
    my.cols <- rev(brewer.pal(nlev, "RdYlBu"))
    points(x, y, col="darkgrey", pch=".")
    contour(z, drawlabels=FALSE, nlevels=nlev, col="black", add=TRUE)
}
png.plot("manuscript/figures/parameter-covariance.png", h=4, w=4)
pairs(sdat[,1:4,with=F], 
      upper.panel=NULL, 
      lower.panel=densfunc, 
      labels = c("N", "Cab", "Cw\n(x1000)", "Cm\n(x1000)"),
      cex.axis = 0.8,
      mar=c(0, 0, 0, 0))
dev.off()
