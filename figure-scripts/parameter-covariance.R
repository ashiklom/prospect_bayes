if(grepl("figure-scripts", getwd())) setwd("..")

library(data.table)
library(PEcAnRTM)
library(coda)
library(ggplot2)
library(gridExtra)
library(gtable)
library(data.table)
library(reshape2)

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

#' Density plots
th.dens <- theme_bw()
dvar <- ggplot(sdat) + geom_density() + th.dens
dplot <- ggplot(sdat) + 
    stat_density2d(aes(fill=..level..), geom="polygon", color="black") +
    scale_fill_continuous(low="white", high="red") +
    guides(fill=FALSE) +
    th.dens
#dplot <- ggplot(sdat) + geom_point(size=1)
lim.N <- sdat[, c(min(N), max(N))]
d1 <- dvar + aes(x=N) + xlim(lim.N)
d2 <- dvar + aes(x=Cab)
d3 <- dvar + aes(x=Cw)
d4 <- dvar + aes(x=Cm) + coord_flip()
d12 <- dplot + aes(x=N, y=Cab) + xlim(lim.N)
d13 <- dplot + aes(x=N, y=Cw)
d14 <- dplot + aes(x=N, y=Cm)
d23 <- dplot + aes(x=Cab, y=Cw)
d24 <- dplot + aes(x=Cab, y=Cm)
d34 <- dplot + aes(x=Cw, y=Cm)
blank2 <- grid.rect(gp=gpar(col="white"))
blank <- rectGrob()

d1g <- ggplot_gtable(ggplot_build(d1))
d2g <- ggplot_gtable(ggplot_build(d2))
d12g <- ggplot_gtable(ggplot_build(d12))
mat.lay <- matrix(list(d1g, blank, blank, blank,
             d12g, d2g, blank, blank),
             nrow=2, byrow=TRUE)
g <- gtable_matrix("try1", mat.lay, widths=rep(1,4), heights=rep(1,2))
plot(g)


d1g <- ggplotGrob(d1)
d12g <- ggplotGrob(d12)
g <- gtable:::rbind_gtable(d1g, blank, blank d12g, "first")
grid.newpage()
grid.draw(g)

grid.arrange(d1, blank, blank, blank,
             d12, d2, blank, blank,
             d13, d23, d3, blank,
             d14, d24, d34, d4,
             nrow=4)

