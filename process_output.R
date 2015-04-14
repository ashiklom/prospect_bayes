## Script to process FFT individual leaf results
library(data.table)
load("run_results/FFTdata.Rdata")

## Basic statistics on individuals
leaf.means <- all.r[,lapply(.SD, mean), by = Spectra]
leaf.sd <- all.r[,lapply(.SD, sd), by=Spectra]

## Log statistics on individuals
all.r.log <- data.table(Spectra = all.r$Spectra)
for(cn in names(all.r)[2:5]) all.r.log[, cn := log(all.r[,cn, with=F]), with=F]
leaf.logmean <- all.r.log[,lapply(.SD, mean), by=Spectra]
leaf.logsd <- all.r.log[,lapply(.SD, mean), by=Spectra]

rename <- function(dt, s)
        setnames(dt, names(dt)[-1], sprintf("%s.%s", names(dt)[-1], s))

rename(leaf.means, "mu")
rename(leaf.sd, "sd")
rename(leaf.logmean, "lmu")
rename(leaf.logsd, "lsd")

## Merge data tables
setkey(leaf.means, Spectra)
setkey(leaf.sd, Spectra)
setkey(leaf.logmean, Spectra)
setkey(leaf.logsd, Spectra)
fft.data <- leaf.means[leaf.sd][leaf.logmean][leaf.logsd]
