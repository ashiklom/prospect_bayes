# Generate scripts for inversions on known parameters
library(data.table)
library(PEcAnRTM)
load("../data/FFT.processed.RData")

# Get quantiles of FFT inversion
params <- c("N", "Cab", "Car", "Cw", "Cm")
params.mu <- paste0(params, ".mu")
quants <- fft.f[, lapply(.SD, quantile,
                         c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975),
                         na.rm=TRUE), .SDcols=params.mu]
qmat <- as.matrix(quants)

# Generate full parameter inversion matrix
full.mat.list <- list()
for(p in 1:5){
    parmat <- matrix(NA, 7*7, 5)
    for(h in 1:7){
        for(a in 1:7){
            matrow <- c(qmat[h,p], qmat[a,-p])
            matrow <- matrow[params.mu]
            parmat[a+7*(h-1),] <- matrow
        }
    }
    full.mat.list[[p]] <- parmat
}
full.mat <- do.call(rbind, full.mat.list)

# Write submission script
n.sensor <- length(sensor.list)
N <- rep(full.mat[,1], n.sensor)
Cab <- rep(full.mat[,2], n.sensor)
Car <- rep(full.mat[,3], n.sensor)
Cw <- rep(full.mat[,4], n.sensor)
Cm <- rep(full.mat[,5], n.sensor)
sensor <- rep(sensor.list, each=nrow(full.mat))
ngibbs <- 50000
run.param <- rep(params, each=7*7)
run.pnum <- rep(1:7, each=7)
run.qnum <- rep(1:7, 7)
run.par <- sprintf("%s%d.q%d", run.param, rep(run.pnum, 5), rep(run.qnum, 5))
runname <- sprintf("%s.%s", rep(run.par, n.sensor), sensor)
submit.string <- sprintf('qsub -V -v N=%f,Cab=%f,Car=%f,Cw=%f,Cm=%f,sensor=%s,ngibbs=%d,runname=%s -N "%s" submit.simulated.qsub',
                         N, Cab, Car, Cw, Cm, sensor, ngibbs, runname, runname)

# Create file
fname <- "run-scripts/run.simulation.sh"
write("#!/bin/bash", file=fname)
write(submit.string, file=fname, append=TRUE)
