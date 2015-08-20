# Spectra simulated using PROSPECT

# Generate scripts for inversions on known parameters
library(data.table)
library(PEcAnRTM)
load("../data/FFT.processed.RData")
vars <- c("N.mu", "Cab.mu", "Car.mu", "Cw.mu", "Cm.mu") 
fft.sub <- fft.f[!is.na(N.mu)][sensor=="identity"][,vars,with=F]
full.mat <- as.matrix(fft.sub)
full.mat <- full.mat[sample(1:nrow(full.mat)), ]

# Write submission script
n.sensor <- length(sensor.list)
N <- rep(full.mat[,1], each=n.sensor)
Cab <- rep(full.mat[,2], each=n.sensor)
Car <- rep(full.mat[,3], each=n.sensor)
Cw <- rep(full.mat[,4], each=n.sensor)
Cm <- rep(full.mat[,5], each=n.sensor)
sensor <- rep(sensor.list, nrow(full.mat))
ngibbs <- 50000
runname <- sprintf("%s.%d", sensor, rep(1:nrow(full.mat), each=n.sensor))
submit.string <- sprintf('qsub -V -v N=%f,Cab=%f,Car=%f,Cw=%f,Cm=%f,sensor=%s,ngibbs=%d,runname=%s -N "%s" submit.simulated.qsub',
                         N, Cab, Car, Cw, Cm, sensor, ngibbs, runname, runname)

# Create file
fname <- "run-scripts/run.simulation.sh"
write("#!/bin/bash", file=fname)
write(submit.string, file=fname, append=TRUE)
