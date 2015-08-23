# Spectra simulated using PROSPECT

# Generate scripts for inversions on known parameters
library(PEcAnRTM)
load("../data/FFT.processed.RData")
vars <- sprintf("%s.mu", params.prospect5)
fft.sub <- fft.f[!is.na(N.mu)][sensor=="identity"][,vars,with=F]
quants <- fft.sub[, lapply(.SD, quantile, c(0.05, 0.95))]
bounds <- function(x, n) (x > quants[[n]][1]) & (x < quants[[n]][2])
fft.trim <- fft.sub[bounds(N.mu, "N.mu")
                    ][bounds(Cab.mu, "Cab.mu")
                    ][bounds(Car.mu, "Car.mu")
                    ][bounds(Cw.mu, "Cw.mu")
                    ][bounds(Cm.mu, "Cm.mu")]
par.mat <- as.matrix(fft.trim)
npar <- nrow(par.mat)
noise.mat <- sapply(1:npar, function(x) generate.noise())
spec.mat <- apply(par.mat, 1, function(x) prospect(x, 5)[,1])
obs.mat <- spec.mat + noise.mat
save(par.mat, spec.mat, noise.mat, obs.mat, file="run-scripts/simulation.inputs.RData")

# Write submission script
n.sensor <- length(sensor.list)
index <- rep(1:npar, each=n.sensor)
sensor <- rep(sensor.list, npar)
ngibbs <- 50000
runname <- sprintf("%s.%d", sensor, rep(1:npar, each=n.sensor))
submit.string <- sprintf('qsub -V -v index=%d,sensor=%s,ngibbs=%d,runname=%s -N "%s" submit.simulated.qsub',
                         index, sensor, ngibbs, runname, runname)

# Create file
fname <- "run-scripts/run.simulation.sh"
write("#!/bin/bash", file=fname)
write(submit.string, file=fname, append=TRUE)
