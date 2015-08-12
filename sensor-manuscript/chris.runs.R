# Redo CHRIS PROBA runs
library(data.table)
load("../data/simulation.samp.dat.RData")
cp <- simulation.dat[sensor.sensor == "chris.proba"]
cp.bad <- cp[abs(true.param.N.true - N.mu) > 0.05]

N <- cp.bad[,true.param.N.true]
Cab <- cp.bad[,true.param.Cab.true]
Car <- cp.bad[,true.param.Car.true]
Cw <- cp.bad[,true.param.Cw.true]
Cm <- cp.bad[,true.param.Cm.true]
sensor <- "chris.proba"
ngibbs <- 200000
runname <- cp.bad[,fname]
submit.string <- sprintf('qsub -V -v N=%f,Cab=%f,Car=%f,Cw=%f,Cm=%f,sensor=%s,ngibbs=%d,runname=%s -N "%s" submit.simulated.qsub',
                         N, Cab, Car, Cw, Cm, sensor, ngibbs, runname, runname)

fname <- "run-scripts/run.chrisbad.sh"
write("#!/bin/bash", file=fname)
write(submit.string, file=fname, append=TRUE)
