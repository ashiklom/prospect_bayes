# Generate runs for spectral noise
library(PEcAnRTM)

p5 <- prior.defaultvals.prospect()
sim.p5 <- function(){
    r <- rlnorm(5, p5$mu, p5$sigma)
    r[1] <- r[1] + 1
    return(r)
}

nruns <- 1:100
fsd <- 2
fw <- c(3,5,7,9,11,13,15,17,19,21,31,41,51,101,201)
sigma <- c(1e-4, 2.5e-4, 5e-4, 7.5e-4, 1e-3, 1.5e-3,
           2e-3, 3e-3, 5e-3, 1e-2) 

fname <- "run-scripts/run-noise.sh"
submit.string <- 'qsub -V -v N=%.3f,Cab=%.3f,Car=%.3e,Cw=%.3e,Cm=%.3e,fw=%f,sigma=%f,n=%d submit.noise.qsub'
write("#!/bin/bash", file=fname)
for(n in nruns){
    for(f in fw){
        for(s in sigma){
            param <- sim.p5()
            s <- sprintf(submit.string, param[1], param[2], param[3],
                         param[4], param[5], f, s, n)
            write(s, file=fname, append=TRUE)
        }
    }
}

