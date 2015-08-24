# Noise run
library(PEcAnRTM)
arg <- commandArgs(trailingOnly=TRUE)
param <- as.numeric(arg[1:5])
names(param) <- paste0(params.prospect5, ".true")
fw <- as.numeric(arg[6])
sigma <- as.numeric(arg[7])
rnum <- as.numeric(arg[8])

spec <- prospect(param, 5)[,1]
noise <- generate.noise(n=2101, fw=fw, sigma=sigma, fsd=2)
obs <- spec + noise

samples <- default.invert.prospect(obs, "identity", 50000, version=5, 
                                   do.mle=FALSE, quiet=TRUE)
colnames(samples) <- c(params.prospect5, "residual")
samples.sub <- burnin.thin(samples)
parmeans <- colMeans(samples.sub)
names(parmeans) <- paste0(params.prospect5, ".mu")
parsd <- apply(samples.sub, 2, sd)
names(parsd) <- paste0(params.prospect5, ".sigma")

runname <- sprintf("noise.%d.%d.%d", rnum, fw, sigma*10000)
save.p2 <- list(fw=fw, sigma=sigma, samples=samples)
save.list <- c(as.list(parmeans), as.list(parsd),
               as.list(param), save.p2)
assign(runname, save.list)
save(list=runname, file=sprintf("../results-noise/%s.RData", runname))
