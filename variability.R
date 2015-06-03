### Variability
library(rjags)
source("preprocess.fft.R")

arg <- commandArgs(trailingOnly=TRUE)
v <- arg[1]
dat <- get(arg[2])

### Run models
vm <- paste0(v, ".mu")
vsd <- paste0(v, ".sd")
vlm <- paste0(v, ".lmu")
vlsd <- paste0(v, ".lsd")

f <- sprintf("%s ~ Height + PFT + Label + Site + Plot + Family", vm)
design.reg <- model.matrix(formula(f), data=dat)

burnin <- 50
nchains <- 20
nthin <- 100
ngibbs <- 5 * nthin / nchains

vcode <- "
model{
        ### Hierarchical random-effect, global mean
        for(i in 1:n.all){
                Ex[i] <- X[i,] %*% Beta
                y[i] ~ dnorm(Ex[i], tau)
                y.m[i] ~ dnorm(y[i], tau_obs[i])
        }

        ### Priors
        for(i in 1:nfe){Beta[i] ~ dnorm(0, 0.01)}
        tau ~ dgamma(0.01, 0.01)
        sd <- 1/tau
}
"

vdata <- list(X = design.reg,
              nfe = ncol(design.reg),
              n.all = nrow(dat))
vdata$y.m <- dat[[vm]]
vdata$tau_obs <- 1/dat[[vsd]]^2

vmodel <- jags.model(file = textConnection(vcode),
                     data = vdata,
                     n.chains = nchains)

update(vmodel, burnin)
monitors <- c("Beta", "sd")
z <- system.time(vsamples <- coda.samples(vmodel, monitors, ngibbs, nthin))

fname <- sprintf("SAMPLES %s %s %s.RData", Sys.time(), arg[1], arg[2])
save(vsamples, file=fname)
