### Variability
library(rjags)
source("preprocess.fft.R")
design.reg <- model.matrix(N.mu ~ Height + PFT + Label + Site + Plot,
                           data=fft.h)

burnin <- 5000
nchains <- 1
nthin <- 1
ngibbs <- 5000 * nthin / nchains
vcode <- "
model{
        ### Hierarchical random-effect, global mean
        for(i in 1:n.all){
                Ex[i] <- X[i,] %*% Beta
                y.m[i] ~ dnorm(Ex[i], tau)
        }

        ### Priors
        for(i in 1:nfe){Beta[i] ~ dnorm(0, 0.01)}
        tau ~ dgamma(0.01, 0.01)
}
"
vdata <- list(X = design.reg,
              nfe = ncol(design.reg),
              n.all = nrow(fft.h))
vdata$y.m <- fft.spec[,N.mu]
vmodel <- jags.model(file = textConnection(vcode),
                     data = vdata,
                     n.chains = nchains)
update(vmodel, burnin)
monitors <- c("Beta", "tau")
z <- system.time(vsamples <- coda.samples(vmodel, monitors, ngibbs, nthin))
