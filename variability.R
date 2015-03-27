### Variability
library(rjags)
library(data.table)
load("data/FFT_full.Rdata")
fft.spec <- fft.spec[!is.na(Spectra) & !is.na(PFT)]
design.reg <- model.matrix(N ~ Height + PFT + Label, data=fft.spec)
burnin <- 5000
nchains <- 1
thin <- 1
niter <- 5000 * thin / nchains
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
              n.all = nrow(fft.spec))
vdata$y.m <- fft.spec[,N]
vmodel <- jags.model(file = textConnection(vcode),
                     data = vdata,
                     n.chains = nchains)
update(vmodel, burnin)
monitors <- c("Beta", "tau")
z <- system.time(vsamples <- coda.samples(vmodel, monitors, niter, thin))