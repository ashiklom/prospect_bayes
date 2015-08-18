#' Perform inversion (in sequence) of new chlorophyll data
library(PEcAnRTM)
library(mclust)
load("../../data/chl-raw.RData")

nspec <- nrow(chl.raw)
waves <- sprintf("Wave_%d", 400:2500)
ngibbs <- 10000

p5.sample <- function(spec){
    model <- function(param, const) prospect(param, 5)[,1]
    prior.vals <- prior.defaultvals.prospect()
    prior <- priorfunc.prospect(prior.vals$mu, prior.vals$sigma)
    inits <- rlnorm(5, prior.vals$mu, prior.vals$sigma)
    inits[1] <- inits[1] + 1
    samples <- invert.slow(observed = spec,
                           inits = inits,
                           constants = NULL,
                           ngibbs = ngibbs,
                           prior = prior,
                           pm = c(1,0,0,0,0),
                           model = model,
                           do.mle = TRUE,
                           quiet = TRUE)
    return(samples)
}

burnin.thin <- function(samples){
    ngibbs <- nrow(samples)
    burnin <- floor(ngibbs / 2)
    thin <- floor((ngibbs - burnin) / 5000)
    samples.bt <- samples[seq(burnin, ngibbs, by=thin),]
    return(samples.bt)
}

parnames <- c("N", "Cab", "Car", "Cw", "Cm", "residual")
sigmanames <- sprintf("%s-%s.sigma", rep(parnames, 1, each=6),
                      rep(parnames, 6))
samples.summary <- function(samples.bt){
    fit <- mvn("XXX", samples.bt)
    mu <- as.numeric(fit$parameters$mean)
    sigma <- c(fit$parameters$variance$Sigma)
    names(mu) <- sprintf("%s.mu", parnames)
    names(sigma) <- sigmanames
    out.list <- c(as.list(mu), as.list(sigma))
    return(out.list)
}

chl.dat <- chl.raw[,1:29, with=FALSE]
for(i in 1:nspec){
    dat.row <- chl.raw[i]
    specname <- dat.row[, Spectra]
    fname <- sprintf("chl-results/%s.RData", specname)
    print(c(i, specname))
    spec <- unlist(dat.row[,waves, with=FALSE])
    samples <- p5.sample(spec)
    save(samples, file=fname)
    samples.bt <- burnin.thin(samples)
    samples.calc <- samples.summary(samples.bt)
    chl.dat[Spectra == specname, 
            names(samples.calc) := samples.calc, with=FALSE]
}

save(chl.dat, file="../data/chl-dat.RData")
