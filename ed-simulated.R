# Spectra simulated using PROSPECT

# Load FFT data
source("preprocess.fft.R")

# Fit multivariate normal to data and simulate {{{
fit.mv <- function(dat, n=100, logval=FALSE){
    require(mclust)
    require(MASS)
    fit <- mvn("XXX", as.matrix(dat))
    mu <- fit$parameters$mean
    sigma <- fit$parameters$variance$Sigma
    samp <- mvrnorm(n, mu, sigma)
    nl <- 0
    print("Sampling")
    while(TRUE){
        if(logval) samp.conds <- apply(samp, 1, function(x) all(exp(x) > c(1,0,0,0)))
        else samp.conds <- apply(samp, 1, function(x) all(x > c(1,0,0,0)))
        samp <- samp[samp.conds,]
        nl <- dim(samp)[1]
        if(nl == n) break
        nlm1 <- n - nl
        print(nlm1)
        samp2 <- mvrnorm(nlm1, mu, sigma)
        samp <- rbind(samp, samp2)
    }
    return(samp)
}
# }}}

# Generate PROSPECT spectra from samples {{{
#   1:2101 is reflectance, 2102:4202 is transmittance
spec.sim <- function(samp, PAR=400:700, NIR=800:2500){
    require(PEcAnRTM)
    PAR.refl <- PAR-399
    NIR.refl <- NIR-399
    PAR.trans <- PAR-399+2101
    NIR.trans <- NIR-399+2101
    calc.parnir <- function(x) c("PAR.refl"=mean(x[PAR.refl], na.rm=TRUE),
                                 "NIR.refl"=mean(x[NIR.refl], na.rm=TRUE),
                                 "PAR.trans"=mean(x[PAR.trans], na.rm=TRUE),
                                 "NIR.trans"=mean(x[NIR.trans], na.rm=TRUE))
    spec.list <- apply(samp, 1, prospect, version=4)
    pn.list <- apply(spec.list, 2, calc.parnir)
    return(pn.list)
}
#}}}

# Perform calculations {{{
vars <- c("N.mu", "Cab.mu", "Cw.mu", "Cm.mu")
succ <- c("early", "mid", "late")
plant.types <- c("hardwood", "conifer")
PFT <- sprintf("%s %s", rep(succ, 2), rep(plant.types, each=3))

n = 5000
logval = FALSE
pft.list <- list()
spec.list <- list()
summary.func <- function(x) c("mean"=mean(x), "sd"=sd(x), 
                                 "min"=min(x), "max"=max(x), "median"=median(x), 
                                 "q2.5"=quantile(x, 0.025), "q97.5"=quantile(x, 0.975))
for(p in PFT){
    dat <- fft[PFT == p, vars, with=F]
    fit <- fit.mv(dat, n=n)
    spec <- spec.sim(fit)
    spec.list[[p]] <- spec
    pft.list[[p]] <- apply(spec, 1, summary.func)
}
pft.dat <- as.data.frame(t(as.data.frame(pft.list)))
pft.dat$variable <- rownames(pft.dat)
pft.dat <- data.table(pft.dat)
rxp <- paste(rep("([[:alpha:]]+)", 4), collapse="\\.")
pft.dat[, succession := gsub(rxp, "\\1", variable)]
pft.dat[, plant.type := gsub(rxp, "\\2", variable)]
pft.dat[, PFT := paste(succession, plant.type)]
pft.dat[, region := gsub(rxp, "\\3", variable)]
pft.dat[, measure := gsub(rxp, "\\4", variable)]
pft.dat[, variable := NULL]
write.csv(pft.dat, file="other-figures/par-nir-summary.csv", row.names=FALSE)
save(pft.list, spec.list, file="other-figures/par-nir-data.RData")
# }}}

