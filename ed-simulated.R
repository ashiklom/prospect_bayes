# Spectral priors for ED runs

# Convert summary statistics to lognormal parameters
ln.mu <- function(m,v) log(m / sqrt(1 + v/m^2))
ln.sigma <- function(m,v) sqrt(log(1 + v/m^2))

# Calculate PAR and NIR for reflectance and transmittance
parnir <- function(spec){
    parwl.r <- 400:700 - 399
    nirwl.r <- 800:2500 - 399
    parwl.t <- parwl.r + 2101
    nirwl.t <- nirwl.r + 2101
    par.r <- mean(spec[parwl.r])
    nir.r <- mean(spec[nirwl.r])
    par.t <- mean(spec[parwl.t])
    nir.t <- mean(spec[nirwl.t])
    out <- c("par.r" = par.r, "nir.r" = nir.r,
             "par.t" = par.t, "nir.t" = nir.t)
    return(out)
}

# Simulate spectra
simspec <- function(means, vars, n=5000){
    require(PEcAnRTM)
    mpars <- ln.mu(means,vars)
    vpars <- ln.sigma(means,vars)
    samps <- rlnorm(n, mpars, vpars)
    params <- matrix(samps, n, 5, byrow=TRUE)
    spec <- apply(params, 1, prospect, version=5)
    return(spec)
}
    
# Hardwood summary stats
# From Feret et al 2011 RSE 
mean.h <- c(1.7, 32.81, 8.51, 0.0129, 0.0077)
var.h <- c(0.6, 17.87, 3.2, 0.0073, 0.0035)^2

# Conifer summary stats
# N (1), Cw (4), Cm (5)
#   Means from Croft et al 2014 Ecological Complexity 17:119-130
#   SD based on ranges in Croft et al 2015 ISPRS J of Photogrammetry and Remote Sensing 102:85-85
#       (set to SD = range/2)
# Cab (2), Car (3) 
#   From Di Vittorio 2009 RSE 113:948-1966
mean.c <- c(2.6, 68.6, 12.1, 0.001, 0.035)
var.c <- c(2.5, 11.9, 40, 0.02, 0.025)^2

# Run simulations
n = 5000
spec.h <- simspec(mean.h, var.h, n)
pn.h <- apply(spec.h, 2, parnir)

spec.c <- simspec(mean.c, var.c, n)
pn.c <- apply(spec.c, 2, parnir)

# Stat summary function
statsummary <- function(x) c("mean" = mean(x, na.rm=TRUE), "sd" = sd(x, na.rm=TRUE), 
                             "min" = min(x, na.rm=TRUE), "max" = max(x, na.rm=TRUE), 
                             "q" = quantile(x, c(0.025, 0.5, 0.975), na.rm=TRUE))
summary.h <- apply(pn.h, 1, statsummary)
summary.c <- apply(pn.c, 1, statsummary)
print(summary.h)
print(summary.c)
write.csv(summary.h, file="other-figures/hardwood-parnir.csv")
write.csv(summary.c, file="other-figures/conifer-parnir.csv")

# Summary plot
library(ggplot2)
library(reshape2)
library(data.table)
library(gridExtra)
pnd.h <- melt(data.frame(t(pn.h), id=1:nrow(pn.h)), id.vars="id")
plot.h <- ggplot(pnd.h) + aes(x=variable, y=value) + geom_violin() + xlab("hardwood")
pnd.c <- melt(data.frame(t(pn.c), id=1:nrow(pn.c)), id.vars="id")
plot.c <- plot.h %+% pnd.c + xlab("conifer")

png("other-figures/summary-figure.png")
grid.arrange(plot.h, plot.c, nrow=2)
dev.off()

