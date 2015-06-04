library(BayesFactor)
source("preprocess.fft.R")

# Attempt ANOVA with BayesFactor package
# Convert characters to factors (required for anovaBF
to.factor <- function(x) {
    if(is.character(x)) return(factor(x))
    else return(x)
}
fft.f <- fft[, lapply(.SD, to.factor)]

bf.N <- anovaBF(N.mu ~ PFT + Height + Family,# + Label + Site + Plot,
                data=fft.f, progress=TRUE)

bf.samples <- posterior(max(bf.N), iterations=1000)

# Try BANOVA package
library(BANOVA)
banova.N <- BANOVA.Normal(N.mu ~ 1, ~ PFT + Height, fft.f,
                          fft.f$Spectra, burnin=50, sample=500, thin=1)
av.N <- BAnova(banova.N)

