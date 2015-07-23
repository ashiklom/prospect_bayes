#' Run of Bayesian inversion of PROSPECT
#' Arguments: spectra ngibbs sensor

library(PEcAnRTM)
options(warn=2)

#' Parse arguments
args <- commandArgs(trailingOnly=TRUE)
spectra <- args[1]
ngibbs <- as.numeric(args[2])
sensor <- args[3]

#' Load data
FFT.path <- "../../data/FFT_spectra/FFT_spectra_unique.csv"
all.spec <- fread(FFT.path, header=TRUE)
obs <- as.matrix(as.numeric(all.spec[Spectra == spectra, 72:2172, with=FALSE])) 

#' Setup inversion parameters
parnames <- c("N", "Cab", "Car", "Cw", "Cm")
prior.vals <- prior.defaultvals.prospect(sd.inflate = 3)
start.params <- rlnorm(5, prior.vals$mu, prior.vals$sigma)
names(start.params) <- parnames
start.params["N"] <- start.params["N"] + 1
pm <- c(1,0,0,0,0)
model <- function(params, constants) spectral.response(prospect(params, 5)[,1], sensor)
prior <- priorfunc.prospect(prior.vals$mu, prior.vals$sigma)

#' Perform inversion
samples <- invert.slow(observed = obs,
                       inits = start.params,
                       constants = NULL,
                       ngibbs = ngibbs,
                       prior = prior,
                       pm = pm,
                       model = model,
                       do.mle = TRUE,
                       quiet = TRUE)

#' Export and save
run.name <- sprintf("%s.%s", spectra, sensor)
assign(run.name, list(spectra, start.params, sensor, samples))
save(list = run.name, file=sprintf("../results/%s.RData", run.name)

