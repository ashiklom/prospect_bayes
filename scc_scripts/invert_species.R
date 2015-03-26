##' Run of Bayesian inversion of PROSPECT
##'
##' For initializing inversions from the command line, e.g. for submission to cluster.
##' This script takes 8 command line arguments, which are as follows:
##'   [1] Spectype.Species (e.g. SE.QUCH)
##'   [2] Random effects; current options are 'none' and 'leaf'
##'   [3] How to generate initial conditions; OPTIONS: 'random', 'mle' (R optim function for mean observed spectrum), 'guess' (preset values)
##'   [4] Number of MCMC steps.
##'   [5] Filename tag to identify run results
##'   [6] Sub-folder in "run_results" for storing output.

library(PEcAnRTM)
library(data.table)
options(warn=2)

if(exists("TEST")){
        args <- c("ACRU", "100", "ACRU_species_test", "test")
} else {
        args <- commandArgs(trailingOnly=TRUE)
}
species <- args[1]
ngibbs <- as.numeric(args[2])
folder <- args[3]
run_id <- args[4]

filename <- sprintf("../run_results/%s/%s_%s.dat",
                    folder,
                    spectra,
                    run_id)
dir.create(sprintf("../run_results/%s", folder), recursive=TRUE, showWarnings = FALSE)

FFT.path <- "../data/FFT_spectra/FFT_spectra_unique.csv"
all.spec <- fread(FFT.path, header=TRUE)
obs.spec <- t(as.matrix(all.spec[Species == species, 72:2172, with=FALSE]))

result <- invert_prospect_re(obs.spec, ngibbs)
write.csv(result, filename, row.names=FALSE)

