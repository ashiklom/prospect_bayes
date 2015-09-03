#' ---
#' title: Load results of simulation experiment
#' author: Alexey Shiklomanov
#' output_format: html_document
#' ---

#' The script used to perform inversions on the BU cluster returns an RData 
#' file containing a list of the matrix of MCMC samples, the name of the 
#' sensor, the ID of the run, and pre-computed summary statistics. This script 
#' loops over the list of all simulation output files, loads each, extracts the 
#' samples, computes the relevant summary statistics, and adds the results to a 
#' list that is then melded into a `data.table` and saved.

library(PEcAnRTM)
results.path <- "results-simulation"
flist.full <- list.files(results.path, ".*.RData")
results.list <- list()
for(f in flist.full){
    print(f)
    dat.row <- load.from.name(f, filepath=results.path)
    samples <- dat.row$samples
    samples.bt <- burnin.thin(samples)
    samples.summary <- summary.simple(samples.bt)
    dat.sub <- dat.row[c("sensor", "fname")]
    results.list[[f]] <- c(dat.sub, samples.summary)
}

#' This step uses the `rbindlist` function from the `data.table` package for 
#' very efficient row binding of multiple list objects. 
library(data.table)
simulation.dat <- do.call(rbindlist, list(results.list)) 
save(simulation.dat, file="../data/simulation.samp.dat.RData")
