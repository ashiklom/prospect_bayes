#' Script to load all runs
library(PEcAnRTM)
results.path <- "results-simulation"
flist.full <- list.files(results.path, ".*.RData")
results.list <- list()
for(f in flist.full){
# Extract data
    print(f)
    dat.row <- load.from.name(f, filepath=results.path)
    dat.row[c("noise", "obs.raw", "samples")] <- NULL
    results.list[[f]] <- dat.row
}

#' Combine into data.table
library(data.table)
simulation.dat <- do.call(rbindlist, list(results.list)) 
save(simulation.dat, file="../data/simulation.samp.dat.RData")
