#' Script to load all runs
library(PEcAnRTM)
parnames <- c(params.prospect5, "residual")

results.path <- "results-simulation"
flist.full <- list.files(results.path, ".*.RData")
results.list <- list()
for(f in flist.full){
# Extract data
    print(f)
    f.list <- load.from.name(f, filepath=results.path)
    sensor <- list(sensor = f.list$sensor)
    samples <- f.list$samples
    true.param <- as.list(f.list$true.param)
    names(true.param) <- sprintf("%s", parnames[1:5])

# Burnin, thin, and summary statistics
    samples.sub <- burnin.thin(samples)
    samples.summary <- summary.mvnorm(samples.sub)
    dat.row <- c(true.param = true.param,
                 sensor = sensor,
                 fname = f.name,
                 samples.summary)
    results.list[[f]] <- dat.row
}

#' Combine into data.table
library(data.table)
simulation.dat <- do.call(rbindlist, list(results.list)) 
save(simulation.dat, file="../data/simulation.samp.dat.RData")
