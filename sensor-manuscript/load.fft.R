#' Script to load all runs
library(PEcAnRTM)
parnames <- c(params.prospect5, "residual")

results.path <- "results"
flist.full <- list.files(results.path, ".*.RData")
results.list <- list()
for(f in flist.full){
# Extract data
    print(f)
    f.list <- load.from.name(f, filepath=results.path)
    sensor <- f.list$sensor
    samples <- f.list$samples
    spectra <- f.list$spectra
    samples.sub <- burnin.thin(samples,
                               auto=TRUE,
                               burnin.ratio=2)
    
# Burnin, thin, and summary statistics
    samples.summary <- summary.mvnorm(samples.sub)
    dat.row <- c(list(spectra = spectra, 
                      sensor = sensor,
                      fname = f.name),
                 samples.summary)
    results.list[[f]] <- dat.row
}

#' Combine into data.table
library(data.table)
fft.dat <- do.call(rbindlist, list(results.list)) 
save(fft.dat, file="fft.dat.RData")
