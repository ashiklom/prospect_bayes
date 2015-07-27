#' Script to load all runs

parnames <- c("N", "Cab", "Car", "Cw", "Cm", "residual")

results.path <- "results-simulation"
flist.full <- list.files(results.path, ".*.RData")
results.list <- list()
for(f in flist.full){
# Extract data
    print(f)
    f.path <- file.path(results.path, f)
    load(f.path)
    f.name <- gsub("(.*)[.]RData", "\\1", f)
    f.list <- get(f.name)
    rm(list=f.name)
    sensor <- list(sensor = f.list$sensor)
    samples <- f.list$samples
    true.param <- as.list(f.list$true.param)
    names(true.param) <- sprintf("%s.true", parnames[1:5])

# Burn in and thin
    ngibbs <- nrow(samples)
    burnin <- floor(ngibbs/2)
    thin <- (ngibbs - burnin) / 5000
    bt <- seq(burnin, ngibbs, by=thin)
    samples.sub <- samples[bt,]
    
# Calculate summary statistics
    mu <- colMeans(samples.sub, na.rm=TRUE)
    names(mu) <- sprintf("%s.mu", parnames)
    sigma <- apply(samples.sub, 2, sd, na.rm=TRUE)
    names(sigma) <- sprintf("%s.sigma", parnames)
    q25 <- apply(samples.sub, 2, quantile, 0.025, na.rm=TRUE)
    names(q25) <- sprintf("%s.q25", parnames)
    med <- apply(samples.sub, 2, median, na.rm=TRUE)
    names(med) <- sprintf("%s.med", parnames)
    q975 <- apply(samples.sub, 2, quantile, 0.975, na.rm=TRUE)
    names(q975) <- sprintf("%s.q975", parnames)
    dat.row <- c(true.param = true.param,
                 sensor = sensor,
                 fname = f.name,
                 as.list(mu), as.list(sigma),
                 as.list(q25), as.list(med),
                 as.list(q975))
    results.list[[f]] <- dat.row
}

#' Combine into data.table
library(data.table)
fft.dat <- do.call(rbindlist, list(results.list)) 
save(fft.dat, file="../data/simulation.dat.RData")
