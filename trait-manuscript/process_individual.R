### Diagnostics for plot
library(coda)
library(data.table)
f.path <- "run_results/FFT_individuals_0323/"
f.list <- list.files(f.path)
spec.list <- unique(gsub("(.*)_[0-9]{2}\\.dat", "\\1", f.list))
nspec <- length(spec.list)
blank <- double(nspec)
load.spec <- function(specname, burnin=20000, thin=30){
        in.files <- list.files(f.path, specname, full.names=TRUE)
        r.list <- mcmc.list()
        for(i in 1:length(in.files)){
                in.list <- mcmc(fread(in.files[i], header=TRUE))
                r.list[[i]] <- as.matrix(window(in.list, start=burnin, thin=thin))
        }
        r.out <- data.table(do.call(rbind, r.list))
        r.out[, Spectra := specname]
        setcolorder(r.out, c("Spectra", "N", "Cab", "Cw", "Cm", "rsd"))
        return(r.out)
}
### Read all samples
system.time(all.samples <- lapply(spec.list, load.spec))

### Merge Mean and SD from all.summary with
### fitted values from all.table into "results"
as2 <- data.table(t(sapply(all.summary,
                           function(l) l$statistics[,c("Mean", "SD")])),
                  keep.rownames=TRUE)
summary.vars <- c("N.m", "Cab.m", "Cw.m", "Cm.m", "resp",
                  "N.sd", "Cab.sd", "Cw.sd", "Cm.sd", "resp.sd")
setnames(as2, c("leaf", summary.vars))
results <- merge(as2, all.table, by="leaf")
save(results, file = "results_FFT0203.Rdata")