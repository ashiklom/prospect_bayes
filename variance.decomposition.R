library(coda)
library(data.table)
library(reshape2)
library(ggplot2)
load("data/jags-results.RData")

# Unchain mcmc.list's
all.mat <- lapply(all.mcmc, as.matrix)

# Get names of Betas
for(n in names(all.mat)){
    nl <- tolower(n)
    fname <- sprintf("data/design-mat/%s.R", nl)
    source(fname)
    cnames <- sprintf("%s-%s", nl, c(X, "tau"))
    colnames(all.mat[[n]]) <- cnames
}

# Combine into single matrix
all.dat <- do.call(cbind, all.mat)

# Calculate summary statistics
all.means <- colMeans(all.dat)
all.sd <- apply(all.dat, 2, sd)
all.quant <- t(apply(all.dat, 2, quantile, c(0.025, 0.975)))
all.sum <- cbind(all.means, all.sd, all.quant)

# Create data.table
rnames <- rownames(all.sum)
ams.dat <- data.table(all.sum, varname=rnames)
setkey(ams.dat, varname)

# Extract information from varname
vrxp <- "([[:alpha:]]+)-(fft[.]?[hc]?)-(.*)"
ams.dat[, pvar := gsub(vrxp, "\\1", varname)]
ams.dat[, dat := gsub(vrxp, "\\2", varname)]
ams.dat[, fullbeta := gsub(vrxp, "\\3", varname)]
preds <- c("PFT", "Height", "Label", "Site", "Plot", "Family")
preds.rxp <- sprintf("(%s)(.*)", paste(preds, collapse="|"))
ams.dat[, pred := gsub(preds.rxp, "\\1", fullbeta)]

# Perform variance decomposition
var.t <- function(x){
    if(length(x) > 1) return(var(x))
    else return(1/x)
}
partial.variance <- ams.dat[pred != "(Intercept)", list(pv =var.t(all.means)), by=c("pvar", "dat", "pred")]
total.variance <- partial.variance[, list(tv = sum(pv)), by=c("pvar", "dat")]
vd <- merge(partial.variance, total.variance, by=c("pvar", "dat"))
vd[, fv := pv/tv]

# Stacked bar plot
plt <- ggplot(vd) + aes(x=pvar, y=fv, fill=pred) + 
    geom_bar(stat="identity") +
    facet_grid(dat ~ ., scales="free_y")
png("manuscript/figures/variance-decomposition.png", 
    height=5, width=4, units="in", res=300)
plot(plt)
dev.off()


