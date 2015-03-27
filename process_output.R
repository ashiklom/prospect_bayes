## Script to process FFT individual leaf results
library(data.table)
library(rjags)
load("run_results/FFTdata.Rdata")
spec.regex <- "^(.*)_(.*)_(.*)_LC_REFL_(.*)"

## Some statistics based on leaf means
leaf.means <- all.r[,lapply(.SD, mean), by = Spectra]
leaf.means[, Site.plot := gsub(spec.regex, "\\1", Spectra)]
leaf.means[, Species := gsub(spec.regex, "\\2", Spectra)]
leaf.means[, HeightAge := gsub(spec.regex, "\\3", Spectra)]
