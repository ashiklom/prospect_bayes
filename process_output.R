## Script to process FFT individual leaf results
library(data.table)
library(coda)

fft.path <- "run_results/FFT_individuals_0323/"
fft.list <- list.files(fft.path)