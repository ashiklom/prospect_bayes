#' Generate list of FFT spectra
library(data.table)
FFT.path <- "../../data/FFT_spectra/FFT_spectra_unique.csv"
all.spec <- fread(FFT.path, header=TRUE)
spec.list <- all.spec[, Spectra]

write.table(spec.list, file="fft-speclist.dat",
            row.names=FALSE, col.names=FALSE)

