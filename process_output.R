## Script to process FFT individual leaf results
library(data.table)
library(coda)

fft.path <- "run_results/FFT_individuals_0323/"
fft.list <- list.files(fft.path)

bt.list <- list()
for(f in fft.list){
        fname <- sprintf("%s%s", fft.path, f)
        infile <- fread(fname, header=TRUE)
        b <- infile[80000:100000,,with=FALSE]