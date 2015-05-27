### Calculate RMSE
library(PEcAnRTM)

source("preprocess.fft.R")

# Load transmittance data
trans.path <- file.path("~", "Documents", "Dropbox",
                        "NASA_TE_PEcAn-RTM_Project",
                        "Data", "Spectra",
                        "NASA_FFT_Spectra",
                        "NASA_FFT_IS_Tran_Spectra_v4.csv")
trans <- fread(trans.path, header=TRUE)
wl.names <- sprintf("Wave_%d", 350:2500)
t.names <- names(trans)[1:21]

t.names[t.names %in% names(fft)]
setkey(fft, Sample_Name)
setkey(trans, Sample_Name)
rt.big <- fft[trans]

tdiff <- function(pars, obs){
    # pars <- with(dtrow, c(N.mu, Cab.mu, Cw.mu, Cm.mu))
    # obs <- with(dtrow, get(wl.names))
    if(any(is.na(pars))) return(NA)
    tmodel <- prospect(pars, 4)
    tdiff <- tmodel - obs[-50:0]
    rmse <- sqrt(mean(tdiff^2))
    return(rmse)
}

rms <- numeric(nrow(rt.big))
for (i in 1:nrow(rt.big)){
    pars <- unlist(rt.big[i, list(N.mu, Cab.mu, Cw.mu, Cm.mu)])
    obs <- unlist(rt.big[i, wl.names, with=FALSE])
    rms[i] <- try(tdiff(pars, obs),silent=TRUE)
}

rms <- as.numeric(rms)
rt.big[, trans.rmse := rms]

save(rt.big, file="transmittance.RData")

## Manuscript plot
library(ggplot2)
load("transmittance.RData")
th <- theme_bw() + theme()
theme_set(th)

png("manuscript/figures/transmittance.validation.png", width=4, height=3, units="in", res=300)
plot.dat <- rt.big[plant.type %in% c("hardwood", "conifer"),]
r1 <- ggplot(plot.dat) + th +
    aes(x = plant.type, y = trans.rmse) +
    geom_violin() +
    xlab("Plant functional type") +
    ylab("Transmittance RMSE")
plot(r1)
dev.off()
