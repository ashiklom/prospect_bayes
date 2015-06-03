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
wl.names <- sprintf("Wave_%d", 400:2500)
t.names <- names(trans)[1:21]

t.names[t.names %in% names(fft)]
setkey(fft, Sample_Name)
setkey(trans, Sample_Name)
rt.big <- fft[trans]

tdiff <- function(mrow){
    pars <- mrow[1:4]
    obs <- mrow[-4:0]
    if(any(is.na(pars))) return(rep(NA,2101))
    tmodel <- prospect(pars,4)[,2]
    tdiff <- tmodel - obs
    return(tdiff)
}

rt.mat <- as.matrix(rt.big[,c("N.mu", "Cab.mu", "Cw.mu", "Cm.mu", wl.names), with=FALSE])
trans.mat <- matrix(NA, nrow(rt.big), 2101)
for(i in 1:nrow(rt.big))
    trans.mat[i,] <- tdiff(rt.mat[i,])

trans.means <- colMeans(trans.mat, na.rm=TRUE)
trans.q <- apply(trans.mat, 2, quantile, c(0.025,0.1,0.9,0.975), na.rm=TRUE)
trans.dat <- data.table(wavelength = 400:2500, means=trans.means)
trans.dat <- cbind(trans.dat, t(trans.q))
setnames(trans.dat, rownames(trans.q), c("low2", "low1", "high1", "high2"))

## Transmittance error plot
require(ggplot2)
trans.spec <- ggplot(trans.dat) +
    aes(x=wavelength) +
    geom_line(aes(y=means)) +
    geom_ribbon(aes(ymin=low1,ymax=high1), alpha=0.3, fill="red") +
    geom_ribbon(aes(ymin=low2,ymax=high2), alpha=0.3) +
    xlab("Wavelength(nm)") +
    ylab("Mean transmittance error") +
    theme_bw() +
    theme(axis.title = element_text(size=10),
          axis.text = element_text(size=7))
png("manuscript/figures/transmittance-validation.png", 
    width=4, height=3, units="in", res=300)
plot(trans.spec)
dev.off()


## Manuscript plot -- violin
require(ggplot2)
load("transmittance.RData")
th <- theme_bw() + theme()
theme_set(th)

png("manuscript/figures/transmittance-validation.png", width=4, height=3, units="in", res=300)
plot.dat <- rt.big[plant.type %in% c("hardwood", "conifer"),]
r1 <- ggplot(plot.dat) + th +
    aes(x = plant.type, y = trans.rmse) +
    geom_violin() +
    xlab("Plant functional type") +
    ylab("Transmittance RMSE")
plot(r1)
dev.off()
