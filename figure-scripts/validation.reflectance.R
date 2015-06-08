#' Plot for validation based on reflectance error
library(PEcAnRTM)

if(grepl("figure-scripts", getwd())) setwd("..")
source("preprocess.fft.R")

# Load reflectance data
refl.path <- file.path("~", "Documents", "Dropbox",
                        "NASA_TE_PEcAn-RTM_Project",
                        "Data", "Spectra",
                        "NASA_FFT_Spectra",
                        "NASA_FFT_LC_Refl_Spectra_v4.csv")
refl <- fread(refl.path, header=TRUE)
wl.names <- sprintf("Wave_%d", 400:2500)
r.names <- names(refl)
r.names <- r.names[!(grepl("Wave_", r.names))]
r.names <- r.names[r.names %in% names(fft)]

setkey(fft, Sample_Name, Sample_Year)
setkey(refl, Sample_Name, Sample_Year)
rt.big <- fft[refl]

rdiff <- function(mrow){
    pars <- mrow[1:4]
    obs <- mrow[-4:0]
    if(any(is.na(pars))) return(rep(NA,2101))
    rmodel <- prospect(pars,4)[,1]
    rdiff <- rmodel - obs
    return(rdiff)
}

rt.mat <- as.matrix(rt.big[,c("N.mu", "Cab.mu", "Cw.mu", "Cm.mu", wl.names), with=FALSE])
refl.mat <- matrix(NA, nrow(rt.big), 2101)
for(i in 1:nrow(rt.big))
    refl.mat[i,] <- rdiff(rt.mat[i,])

refl.means <- colMeans(refl.mat, na.rm=TRUE)
refl.q <- apply(refl.mat, 2, quantile, c(0.025,0.1,0.9,0.975), na.rm=TRUE)
refl.dat <- data.table(wavelength = 400:2500, means=refl.means)
refl.dat <- cbind(refl.dat, t(refl.q))
setnames(refl.dat, rownames(refl.q), c("low2", "low1", "high1", "high2"))

## Reflectance error plot
require(ggplot2)
refl.spec <- ggplot(refl.dat) +
    aes(x=wavelength) +
    geom_line(aes(y=means)) +
    geom_ribbon(aes(ymin=low1,ymax=high1), alpha=0.3, 
                fill="red", color="black", linetype="dotted", size=0.2) +
    geom_ribbon(aes(ymin=low2,ymax=high2), alpha=0.3, 
                color="black", linetype="dashed", size=0.2) +
    xlab("Wavelength(nm)") +
    ylab("Mean reflectance error") +
    theme_bw() +
    theme(axis.title = element_text(size=10),
          axis.text = element_text(size=7))
png("manuscript/figures/reflectance-validation.png", 
    width=4, height=3, units="in", res=300)
plot(refl.spec)
dev.off()

