### Validation based on spectra comparison
library(PEcAnRTM)
library(ggplot2)
library(gridExtra)

if(grepl("figure-scripts", getwd())) setwd("..")
source("preprocess.fft.R")
source("figure-scripts/common.R")

error.matrix <- function(dat.mod, dat.obs, refltrans=2){# [--
    wl.names <- sprintf("Wave_%d", 400:2500)
    d.names <- names(dat.mod)
    d.names <- d.names[!(grepl("Wave_", d.names))]
    d.names[d.names %in% names(dat.mod)]
    setkey(dat.mod, Sample_Name, Sample_Year)
    setkey(dat.obs, Sample_Name, Sample_Year)
    rt.big <- dat.mod[dat.obs]
    tdiff <- function(mrow){
        pars <- mrow[1:4]
        obs <- mrow[-4:0]
        if(any(is.na(pars))) return(rep(NA,2101))
        tmodel <- prospect(pars,4)[,refltrans]
        tdiff <- tmodel - obs
        return(tdiff)
    }
    rt.mat <- as.matrix(rt.big[,c("N.mu", "Cab.mu", "Cw.mu", "Cm.mu", wl.names), with=FALSE])
    trans.mat <- matrix(NA, nrow(rt.big), 2101)
    for(i in 1:nrow(rt.big)){
        trans.mat[i,] <- tdiff(rt.mat[i,])
    }
    return(trans.mat)
}# --]

error.plot <- function(dat.mod, dat.obs, refltrans){# [--
    err.mat <- error.matrix(dat.mod, dat.obs, refltrans)
    require(ggplot2)
    err.means <- colMeans(err.mat, na.rm=TRUE)
    err.q <- apply(err.mat, 2, quantile, c(0.025,0.1,0.9,0.975), na.rm=TRUE)
    err.dat <- data.table(wavelength = 400:2500, means=err.means)
    err.dat <- cbind(err.dat, t(err.q))
    setnames(err.dat, rownames(err.q), c("low2", "low1", "high1", "high2"))
    err.spec <- ggplot(err.dat) +
        aes(x=wavelength) +
        geom_line(aes(y=means)) +
        geom_ribbon(aes(ymin=low1,ymax=high1), alpha=0.3, 
                    fill="red", color="black", linetype="dotted", size=0.2) +
        geom_ribbon(aes(ymin=low2,ymax=high2), alpha=0.3, 
                    color="black", linetype="dashed", size=0.2) +
        xlab("Wavelength(nm)") +
        theme_bw() +
        theme(axis.title = element_text(size=10),
            axis.text = element_text(size=7))
}# --]

#' Load transmittance data
trans.path <- file.path("~", "Documents", "Dropbox",
                        "NASA_TE_PEcAn-RTM_Project",
                        "Data", "Spectra",
                        "NASA_FFT_Spectra",
                        "NASA_FFT_IS_Tran_Spectra_v4.csv")
trans <- fread(trans.path, header=TRUE)

# Load reflectance data
refl.path <- file.path("~", "Documents", "Dropbox",
                        "NASA_TE_PEcAn-RTM_Project",
                        "Data", "Spectra",
                        "NASA_FFT_Spectra",
                        "NASA_FFT_LC_Refl_Spectra_v4.csv")
refl <- fread(refl.path, header=TRUE)

#' Calculate reflectance error
th.mr <- theme_bw() + theme(axis.title.y = element_blank())
re.all <- error.plot(fft, refl, 1) + ggtitle("All") + ylab("Reflectance error")
re.h <- error.plot(fft.h, refl, 1) + ggtitle("Hardwood") + th.mr
re.c <- error.plot(fft.c, refl, 1) + ggtitle("Conifer") + th.mr

#' Calculate transmittance error
te.all <- error.plot(fft, trans, 2) + ggtitle("All") + ylab("Transmittance error")
te.h <- error.plot(fft.h, trans, 2) + ggtitle("Hardwood") + th.mr
te.c <- error.plot(fft.c, trans, 2) + ggtitle("Conifer") + th.mr

png.plot("manuscript/figures/refltrans-validation.png", w=6, h=6)
grid.arrange(re.all, re.h, re.c, te.all, te.h, te.c, nrow=2)
dev.off()

