# Analysis for ED2 optical parameters
# Load packages {{{
library(data.table)
library(ggplot2)
library(reshape2)
# }}}

# Load data {{{
fft.path <- file.path("~", "Dropbox", "NASA_TE_PEcAn-RTM_Project",
                      "Data", "Spectra", "NASA_FFT_Spectra")
fft.path.refl <- file.path(fft.path, "NASA_FFT_LC_Refl_Spectra_v4.csv")
fft.path.trans <- file.path(fft.path, "NASA_FFT_IS_Tran_Spectra_v4.csv")
fft.path.info <- file.path("~/Dropbox/NASA_TE_PEcAn-RTM_Project", "Data",
                           "FFT_species_info_csv.csv")
fft.refl.raw <- fread(fft.path.refl, header=TRUE)
fft.trans.raw <- fread(fft.path.trans, header=TRUE)
fft.info <- fread(fft.path.info, header=TRUE)
setkey(fft.refl.raw, "Species")
setkey(fft.trans.raw, "Species")
setkey(fft.info, "Label")
fft.refl <- fft.refl.raw[fft.info]
fft.trans <- fft.trans.raw[fft.info, allow.cartesian=TRUE]
# }}}

# Function to isolate VIS and NIR {{{
info <- c("Spectra", "Species", "Height", "Age", "Height_Age",
          "Sample_Name", "Sample_Year")
melt.spec <- function(dat, wl.full=400:2500){
    wl <- sprintf("Wave_%d", wl.full)
    dat <- dat[, c(info, wl), with=FALSE]
    dat.melt <- melt(dat, measure.vars=wl, id.vars=info, variable.name="wavelength")
    numeric.wavelength <- function(x) as.numeric(gsub("Wave_","",x))
    dat.melt[, wavelength := numeric.wavelength(wavelength)]
    dat.melt[value < 0, value := 0]
    return(dat.melt)
}
# }}}

fft.refl.melt <- melt.spec(fft.refl)
fft.trans.melt <- melt.spec(fft.trans)

fft.refl.melt[wavelength %in% vis, nv := "PAR"]
fft.refl.melt[wavelength %in% nir, nv := "NIR"]
fft.trans.melt[wavelength %in% vis, nv := "PAR"]
fft.trans.melt[wavelength %in% nir, nv := "NIR"]

refl.vn <- fft.refl.melt[, mean(value), by = c("nv", info)]

nir.break <- 680
vis <- 400:(nir.break-1)
nir <- nir.break:2500
refl.sum <- ggplot(fft.refl.melt, aes(x=wavelength, y=value)) +
    stat_summary(fun.y="mean", geom="line") + 
    stat_summary(fun.ymax=function(x) quantile(x, 0.975), 
                 fun.ymin=function(x) quantile(x, 0.025),
                 geom="ribbon", color="grey", alpha=0.4) +
    annotate("rect", xmin=min(vis), xmax=max(vis), ymin=0, ymax=0.05,
             fill="green", alpha=0.2) +
    annotate("rect", xmin=min(nir), xmax=max(nir), ymin=0, ymax=0.05,
             fill="red", alpha=0.2) +
    ylab("Reflectance")
trans.sum <- refl.sum %+% fft.trans.melt + ylab("Transmittance")
require(gridExtra)
png("other-figures/summary-spectrum.png")
grid.arrange(refl.sum, trans.sum, nrow=2)
dev.off()

plot.refl <- ggplot(fft.refl.melt, aes(x=nv, y=value, fill=nv)) +
    geom_violin() + coord_flip()
plot.trans <- plot.refl %+% fft.trans.melt
png("other-figures/violin-nv.png")
grid.arrange(plot.refl, plot.trans, nrow=2)
dev.off()

summary.stats <- function(x)
    list(mean = mean(x), sd = sd(x),
         min = min(x), med = median(x),
         max = max(x), q2.5 = quantile(x, 0.025),
         q97.5 = quantile(x, 0.975))
refl.summary <- fft.refl.melt[, summary.stats(value), by=nv]
trans.summary <- fft.trans.melt[, summary.stats(value), by=nv]
