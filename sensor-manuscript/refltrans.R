### Validation based on spectra comparison
# Libraries and sources {{{
library(PEcAnRTM)
library(ggplot2)
library(gridExtra)
library(grid)
# }}}

# Load data {{{ 
# Load FFT data
load("../data/FFT.processed.RData")
fft.f <- fft.f[sensor=="identity"][!is.na(N.mu)][plant.type %in% c("hardwood", "conifer")]
rm(fft.h, fft.c)

# Load transmittance data
keep.wl <- sprintf("Wave_%d", 400:2500)
keep.other <- c("Spectra", "Species", "Sample_Name", "Sample_Year")
keep.cols <- c(keep.wl, keep.other)
dropbox.path <- "~/Dropbox"
trans.path <- file.path(dropbox.path,
                        "NASA_TE_PEcAn-RTM_Project",
                        "Data", "Spectra",
                        "NASA_FFT_Spectra",
                        "NASA_FFT_IS_Tran_Spectra_v4.csv")
trans <- fread(trans.path, header=TRUE, select=keep.cols)
remove.negatives <- function(x){
    if(is.numeric(x)) x[x < 0] <- NA
    return(x)
}
trans <- trans[, lapply(.SD, remove.negatives)]

# Load reflectance data
refl.path <- file.path(dropbox.path,
                        "NASA_TE_PEcAn-RTM_Project",
                        "Data", "Spectra",
                        "NASA_FFT_Spectra",
                        "NASA_FFT_LC_Refl_Spectra_v4.csv")
refl <- fread(refl.path, header=TRUE, select=keep.cols)
# }}}

# Error matrix function {{{
error.matrix <- function(dat.mod, dat.obs, refltrans=2){
    setkey(dat.mod, Sample_Name, Sample_Year)
    setkey(dat.obs, Sample_Name, Sample_Year)
    rt.big <- dat.mod[dat.obs][!is.na(N.mu)]
    rt.mat <- as.matrix(rt.big[,c("N.mu", "Cab.mu", "Car.mu", "Cw.mu", "Cm.mu", keep.wl), with=FALSE])
    tdiff <- function(mrow){
        pars <- mrow[1:5]
        obs <- mrow[-5:0]
        if(any(is.na(pars))) return(rep(NA,2101))
        tmodel <- prospect(pars,5)[,refltrans]
        tdiff <- tmodel - obs
        return(tdiff)
    }
    error.mat <- apply(rt.mat, 1, tdiff)
    #error.mat <- error.mat[!is.na(error.mat[,1]),]
    return(error.mat)
}
# }}}

# Calculate reflectance and transmittance error {{{
rem.all.raw <- error.matrix(fft.f, refl, 1)
rem.h.raw <- error.matrix(fft.f[plant.type=="hardwood"], refl, 1)
rem.c.raw <- error.matrix(fft.f[plant.type=="conifer"], refl, 1)
tem.all.raw <- error.matrix(fft.f, trans, 2)
tem.h.raw <- error.matrix(fft.f[plant.type=="hardwood"], trans, 2)
tem.c.raw <- error.matrix(fft.f[plant.type=="conifer"], trans, 2)
# }}}

# Error plot function {{{
error.plot <- function(err.mat){
    require(ggplot2)
    err.mat <- t(err.mat)
    err.means <- colMeans(err.mat, na.rm=TRUE)
    err.q <- apply(err.mat, 2, quantile, c(0.025,0.10,0.90,0.975), na.rm=TRUE)
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
}
# }}}

# Initialize error plots {{{
re.all.raw <- error.plot(rem.all.raw)
re.h.raw <- error.plot(rem.h.raw)
re.c.raw <- error.plot(rem.c.raw)
te.all.raw <- error.plot(tem.all.raw)
te.h.raw <- error.plot(tem.h.raw)
te.c.raw <- error.plot(tem.c.raw)
# }}}

# Modify and draw plots {{{
th.all <- theme_bw() +
    theme(text = element_text(size=11),
          axis.text = element_text(size=rel(0.6)),
          axis.title = element_text(size=rel(0.8)),
          plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "lines"))
th.mr <- theme(axis.title.y = element_blank()) 
re.ylims <- ylim(-0.04, 0.04)
te.ylims <- ylim(-0.28, 0.14)
re.all <- re.all.raw + ggtitle("All") + ylab("Reflectance error (Model - Obs.)") + th.all + re.ylims
re.h <- re.h.raw + ggtitle("Hardwood") + th.all + th.mr + re.ylims
re.c <- re.c.raw + ggtitle("Conifer") + th.all + th.mr + re.ylims
te.all <- te.all.raw + ylab("Transmittance error (Model - Obs.)") + th.all + te.ylims
te.h <- te.h.raw + th.all + th.mr + te.ylims
te.c <- te.c.raw + th.all + th.mr + te.ylims
png("manuscript/figures/refltrans-validation.png", width=6, height=6, units="in", res=300)
grid.arrange(re.all, re.h, re.c, te.all, te.h, te.c, nrow=2)
dev.off()
# }}} 

# Statistics function {{{
vis.wl <- 400:800
nir.wl <- 801:2500
rmse.wl <- function(mat){
    mat <- t(mat)
    i.vis <- vis.wl - 399
    i.nir <- nir.wl - 399
    rmse <- apply(mat, 2, function(xd) sqrt(mean(xd^2, na.rm=TRUE)))
    rmse.vis <- mean(rmse[i.vis], na.rm=TRUE)
    rmse.nir <- mean(rmse[i.nir], na.rm=TRUE)
    bias <- apply(mat, 2, mean, na.rm=TRUE)
    bias.vis <- mean(bias[i.vis], na.rm=TRUE)
    bias.nir <- mean(bias[i.nir], na.rm=TRUE)
    err.c <- t(apply(mat, 1, "+", -bias))
    sepc <- apply(err.c, 2, function(x) sqrt(mean(x^2, na.rm=TRUE)))
    sepc.vis <- mean(sepc[i.vis], na.rm=TRUE)
    sepc.nir <- mean(sepc[i.nir], na.rm=TRUE)
    out <- c("rmse.vis"=rmse.vis,
             "bias.vis"=bias.vis,
             "sepc.vis"=sepc.vis,
             "rmse.nir"=rmse.nir,
             "bias.nir"=bias.nir,
             "sepc.nir"=sepc.nir)
    return(out)
}
# }}}

# Apply statistics function {{{
rtlist <- list(rem.all.raw, rem.h.raw, rem.c.raw,
               tem.all.raw, tem.h.raw, tem.c.raw)
cnames <- c("All-R", "Hardwood-R", "Conifer-R",
            "All-T", "Hardwood-T", "Conifer-T")
names(rtlist) <- cnames
rmse.full<- lapply(rtlist, rmse.wl)
sumtab <- do.call(cbind, rmse.full)
rownames(sumtab) <- c("VIS-RMSE", "VIS-BIAS", "VIS-SEPC",
                      "IR-RMSE", "IR-BIAS", "IR-SEPC")
tsumtab <- t(sumtab)
print(tsumtab)
# }}}

# Print summary table {{{
library(xtable)
cap <- "
Modeled reflectance (R) and transmittance (T) error statistics.
"
cap <- gsub("\\n", " ", cap)
out.tab <- xtable(tsumtab, caption=cap, label="tab:refltrans", digits=4)
# Post processing
out.tab.pre <- print(out.tab, file="", include.rownames=TRUE)
out.tab.post <- out.tab.pre 
cat(out.tab.post, file="manuscript/tables/error-stats.tex")
# }}}

# vim: set foldlevel=0 :
