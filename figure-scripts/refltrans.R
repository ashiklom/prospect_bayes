### Validation based on spectra comparison
# Libraries and sources {{{
library(PEcAnRTM)
library(ggplot2)
library(gridExtra)

if(grepl("figure-scripts", getwd())) setwd("..")
source("preprocess.fft.R")
source("figure-scripts/common.R")
# }}}

# Error matrix function {{{
error.matrix <- function(dat.mod, dat.obs, refltrans=2, relative=FALSE){
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
        if(relative) tdiff <- tdiff / tmodel
        return(tdiff)
    }
    rt.mat <- as.matrix(rt.big[,c("N.mu", "Cab.mu", "Cw.mu", "Cm.mu", wl.names), with=FALSE])
    trans.mat <- matrix(NA, nrow(rt.big), 2101)
    for(i in 1:nrow(rt.big)){
        trans.mat[i,] <- tdiff(rt.mat[i,])
    }
    trans.mat <- trans.mat[!is.na(trans.mat[,1]),]
    return(trans.mat)
}
# }}}
# Error plot function {{{
error.plot <- function(err.mat){
    #err.mat <- error.matrix(dat.mod, dat.obs, refltrans)
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
}
# }}}

# Load data {{{ 
#' Load transmittance data
trans.path <- file.path("~", "Documents", "Dropbox",
                        "NASA_TE_PEcAn-RTM_Project",
                        "Data", "Spectra",
                        "NASA_FFT_Spectra",
                        "NASA_FFT_IS_Tran_Spectra_v4.csv")
trans <- fread(trans.path, header=TRUE)
remove.negatives <- function(x){
    if(is.numeric(x)) x[x < 0] <- NA
    return(x)
}
trans <- trans[, lapply(.SD, remove.negatives)]
# Remove transmittance data with negative values in visible
transmat <- as.matrix(trans[,-(1:71), with=F])
bad.trans <- which(apply(trans.mat[,1:1901], 1, function(x) any(x < 0)))
trans <- trans[-bad.trans,]

# Load reflectance data
refl.path <- file.path("~", "Documents", "Dropbox",
                        "NASA_TE_PEcAn-RTM_Project",
                        "Data", "Spectra",
                        "NASA_FFT_Spectra",
                        "NASA_FFT_LC_Refl_Spectra_v4.csv")
refl <- fread(refl.path, header=TRUE)
# }}}

# Calculate reflectance and transmittance error {{{ 
#fft <- fft[plant.type %in% c("hardwood", "conifer"),]
#rem.all.raw <- error.matrix(fft, refl, 1)
#rem.h.raw <- error.matrix(fft.h, refl, 1)
#rem.c.raw <- error.matrix(fft.c, refl, 1)
#tem.all.raw <- error.matrix(fft, trans, 2)
#tem.h.raw <- error.matrix(fft.h, trans, 2)
#tem.c.raw <- error.matrix(fft.c, trans, 2)
#remr.all.raw <- error.matrix(fft, refl, 1, relative=TRUE)
#remr.h.raw <- error.matrix(fft.h, refl, 1, relative=TRUE)
#remr.c.raw <- error.matrix(fft.c, refl, 1, relative=TRUE)
#temr.all.raw <- error.matrix(fft, trans, 2, relative=TRUE)
#temr.h.raw <- error.matrix(fft.h, trans, 2, relative=TRUE)
#temr.c.raw <- error.matrix(fft.c, trans, 2, relative=TRUE)
#save(rem.all.raw, rem.h.raw, rem.c.raw,
     #tem.all.raw, tem.h.raw, tem.c.raw,
     #remr.all.raw, remr.h.raw, remr.c.raw,
     #temr.all.raw, temr.h.raw, temr.c.raw,
     #file="data/refltrans-error.RData")

load("data/refltrans-error.RData")
re.all.raw <- error.plot(rem.all.raw)
re.h.raw <- error.plot(rem.h.raw)
re.c.raw <- error.plot(rem.c.raw)
te.all.raw <- error.plot(tem.all.raw)
te.h.raw <- error.plot(tem.h.raw)
te.c.raw <- error.plot(tem.c.raw)
# }}}

# Modify plots {{{
th.all <- theme_bw() +
    theme(text = element_text(size=11),
          axis.text = element_text(size=rel(0.6)),
          axis.title = element_text(size=rel(0.8)),
          plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "lines"))
th.mr <- theme(axis.title.y = element_blank()) 
re.ylims <- ylim(-0.05, 0.12)
te.ylims <- ylim(-0.28, 0.16)
re.all <- re.all.raw + ggtitle("All") + ylab("Reflectance error (Model - Obs.)") + th.all + re.ylims
re.h <- re.h.raw + ggtitle("Hardwood") + th.all + th.mr + re.ylims
re.c <- re.c.raw + ggtitle("Conifer") + th.all + th.mr + re.ylims
te.all <- te.all.raw + ylab("Transmittance error (Model - Obs.)") + th.all + te.ylims
te.h <- te.h.raw + th.all + th.mr + te.ylims
te.c <- te.c.raw + th.all + th.mr + te.ylims
png.plot("manuscript/figures/refltrans-validation.png", w=6, h=6)
grid.arrange(re.all, re.h, re.c, te.all, te.h, te.c, nrow=2)
dev.off()
# }}} 

# Calculate comparison statistics {{{
rmse.wl <- function(xd){
    nx <- sum(!is.na(xd))
    rmse <- sqrt(sum(xd^2, na.rm=TRUE)/nx)
    bias <- mean(xd, na.rm=TRUE)
    sepc <- sqrt(sum((xd-bias)^2, na.rm=TRUE)/nx)
    cv <- sepc / mean(xd, na.rm=TRUE)
    out <- c("rmse"=rmse, "bias"=bias, "sepc"=sepc, "cv"=cv)
    return(out)
}
rtlist <- list(rem.all.raw, rem.h.raw, rem.c.raw,
               tem.all.raw, tem.h.raw, tem.c.raw)
cnames <- c("All-R", "Hardwood-R", "Conifer-R",
            "All-T", "Hardwood-T", "Conifer-T")
names(rtlist) <- cnames
rmse.full<- lapply(rtlist, apply, 2, rmse.wl)
rmse.vis <- lapply(rmse.full, function(x) rowMeans(x[,(400:800)-399]))
rmse.ir <- lapply(rmse.full, function(x) rowMeans(x[,(801:2500)-399]))
sumtab <- do.call(cbind, c(rmse.vis, rmse.ir))
sumtab.names <- sprintf("%s-%s", rep(cnames,2), rep(c("VIS","IR"), each=length(cnames)))
colnames(sumtab) <- sumtab.names
rmse.simple <- function(x) sqrt(sum(x^2, na.rm=TRUE) / sum(!is.na(x)))
rerlist <- list(remr.all.raw, remr.h.raw, remr.c.raw,
                temr.all.raw, temr.h.raw, temr.c.raw)
names(rerlist) <- cnames
rmspe.full <- lapply(rerlist, apply, 2, rmse.simple)
rmspe.vis <- lapply(rmspe.full, function(x) mean(x[(400:800)-399]))
rmspe.ir <- lapply(rmspe.full, function(x) mean(x[(800:2300)-399]))
rmspe <- c(rmspe.vis, rmspe.ir)
sumtab <- rbind(sumtab, rmspe)
statnames <- c("RMSE", "BIAS", "SEPC", "CV", "RMSPE")
rownames(sumtab) <- statnames
print(sumtab)
# }}}

# Print summary table {{{
library(xtable)
cap <- "
Modeled reflectance (R) and transmittance (T) error statistics.
"
cap <- gsub("\\n", " ", cap)
tsumtab <- t(sumtab)
out.tab <- xtable(tsumtab, caption=cap, label="tab:error-stats",
                  digits=-3)
# Post processing
out.tab.pre <- print(out.tab, file="", include.rownames=TRUE)
out.tab.post <- out.tab.pre 
cat(out.tab.post, file="manuscript/tables/error-stats.tex")
# }}}

# vim: set foldlevel=0 :
