#' ---
#' title: Validation based on spectral comparison
#' author: Alexey Shiklomanov
#' output_format: html_document
#' ---

#' # Introduction
#' One way to validate a model inversion is to compare the model output in 
#' forward mode to the original measurements used to perform the inversion.  
#' This script does just that by running PROSPECT 5 with parameters from the 
#' FFT database inversion and comparing the resultign simulated reflectance and 
#' transmittance with the original measurements. The results of this validation 
#' are presented as a figure displaying the 90th and 95th percentiles of the 
#' error (bias) of the simulated spectra as well as a table displaying relevant 
#' summary statistics (see manuscript text for more details.)

#' # Setup
#' First, we load required packages. `PEcAnRTM` is used to perform simulations 
#' and contains useful information vectors. `ggplot2`, `gridExtra`, and `grid` 
#' are used to draw and arrange the figure. `xtable` is used to generate a 
#' pretty LaTeX table from the R output.
library(PEcAnRTM)
library(ggplot2)
library(gridExtra)
library(grid)
library(xtable)

#' The data for generating simulated spectra come from the `fft.f` object from 
#' the `FFT.processed.RData` file (for details, see `fft.load.R` and 
#' `fft.preprocess.R`). We then subset `fft.f` to only the full spectra, where 
#' inversion parameter estimates are not `NA`, and only for tree species (i.e.  
#' excluding shrubs and grasses).

load("../data/FFT.processed.RData")
fft.f <- fft.f[sensor=="identity"][!is.na(N.mu)][plant.type %in% c("hardwood", "conifer")]
rm(fft.h, fft.c)

#' We then load the observed transmittance data. These data are stored in a 
#' private Dropbox and are available on request. To accelerate and facilitate 
#' loading and pre-processing, we only load a subset of the columns, indicated 
#' by `keep.cols`.

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

#' Because they are a physical impossibility and are indicative of measurement 
#' error, we convert all negative transmittance values to `NA`.

remove.negatives <- function(x){
    if(is.numeric(x)) x[x < 0] <- NA
    return(x)
}
trans <- trans[, lapply(.SD, remove.negatives)]

#' We then load the reflectance data, stored in the same place.

refl.path <- file.path(dropbox.path,
                        "NASA_TE_PEcAn-RTM_Project",
                        "Data", "Spectra",
                        "NASA_FFT_Spectra",
                        "NASA_FFT_LC_Refl_Spectra_v4.csv")
refl <- fread(refl.path, header=TRUE, select=keep.cols)

#' # Calculating the error
#' We define a function that takes data tables of parameter inversion estimates 
#' (`dat.mod`) and measured spectra (`dat.obs`) as input and returns a matrix 
#' of error values by wavelength. The argument `refltrans` dictates whether 
#' this is comparing reflectance (1) or transmittance (2) by selecting the 
#' appropriate column from the PROSPECT model output.  First, this function 
#' merges the two input data tables into a single large table (`rt.big`) with 
#' correctly aligned rows based on the sample name and year. We then subset 
#' `rt.big` to only the mean parameter estimates and the measured spectra and 
#' convert the resulting object to a matrix (`rt.mat`). We define a `tdiff` 
#' function that operates on a single row of this matrix, simulating the 
#' spectra based on the PROSPECT parameters (first 5 columns) and subtracting 
#' this from the spectra (remaining columns).  Finally, we apply this function 
#' to the matrix row-by-row, creating a matrix of error values.

error.matrix <- function(dat.mod, dat.obs, refltrans=2){
    setkey(dat.mod, Sample_Name, Sample_Year)
    setkey(dat.obs, Sample_Name, Sample_Year)
    rt.big <- dat.mod[dat.obs][!is.na(N.mu)]
    parnames <- sprintf("%s.mu", params.prospect5)
    rt.mat <- as.matrix(rt.big[,c(parnames, keep.wl), with=FALSE])
    tdiff <- function(mrow){
        pars <- mrow[1:5]
        obs <- mrow[-5:0]
        if(any(is.na(pars))) return(rep(NA,2101))
        tmodel <- prospect(pars,5)[,refltrans]
        tdiff <- tmodel - obs
        return(tdiff)
    }
    error.mat <- apply(rt.mat, 1, tdiff)
    return(error.mat)
}

#' With the `error.matrix` function defined, we apply it to reflectance and 
#' transmittance data. For each, we compute the error for all available data 
#' and again for hardwood and conifer species separately, due to the large 
#' systematic differences in both measurements and inversion parameter 
#' estimates for both. Here and consequently, we use the `re` prefix to 
#' indicate reflectance and `te` to indicate transmittance.

rem.all.raw <- error.matrix(fft.f, refl, 1)
rem.h.raw <- error.matrix(fft.f[plant.type=="hardwood"], refl, 1)
rem.c.raw <- error.matrix(fft.f[plant.type=="conifer"], refl, 1)
tem.all.raw <- error.matrix(fft.f, trans, 2)
tem.h.raw <- error.matrix(fft.f[plant.type=="hardwood"], trans, 2)
tem.c.raw <- error.matrix(fft.f[plant.type=="conifer"], trans, 2)

#' # Plotting the error
#' Here, we define a function that takes an error matrix as input and produces 
#' a nicely formatted `ggplot` graphic as output. The function first calculates 
#' the mean and 80th and 95th percentiles of the error by wavelength. It then 
#' draws the mean as a solid black line and the percentile regions as lightly 
#' shaded regions bounded by dashed lines. We then apply this function to each 
#' of the error matrices computed above.

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

re.all.raw <- error.plot(rem.all.raw)
re.h.raw <- error.plot(rem.h.raw)
re.c.raw <- error.plot(rem.c.raw)
te.all.raw <- error.plot(tem.all.raw)
te.h.raw <- error.plot(tem.h.raw)
te.c.raw <- error.plot(tem.c.raw)

#' With all six plots generated, we then post-process these plots to create a 
#' single multi-panel figure. To emphasize differences between hardwoods and 
#' conifers, we set common y-axis limits for reflectance and again for 
#' transmittance. We also eliminate redundant axes where appropriate.

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
pdf("manuscript/figures/refltrans-validation.pdf", width=6, height=6)
grid.arrange(re.all, re.h, re.c, te.all, te.h, te.c, nrow=2)
dev.off()

#' # Error statistics table
#' To facilitate comparison with other studies and due to their inherent 
#' differences, we compute statistics separately for the visible (VIS) and near 
#' infrared (NIR) regions of the spectrum, defined here.

vis.wl <- 400:800
nir.wl <- 801:2500

#' We define a function that takes a matrix as input and returns a vector 
#' containing the values of error statistics as output. For each statistic, we 
#' first compute the statistic for each wavelength and then calculate a mean 
#' according the spectral region of interest.

rmse.wl <- function(mat){
    mat <- t(mat)
    i.vis <- vis.wl - 399
    i.nir <- nir.wl - 399
    rmse <- apply(mat, 2, function(x) sqrt(mean(x^2, na.rm=TRUE)))
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

#' We then apply the above function to each error matrix. For conciseness and 
#' to facilitate merging the results into a single data table, we first create 
#' a list containing all the error matrices, add informative names to the list, 
#' and then use `lapply` to apply the error statistic function to each item of 
#' the list. We then `cbind` the resulting list of vectors into a matrix.

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

#' Finally, we use `xtable` combined with some post-processing to generate a 
#' pretty LaTeX table of the results. Values are printed to 4 decimal places, 
#' and the `centerline` TeX environment is used to center the table across the 
#' entire page.

cap <- "
Modeled reflectance (R) and transmittance (T) error statistics.
"
cap <- gsub("\\n", " ", cap)
out.tab <- xtable(tsumtab, caption=cap, label="tab:refltrans", digits=4)
# Post processing
out.tab.pre <- print(out.tab, file="", include.rownames=TRUE)
out.tab.post <- out.tab.pre
out.tab.post <- gsub("centering", "centerline{", out.tab.post)
out.tab.post <- gsub("(end\\{tabular\\})", "\\1\n\\}", out.tab.post)
cat(out.tab.post, file="manuscript/tables/refltrans.tex")
