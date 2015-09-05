#' ---
#' title: Sensor covariance plot
#' author: Alexey Shiklomanov
#' output_format: html_document
#' ---

#' # Introduction
#' One of the main goals of this manuscript is to demonstrate how covariance 
#' and uncertainty of RTM parameter retrieval changes with spectral resolution.  
#' This script generates a figure comparing the joint probability distributions 
#' of the original spectra and three of the sensors investigated. The 
#' distributions are shown using contour lines, and each sensor is shown in a 
#' different color.

#' # Load packages and data
#' First, we load the packages. PEcAnRTM is used for data processing functions 
#' and some built-in descriptive character vectors. MASS is used for the 2D 
#' kernel density estimates.

library(PEcAnRTM)
library(MASS)

#' The data for this script are stored with a filename such as 
#' `modis.38.RData`, where `modis` is the sensor tag and `38` is the ID of the 
#' starting parameters used (i.e. `modis.38.RData`, `landsat8.38.RData`, etc.  
#' are all based on the same observation spectrum with the same set of 
#' parameters).  The value `datnum` specifies which data to load. Prior to 
#' running this script, all the required RData files (for all 4 
#' sensors--defined by `sensor.sub`--with the correct ID) must be in the 
#' directory set by `path`.

datnum <- 77
sensor.sub <- c("identity", "aviris.ng", "landsat8", "modis")
path <- "some-simulations"

#' The following loop loads the RData file for each sensor, extracts the matrix 
#' of MCMC samples, performs automatic burn-in and thinning (see PEcAnRTM 
#' documentation), selects only the PROSPECT 5 parameter MCMC samples (omitting 
#' the residual), and saves this matrix to a list. The sensor is extracted as 
#' an integer code corresponding to the `sensor.sub` vector (via R's `factor` 
#' function) and is stored as the sixth column of the samples matrix.

datlist <- list()
for(s in sensor.sub){
    f.name <- sprintf("%s.%d.RData", s, datnum)
    f.list <- load.from.name(f.name, path)
    sensor <- as.integer(factor(s, levels=sensor.sub))
    sensor.factor <- factor(sensor, levels=sensor.sub)
    samples.full <- f.list$samples
    samples.sub <- burnin.thin(samples.full)
    samples.pairs <- samples.sub[,1:5]
    samples <- cbind(samples.pairs, sensor)
    colnames(samples) <- c(params.prospect5, "sensor")
    datlist[[s]] <- samples
}

#' This matrix list is then converted into a single large matrix 
#' (`samples.all`) for subsequent analysis.

samples.all <- do.call(rbind, datlist)

#' For validation, we also extract the true parameter values from the last 
#' loaded object (which one we do doesn't actually matter because these values 
#' are the same within an ID value). 
true.param <- unlist(f.list[params.prospect5])

#' # Generating the plots
#' For convenience, we define the colors for the plots here and reference this 
#' vector in the function.
col.sensor <- c("identity" = "orange",
                "aviris.ng" = "dark blue",
                "landsat8" = "dark red",
                "modis" = "dark green")

#' By default, R's plot function sets axes based on the entire width of the 
#' data. This doesn't actually work very well for contour plots because the 
#' lines don't include all of the data, instead only going up to a certain 
#' quantile. These two facts mean that default axes for contour plots are too 
#' zoomed out and difficult to read. To fix this, we set the axis limits based 
#' on quantiles of the data. This results in the outermost contour lines being 
#' cut off in some places, but the gain in zoom makes the plots easier to 
#' interpret.
quant <- c(0.015, 0.985)

#' The following function generates a generic density plot. It first computes 
#' the 2D kernel density estimate (`kde2d` from the `MASS` package).  It then 
#' partitions this into 7 groups based on the range the density values and 
#' plots the resulting contour lines. Line color and width can be passed as 
#' optional arguments. 
densplot <- function(x, y, color="black", lwd=1){
    z <- kde2d(x, y, n=50)
    nlev <- 7
    zlim <- range(z$z, finite=TRUE)
    levs <- seq(zlim[1], zlim[2], length.out=nlev)
    contour(z, drawlabels=FALSE, levels=levs, add=TRUE, col=color, lwd=lwd)
}

#' To use space in the covariance figure most efficiently, we decided to use 
#' the panel below the diagonal for a comparison of the four sensors and to use 
#' the panel above the diagonal to zoom in on the covariance in the full 
#' spectra (since they are much smaller than the covariances for the other 
#' sensors and therefore hard to pick out on the sensor plot). Such a figure is 
#' impossible to generate using R's built-in `pairs` function, so we had to set 
#' the figure up manually, using column indices in the `samples` matrix. 

#' The following two functions are used to generate each of the lower and upper 
#' panels, respectively. Axes are determined based on the quantiles of the data 
#' (set above). Axis labels are placed only for the outer plots to improve 
#' readability.
#' 
#' For the sensor figure, sensors are plotted in reverse order from 
#' `sensor.sub` such that progressively tighter estimates appear above wider 
#' ones. In addition, the line widths are progressively increased from widest 
#' to narrowest to improve contrast.

sensor.plot <- function(nx, ny){
    xrange <- quantile(samples.all[,nx], quant)
    yrange <- quantile(samples.all[,ny], quant)
    plot(xrange, yrange, type='n', xaxt='n', yaxt='n')
    if(ny == 5) axis(1, at=pretty(xrange), labels=pretty(xrange))
    if(nx == 1) axis(2, at=pretty(yrange), labels=pretty(yrange))
    for(i in 4:1){
        dat <- samples.all[samples.all[,"sensor"] == i,]
        s <- col.sensor[sensor.sub[i]]
        densplot(dat[,nx], dat[,ny], s, lwd=i/2)
    }
    abline(v=true.param[nx], lwd=1.5, lty=2)
    abline(h=true.param[ny], lwd=1.5, lty=2)
}

identity.plot <- function(nx, ny){
    dat <- samples.all[samples.all[,"sensor"] == 1,]
    xrange <- quantile(dat[,nx], quant)
    yrange <- quantile(dat[,ny], quant)
    x <- dat[, nx]
    y <- dat[, ny]
    plot(xrange, yrange, type='n', xaxt='n', yaxt='n')
    if(ny == 1) axis(3, at=pretty(xrange), labels=pretty(xrange))
    if(nx == 5) axis(4, at=pretty(yrange), labels=pretty(yrange))
    densplot(x, y)
    abline(v=true.param[nx], lwd=1.5, lty=2)
    abline(h=true.param[ny], lwd=1.5, lty=2)
}

#' Finally, we draw the panels. The list `plt.list` contains information about 
#' the location and content of each panel in the final figure, with the first 
#' number indicating the variable on the X axis, the second number indicating 
#' the variable on the Y axis, and the third number indicating whether the 
#' panel is an upper panel (1) or lower panel (2). Along the diagonal are the 
#' parameter labels. We then loop over this list and display the content 
#' accordingly.

#+ fig.width=7, fig.height=5
plt.list <- list("N", c(2, 1, 1), c(3, 1, 1), c(4, 1, 1), c(5, 1, 1),
                 c(1, 2, 2), "Cab", c(3, 2, 1), c(4, 2, 1), c(5, 2, 1),
                 c(1, 3, 2), c(2, 3, 2), "Car", c(4, 3, 1), c(5, 3, 1),
                 c(1, 4, 2), c(2, 4, 2), c(3, 4, 2), "Cw", c(5, 4, 1),
                 c(1, 5, 2), c(2, 5, 2), c(3, 5, 2), c(4, 5, 2), "Cm")

pdf(file="manuscript/drive-folder/pairs-4.pdf", width=7, height=5)
par(mfrow=c(5,5), mar=c(1,1,1,1), oma=c(2,2,2,2))
for(p in plt.list){
    if(is.character(p)){
# Diagonal -- write parameter
        plot(c(0,1), c(0,1), ann=F, type='n', xaxt='n', yaxt='n')
        text(x=0.5, y=0.5, p, cex=2)
    } else if(p[3] == 1) {
# Upper panel -- plot full spectra densplot
        identity.plot(p[1], p[2])
    } else if(p[3] == 2) {
# Lower panel -- plot density for each sensor
        sensor.plot(p[1], p[2])
    }
}
par(xpd=TRUE)
legend(-0.09, 0.28, c("Full", "AVIRIS", "Landsat8", "MODIS"), 
       pch=16, col=col.sensor, 
       bty='n', ncol=2, x.intersp=0.5, text.width=0.33)
dev.off()

