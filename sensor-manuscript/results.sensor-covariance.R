# Sensor covariance plot
library(PEcAnRTM)
library(MASS)

datnum <- 77
sensor.sub <- c("identity", "aviris.ng", "landsat8", "modis")
# {{{ Load data
# Data identifier -- Each ID corresponds to a single set of true values
datlist <- list()
path <- "some-simulations"
for(s in sensor.sub){
    f.name <- sprintf("%s.%d.RData", s, datnum)
    f.list <- load.from.name(f.name, path)
    sensor <- list(sensor = f.list$sensor)
    samples <- f.list$samples
    samples.sub <- burnin.thin(samples)
    samples.pairs <- samples.sub[,1:5]
    colnames(samples.pairs) <- params.prospect5
    datlist[[s]] <- samples.pairs
}
true.param <- unlist(f.list[params.prospect5])
samples.all <- do.call(rbind, datlist)
# }}}

quant.val <- 97
# Identity, AVIRIS NG, MODIS, Landsat
pal.sensor <- c("Greens", "Reds", "Blues", "Oranges")
col.sensor <- c("orange", "dark blue", "dark red", "dark green")

# {{{ Density plot functions
quant <- c((100-quant.val)/200, (quant.val+(100-quant.val)/2)/100)
# Generic density plot
densplot <- function(x, y, pal, color, lwd=1){
    z <- kde2d(x, y, n=50)
    nlev <- 7
    my.cols <- brewer.pal(nlev, pal)
    zlim <- range(z$z, finite=TRUE)
    levs <- seq(zlim[1], zlim[2], length.out=nlev)
    #.filled.contour(z$x, z$y, z$z, levs, my.cols)
    contour(z, drawlabels=FALSE, levels=levs, add=TRUE, col=color, lwd=lwd)
}

# Density plot for just identity
identity.plot <- function(nx, ny){
    dat <- samples.all[1:5001,]
    xrange <- quantile(dat[,nx], quant)
    #xrange[1] <- min(xrange[1], true.param[nx])
    #xrange[2] <- max(xrange[2], true.param[nx])
    yrange <- quantile(dat[,ny], quant)
    #yrange[1] <- min(yrange[1], true.param[nx])
    #yrange[2] <- max(yrange[2], true.param[nx])
    x <- dat[, nx]
    y <- dat[, ny]
    plot(xrange, yrange, type='n', xaxt='n', yaxt='n')
    if(ny == 1) axis(3, at=pretty(xrange), labels=pretty(xrange))
    if(nx == 5) axis(4, at=pretty(yrange), labels=pretty(yrange))
    densplot(x, y, "Oranges", "black")
    abline(v=true.param[nx], lwd=1.5, lty=2)
    abline(h=true.param[ny], lwd=1.5, lty=2)
}

# Density plot for all sensors
sensor.plot <- function(nx, ny){
    xrange <- quantile(samples.all[,nx], quant)
    yrange <- quantile(samples.all[,ny], quant)
    plot(xrange, yrange, type='n', xaxt='n', yaxt='n')
    if(ny == 5) axis(1, at=pretty(xrange), labels=pretty(xrange))
    if(nx == 1) axis(2, at=pretty(yrange), labels=pretty(yrange))
    for(i in 4:1){
        imin <- (i-1)*5001 + 1
        imax <- i * 5001
        x <- samples.all[imin:imax,nx]
        y <- samples.all[imin:imax,ny]
        densplot(x, y, pal.sensor[i], col.sensor[i], lwd=i/2)
    }
    abline(v=true.param[nx], lwd=1.5, lty=2)
    abline(h=true.param[ny], lwd=1.5, lty=2)
}
# }}}        

# {{{ Draw density plots
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
        sensor.plot(p[1], p[2])
    }
}
dev.off()
# }}}

