if(grepl("figure-scripts", getwd())) setwd("..")

library(PEcAnRTM)
library(colorRamps)
source("figure-scripts/common.R")

nl <- 20
N.seq <- seq(1, 4, length.out=nl)
Cab.seq <- seq(10, 100, length.out=nl)
Cw.seq <- c(seq(0.004, 0.01, length.out=nl/2), 
            seq(0.015, 0.05, length.out=nl/2))
Cm.seq <- c(seq(0.001, 0.005, length.out=nl/2),
            seq(0.006, 0.015, length.out=nl/2))
seq.l <- list(N=N.seq, Cab=Cab.seq, Cw=Cw.seq, Cm=Cm.seq)

defpars <- c("N"=1.4, "Cab"=40, "Cw"=0.01, "Cm"=0.005)

mat <- matrix(NA, 2101, 20)
mat.lr <- list(N=mat, Cab=mat, Cw=mat, Cm=mat)
mat.lt <- mat.lr
mat.names <- names(mat.l)
for(n in mat.names){
    colnames(mat.lr[[n]]) <- seq.l[[n]]
    colnames(mat.lt[[n]]) <- seq.l[[n]]
    for(i in 1:nl){
        params <- defpars
        val <- seq.l[[n]][i]
        params[n] <- val
        RT <- prospect(params, 4)
        mat.lr[[n]][,i] <- RT[,1]
        mat.lt[[n]][,i] <- RT[,2]
    }
}

lb <- c("N"="(a)", "Cab"="(b)", "Cw"="(c)", "Cm"="(d)")
pal <- blue2green(nl)
png.plot("manuscript/figures/sensitivity.png", h=4, w=4)
par(mfrow=c(2,2), mar=c(2.5, 2.5, 1.5, 0.6), cex=0.7, cex.axis=0.85)
for(n in mat.names) {
    matplot(mat.lr[[n]], type='l', col=pal, lty=1,
            xlab="Reflectance", ylab="Wavelength",
            ylim=c(0, 0.7), mgp=c(1.5,0.6,0))
    text(1800, 0.65, lb[n])
}
dev.off()
