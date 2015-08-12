# Test reflectance and transmittance error for simulated spectra
library(PEcAnRTM)

nl <- 5
N.seq <- seq(1, 4, length.out=nl)
Cab.seq <- seq(10, 100, length.out=nl)
Car.seq <- seq(2, 60, length.out=nl)
Cw.seq <- c(seq(0.004, 0.01, length.out=nl/2), 
            seq(0.015, 0.05, length.out=nl/2))
Cm.seq <- c(seq(0.001, 0.005, length.out=nl/2),
            seq(0.006, 0.015, length.out=nl/2))
seq.l <- list(N=N.seq, Cab=Cab.seq, Car=Car.seq, Cw=Cw.seq, Cm=Cm.seq)

defpars <- c("N"=1.4, "Cab"=40, "Car"=8, "Cw"=0.01, "Cm"=0.005)

mat <- matrix(NA, 2101, nl)
mat.lr <- list(N=mat, Cab=mat, Car=mat, Cw=mat, Cm=mat)
mat.lt <- mat.lr
mat.names <- names(mat.lr)
for(n in mat.names){
    colnames(mat.lr[[n]]) <- seq.l[[n]]
    colnames(mat.lt[[n]]) <- seq.l[[n]]
    for(i in 1:nl){
        params <- defpars
        val <- seq.l[[n]][i]
        params[n] <- val
        RT <- prospect(params, 5)
        mat.lr[[n]][,i] <- RT[,1]
        mat.lt[[n]][,i] <- RT[,2]
    }
}

mat.r <- do.call(cbind, mat.lr)
mat.t <- do.call(cbind, mat.lt)
