## Sensitivity analysis
library(ggplot2)
library(reshape2)
library(PEcAnRTM)

data(prospect4)
wl <- 400:2500
plot(wl, P4data[,1], type='l',
     xlab = "Wavelength",
     ylab = "Relative absorption",
     main = "Chlorophyll absorption features")
plot(wl, P4data[,2], type='l',
     xlab = "Wavelength",
     ylab = "Relative absorption",
     main = "Water absorption features")
plot(wl, P4data[,3], type='l',
     xlab = "Wavelength",
     ylab = "Relative absorption",
     main = "Dry matter absorption features")

ns <- 30
N <- seq(1, 3, length.out = ns)
Ns <- sapply(N, function(x) prospect4(c(x, 30, 0.01, 0.01)))
Cab <- seq(10, 70, length.out = ns)
Cabs <- sapply(Cab, function(x) prospect4(c(1.4, x, 0.01, 0.01)))
Cw <- seq(1e-5, 1e-2, length.out = ns)
Cws <- sapply(Cw, function(x) prospect4(c(1.4, 30, x, 0.01)))
Cm <- seq(1e-6, 1e-2, length.out = ns)
Cms <- sapply(Cm, function(x) prospect4(c(1.4, 30, 0.01, x)))

color <- scale_color_continuous(low="red", high="green")
Nsr <- melt(Ns)
Nsr$Var2 <- rep(N, each=2101)
Np <- ggplot(Nsr) + aes(x=Var1+399, y = value, col=Var2) + geom_line() + 
        xlab('Wavelength') + ylab('Reflectance') + color
Cabsr <- melt(Cabs)
Cabsr$Var2 <- rep(Cab, each=2101)
Cabp <- ggplot(Cabsr) + aes(x=Var1+399, y = value, col=Var2) + geom_line() + 
        xlab('Wavelength') + ylab('Reflectance') + color
Cwsr <- melt(Cws)
Cwsr$Var2 <- rep(Cw, each=2101)
Cwp <- ggplot(Cwsr) + aes(x=Var1+399, y = value, col=Var2) + geom_line() + 
        xlab('Wavelength') + ylab('Reflectance') + color
Cmsr <- melt(Cms)
Cmsr$Var2 <- rep(Cm, each=2101)
Cmp <- ggplot(Cmsr) + aes(x=Var1+399, y = value, col=Var2) + geom_line() + 
        xlab('Wavelength') + ylab('Reflectance') + color

plot(Np)
plot(Cabp)
plot(Cwp)
plot(Cmp)