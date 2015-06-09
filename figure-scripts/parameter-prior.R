if(grepl("figure-scripts", getwd())) setwd("..")

library(ggplot2)
library(gridExtra)

load("data/p4.defaults.RData")
source("figure-scripts/common.R")

dlnormt <- function(lx, mu, sigma, MIN){
    x <- log(lx)
    d <- numeric(length(lx))
    #mu <- exp(lmu + lsigma^2/2)
    #sigma <- sqrt((exp(sigma^2)-1) * exp(2*lmu + lsigma^2))
    d[lx > MIN] <- dnorm(x, mu, sigma) / (1 - pnorm(MIN, mu, sigma, 1, 0))
    return(d)
}

pmu <- c("N"=log(1.4), "Cab"=log(25), "Cw"=log(0.011), "Cm"=log(0.005))
psd <- c("N"=1.1, "Cab"=0.9, "Cw"=0.9, "Cm"=0.9)
xs <- function(a, b) seq(a, b, length.out=100)
N.x <- xs(1, 5)
Cab.x <- xs(0, 100)
Cw.x <- xs(0, 0.05)
Cm.x <- xs(0, 0.04)
xx <- data.frame(N.x, Cab.x, Cw.x, Cm.x)

th.prior <- theme_bw() + 
    theme(text = element_text(size=11),
          axis.text = element_text(size=rel(0.7)),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
th.noY <- theme(axis.title.y = element_blank())
lab.N <- xlab("N")
lab.Cab <- xlab(expression("Cab ("*mu*"g cm"^-2*")"))
lab.Cw <- xlab("Cw (cm)")
lab.Cm <- xlab(expression("Cm (g cm"^-2*")"))
pN <- ggplot(xx) + aes(x=N.x) + xlab("N") + ylab("Density") +
    stat_function(fun=dlnormt, arg=list(pmu["N"], psd["N"], 1)) + th.prior
pCab <- ggplot(xx) + aes(x=Cab.x) + lab.Cab +
    stat_function(fun=dlnormt, arg=list(pmu["Cab"], psd["Cab"], 0)) + th.prior + th.noY
pCw <- ggplot(xx) + aes(x=Cw.x) + lab.Cw + ylab("Density") +
    stat_function(fun=dlnormt, arg=list(pmu["Cw"], psd["Cw"], 0)) + th.prior
pCm <- ggplot(xx) + aes(x=Cm.x) + lab.Cm + th.noY +
    stat_function(fun=dlnormt, arg=list(pmu["Cm"], psd["Cm"], 0)) + th.prior + th.noY
png.plot("manuscript/figures/priors.png", h=4, w=4)
grid.arrange(pN, pCab, pCw, pCm, nrow=2)
dev.off()
