# Analyze noise
library(PEcAnRTM)
library(ggplot2)
load("../data/noise.dat.RData")
setnames(noise.dat, 6, "residual.mu")
setnames(noise.dat, 12, "residual.sigma")

# {{{ Calculate absolute error
noise.dat[, N.error.abs := N.mu - N.true]
noise.dat[, Cab.error.abs := Cab.mu - Cab.true]
noise.dat[, Car.error.abs := Car.mu - Car.true]
noise.dat[, Cw.error.abs := Cw.mu - Cw.true]
noise.dat[, Cm.error.abs := Cm.mu - Cm.true]
# }}}

# {{{ Calculate relative error
noise.dat[, N.error := N.error.abs / N.true]
noise.dat[, Cab.error := Cab.error.abs / Cab.true]
noise.dat[, Car.error := Car.error.abs / Car.true]
noise.dat[, Cw.error := Cw.error.abs / Cw.true]
noise.dat[, Cm.error := Cm.error.abs / Cm.true]
# }}}

# {{{ Calculate true value z score
noise.dat[, N.z := N.error.abs / N.sigma]
noise.dat[, Cab.z := Cab.error.abs / Cab.sigma]
noise.dat[, Car.z := Car.error.abs / Car.sigma]
noise.dat[, Cw.z := Cw.error.abs / Cw.sigma]
noise.dat[, Cm.z := Cm.error.abs / Cm.sigma]
# }}}

# {{{ Calculate relative uncertainty
noise.dat[, N.rsd := N.sigma / N.mu]
noise.dat[, Cab.rsd := Cab.sigma / Cab.mu]
noise.dat[, Car.rsd := Car.sigma / Car.mu]
noise.dat[, Cw.rsd := Cw.sigma / Cw.mu]
noise.dat[, Cm.rsd := Cm.sigma / Cm.mu]
# }}}

# Aggregated means
sdcols <- c(sprintf("%s.%s", params.prospect5, rep(c("error", "z"), each=5)))
sigma.dat <- noise.dat[, lapply(.SD, mean), by=c("sigma", "fw"), .SDcols=sdcols]

# True values outside confidence limits
out.func <- function(x) sum(as.numeric(abs(x) > qnorm(0.99)))/length(x)
out.dat <- noise.dat[, lapply(.SD, out.func), by=c("sigma", "fw"), .SDcols=sdcols]

# Plots
p.out <- ggplot(out.dat[sigma <= 0.001]) + aes(x=sigma, y=Cab.z) +
    geom_line() + facet_wrap("fw")
plot(p.out)

p1 <- ggplot(noise.dat) + aes(x=factor(sigma), y=N.z) +
    geom_violin() + facet_wrap("fw", scales="free") + 
    geom_hline(y=c(-1.96,1.96), color="red", linetype="dashed")
plot(p1)

p.z <- ggplot(sigma.dat) + aes(x=sigma, y=abs(N.z)) +
    geom_point() + 
    facet_wrap("fw")
plot(p.z)

