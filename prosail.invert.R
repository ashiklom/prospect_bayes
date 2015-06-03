library(PEcAnRTM)

## Sail parameters
#param.order <- c("N", "Cab", "Car", "Cbrown", "Cw", "Cm",
#                 "LIDFa", "LIDFb", "LIDFtype", "LAI", "q",
#                 "tts", "tto", "psi", "psoil")
# Parameter values
# N         1 - 4
# Cab       10 - 120
# Car       5 - 30
# Cbrown    0
# Cw        0.005 - 0.04
# Cm        0.002 - 0.014 
# LAI       0 - 10
# LIDFa     10 - 80 (degrees)
# q         0.2
# psoil     0 - 1
# tts       30
# tto       0
# psi       0

## Priors
lai.mu <- log(1.5)
lai.sd <- log(4)
Cab.mu <- log(30)
Cab.sd <- log(10)
Car.mu <- log(20)
Car.sd <- log(10)
Cw.mu <- -6.377
Cw.sd <- 0.5
inv.priors <- list()
inv.priors[[1]] <- function(LAI) dnorm(log(LAI), lai.mu, lai.sd, 1)
inv.priors[[2]] <- function(Cab) dnorm(log(Cab), Cab.mu, Cab.sd, 1)
inv.priors[[3]] <- function(Car) dnorm(log(Car), Car.mu, Car.sd, 1)
inv.priors[[4]] <- function(Cw) dnorm(log(Cw), Cw.mu, Cw.sd, 1)
inv.pmin <- c("LAI" = 0, "Cab" = 0, "Car" = 0, "Cw" = 0)
inv.s.low <- c("LAI" = 1, "Cab" = 10, "Car" = 5, "Cw" = 0.005)
inv.s.mid <- c("LAI" = 5, "Cab" = 40, "Car" = 20, "Cw" = 0.017)
inv.s.high <- c("LAI" = 10, "Cab" = 120, "Car" = 30, "Cw" = 0.04)

## Basic inversion
par.1 <- c("LAI" = 2, "Cab" = 35, "Car" = 15, "Cw" = 0.009)
cst.1 <- sail.constants(par.1)
pri.1 <- inv.priors
pmn.1 <- inv.pmin
obs.1 <- pro4sail(par.1, cst.1)

## Invert with different initial conditions
ng <- 1000
inv.1.low <- invert.sail(obs.1, inv.s.low, cst.1, ng, 
                         inv.priors, inv.pmin, pro4sail)
inv.1.mid <- invert.sail(obs.1, inv.s.mid, cst.1, ng, 
                         inv.priors, inv.pmin, pro4sail)
inv.1.high <- invert.sail(obs.1, inv.s.high, cst.1, ng, 
                          inv.priors, inv.pmin, pro4sail)
save(inv.1.low, inv.1.mid, inv.1.high, file = "inv1.Rdata")

# Diagnostic plots
png("inv.1.png")
par(mfrow=c(4,2))
burnin <- -1000:0
for(i in 1:4){
    plot(inv.1.low[,i], type='l')
    lines(inv.1.mid[,i], col=2)
    lines(inv.1.high[,i], col=3)
    plot(density(c(inv.1.mid[burnin,i], inv.1.high[burnin,i])),
         main = names(par.1)[i])
    abline(v = par.1[i])
}
dev.off()


## Basic AVIRIS inversion
obs.2 <- ps.aviris(par.1, cst.1)
inv.2.low <- invert.sail(obs.2, inv.s.low, cst.1, ng, 
                         inv.priors, inv.pmin, ps.aviris)
inv.2.mid <- invert.sail(obs.2, inv.s.mid, cst.1, ng, 
                         inv.priors, inv.pmin, ps.aviris)
inv.2.high <- invert.sail(obs.2, inv.s.high, cst.1, ng, 
                          inv.priors, inv.pmin, ps.aviris)
save(inv.2.low, inv.2.mid, inv.2.high, file = "inv2.Rdata")

# Diagnostic plots
png("inv.2.png")
par(mfrow=c(4,2))
burnin <- -1000:0
for(i in 1:4){
    plot(inv.2.low[,i], type='l')
    lines(inv.2.mid[,i], col=2)
    lines(inv.2.high[,i], col=3)
    plot(density(c(inv.2.low[burnin,i],
                   inv.2.mid[burnin,i], 
                   inv.2.high[burnin,i])),
         main = names(par.1)[i])
    abline(v = par.1[i])
}
dev.off()


## Inversion over range of LAI values
pm.lai <- matrix(rep(par.1, 20), nrow=20, byrow=T)
colnames(pm.lai) <- c("LAI", "Cab", "Car", "Cw")
pm.lai[,"LAI"] <- seq(0.5, 10, length.out=20)

inv.lai <- list()
for(l in 1:nrow(pm.lai)){
    obs <- ps.aviris(pm.lai[l,], cst.1)
    ch <- list()
    ch[[1]] <- invert.sail(obs, inv.s.low, cst.1, ng, 
                             inv.priors, inv.pmin, ps.aviris)
    ch[[2]] <- invert.sail(obs, inv.s.mid, cst.1, ng, 
                             inv.priors, inv.pmin, ps.aviris)
    ch[[3]] <- invert.sail(obs, inv.s.high, cst.1, ng, 
                              inv.priors, inv.pmin, ps.aviris)
    inv.lai[[l]] <- ch
}
save(inv.lai, file="inv.lai.Rdata")


## Inversion over a range of chlorophyll concentrations
pm.cab <- matrix(rep(par.1, 20), nrow=20, byrow=T)
colnames(pm.lai) <- c("LAI", "Cab", "Car", "Cw")
pm.cab[,"Cab"] <- seq(0.5, 10, length.out=20)

inv.cab <- list()
for(l in 1:nrow(pm.cab)){
    obs <- ps.aviris(pm.cab[l,], cst.1)
    ch <- list()
    ch[[1]] <- invert.sail(obs, inv.s.low, cst.1, ng, 
                             inv.priors, inv.pmin, ps.aviris)
    ch[[2]] <- invert.sail(obs, inv.s.mid, cst.1, ng, 
                             inv.priors, inv.pmin, ps.aviris)
    ch[[3]] <- invert.sail(obs, inv.s.high, cst.1, ng, 
                              inv.priors, inv.pmin, ps.aviris)
    inv.cab[[l]] <- ch
}
save(inv.cab, file="inv.cab.Rdata")


## Inversion over a range of carotenoid concentrations
pm.car <- matrix(rep(par.1, 20), nrow=20, byrow=T)
colnames(pm.lai) <- c("LAI", "Cab", "Car", "Cw")
pm.car[,"Car"] <- seq(0.5, 10, length.out=20)

inv.car <- list()
for(l in 1:nrow(pm.car)){
    obs <- ps.aviris(pm.car[l,], cst.1)
    ch <- list()
    ch[[1]] <- invert.sail(obs, inv.s.low, cst.1, ng, 
                             inv.priors, inv.pmin, ps.aviris)
    ch[[2]] <- invert.sail(obs, inv.s.mid, cst.1, ng, 
                             inv.priors, inv.pmin, ps.aviris)
    ch[[3]] <- invert.sail(obs, inv.s.high, cst.1, ng, 
                              inv.priors, inv.pmin, ps.aviris)
    inv.car[[l]] <- ch
}
save(inv.car, file="inv.car.Rdata")

