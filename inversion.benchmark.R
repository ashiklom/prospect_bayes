source.path <- "~/Documents/Unsynced/PEcAnRTM/R/prosail.R"
source(source.path)

# Setup parameters
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

ng <- 1000
Rprof()
inv <- invert.sail(obs.1, inv.s.low, cst.1, ng, 
                   inv.priors, inv.pmin, pro4sail)
Rprof(NULL)
print(summaryRprof())
