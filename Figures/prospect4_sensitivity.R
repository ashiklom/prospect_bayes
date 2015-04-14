## Sensitivity analysis
library(ggplot2)
library(reshape2)
library(PEcAnRTM)

data(prospect4)
vars <- c("N", "Cab", "Cw", "Cm")
wl <- 400:2500
ns <- 20
minmax <- list(N = c(1,3),
               Cab = c(10,70),
               Cw = c(1e-5, 1e-2),
               Cm = c(1e-6, 1e-2))

seq.log <- function(x) seq(x[1], x[2], length.out = ns)
seq.list <- lapply(minmax, seq.log)

Ns <- sapply(seq.list[["N"]], function(x) prospect(4,c(x, 30, 0.01, 0.01)))
Cabs <- sapply(seq.list[["Cab"]], function(x) prospect(4,c(1.4, x, 0.01, 0.01)))
Cws <- sapply(seq.list[["Cw"]], function(x) prospect(4,c(1.4, 30, x, 0.01)))
Cms <- sapply(seq.list[["Cm"]], function(x) prospect(4,c(1.4, 30, 0.01, x)))
sens.spec <- list(N = Ns, Cab = Cabs, Cw = Cws, Cm = Cms)

sens.plot <- ggplot() +
        aes(x = Var1 + 399, y = value, color = Var2, group = Var2) + 
        geom_line() +
        scale_color_continuous(low="red", high="blue", trans="sqrt") +
        xlab("Wavelength") + 
        ylab("Reflectance") +
        guides(color=F)

p <- lapply(vars, function(x) {
        m <- melt(sens.spec[[x]])
        m$Var2 <- rep(seq.list[[x]], each=2101)
        sens.plot %+% m
})

for(pp in p) plot(pp)
