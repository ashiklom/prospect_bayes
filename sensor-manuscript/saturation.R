# Parameter saturation analysis

# Load packages {{{
library(data.table)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
load("../data/simulation.samp.dat.RData")
# }}}

# Parse and format data {{{
run.regex <- "(.*)[.](q[1-7])[.](.*)"
fft.dat <- simulation.dat
rm(simulation.dat)
fft.dat[, par.level := gsub(run.regex, "\\1", fname)]
fft.dat[, q.level := gsub(run.regex, "\\2", fname)]
par.regex <- "(N|Cab|Car|Cw|Cm)([1-7])"
fft.dat[, par.vary := gsub(par.regex, "\\1", par.level)]
fft.dat[, par.lnum := gsub(par.regex, "\\2", par.level)]
setnames(fft.dat, "sensor.sensor", "sensor")
sensor.order <- c("identity", "aviris.ng", "aviris.classic", "hyperion",
                  "chris.proba", "landsat5", "landsat7", "landsat8", "viirs",
                  "modis", "avhrr")
fft.dat[, sensor := factor(sensor, levels=sensor.order)]
fft.dat <- fft.dat[!is.na(sensor)]
# }}}

# Get relative uncertainty {{{
fft.dat[, N.rsd := N.sigma / N.mu]
fft.dat[, Cab.rsd := Cab.sigma / Cab.mu]
fft.dat[, Car.rsd := Car.sigma / Car.mu]
fft.dat[, Cw.rsd := Cw.sigma / Cw.mu]
fft.dat[, Cm.rsd := Cm.sigma / Cm.mu]
fft.dat[, N.cv := (N.mu - true.param.N.true)/true.param.N.true]
fft.dat[, Cab.cv := (Cab.mu - true.param.Cab.true)/true.param.Cab.true]
fft.dat[, Car.cv := (Car.mu - true.param.Car.true)/true.param.Car.true]
fft.dat[, Cw.cv := (Cw.mu - true.param.Cw.true)/true.param.Cw.true]
fft.dat[, Cm.cv := (Cm.mu - true.param.Cm.true)/true.param.Cm.true]
fft.dat[is.infinite(Cm.cv), Cm.cv := NA]
# }}}

# Isolate and melt relative SD and error {{{
params <- c("N", "Cab", "Car", "Cw", "Cm")
tnames <- sprintf("true.param.%s.true", params)
rnames <- sprintf("%s.rsd", params)
cnames <- sprintf("%s.cv", params)
idvars <- c("sensor", "par.vary", "par.lnum", "q.level")
sensor.table <- fft.dat[, lapply(.SD, mean, na.rm=TRUE), by=sensor, .SDcols=c(rnames, cnames)]
library(xtable)
cap <- "
Mean parameter relative standard deviation (RSD = SD / Mean) and relative error (CV = Estimate - True/ True) by sensor.
"
cap <- gsub("\\n", " ", cap)
out.tab <- xtable(sensor.table, caption=cap, digits=4, label="tab:sensor")
print(out.tab, file="manuscript/tables/tab-sensor.tex")

#fft.rsd <- fft.dat[, c(tnames, rnames, idvars), with=FALSE]
#fft.rm <- data.table(melt(fft.rsd, id.vars=c(tnames, idvars)))
#fft.rm <- fft.rm[mapply(grepl, par.vary, variable)]
#fft.err <- fft.dat[, c(tnames, cnames, idvars), with=FALSE]
#fft.em <- data.table(melt(fft.err, id.vars=c(tnames, idvars)))
#fft.em <- fft.em[mapply(grepl, par.vary, variable)]
# }}}

# Plot all {{{
p.theme <- theme_bw() + theme(axis.text.x = element_text(angle=90))
p.unc <- ggplot(fft.rm) + aes(x=sensor, y=value, fill=variable) + 
    facet_grid(q.level ~ par.lnum) + geom_bar(stat="identity", position="dodge") +
    p.theme + ylab("Parameter relative SD")
png("manuscript/figures/uncertainty.all.png", width=10, height=7, units="in", res=300)
plot(p.unc)
dev.off()

png("manuscript/figures/error.all.png", width=10, height=7, units="in", res=300)
p.err <- p.unc %+% fft.em + ylab("Parameter relative error")
plot(p.err)
dev.off()
#}}}

# Plot by parameter {{{
for (p in params) {
    p.unc <- ggplot(fft.rm[par.vary == p]) + aes(x=sensor, y=value) +
        facet_grid(q.level ~ par.lnum) + geom_bar(stat="identity") +
        p.theme + ylab(paste0(p, " relative SD"))
    fname <- sprintf(c("manuscript/figures/uncertainty.%s.png", "manuscript/figures/error.%s.png"), p)
    png(fname[1], width=10, height=7, units="in", res=300)
    plot(p.unc)
    dev.off()
    p.err <- p.unc %+% fft.em[par.vary == p] + ylab(paste0(p, " relative error"))
    png(fname[2], width=10, height=7, units="in", res=300)
    plot(p.err)
    dev.off()
}
# }}}

# Effect of sensor on error and uncertainty {{{
sensor.theme <- theme_bw() + 
    theme(text = element_text(size=8),
          axis.title = element_text(size=rel(0.85)),
          axis.text = element_text(size=rel(0.7)),
          axis.text.x = element_text(angle=90))

fft.sensor.error <- fft.dat[, lapply(.SD, mean, na.rm=TRUE), by=sensor, .SDcols = cnames]
fft.sensor.error.m <- melt(fft.sensor.error)
sensor.error.plot <- ggplot(fft.sensor.error.m) + aes(x=sensor, y=value) + geom_bar(stat="identity") +
    facet_grid(variable~., scales="free") + ylab("Mean relative parameter error (CV)") + 
    sensor.theme

fft.sensor.uncertainty <- fft.dat[, lapply(.SD, mean, na.rm=TRUE), by=sensor, .SDcols = rnames]
fft.sensor.uncertainty.m <- melt(fft.sensor.uncertainty)
sensor.uncertainty.plot <- ggplot(fft.sensor.uncertainty.m) + aes(x=sensor, y=value) + geom_bar(stat="identity") +
    facet_grid(variable~., scales="free") + ylab("Mean relative parameter uncertainty (CV)") + 
    sensor.theme

png("manuscript/figures/sensor.png", width=6, height=4, units="in", res=300)
grid.arrange(sensor.error.plot, sensor.uncertainty.plot, nrow=1)
dev.off()
# }}}

# Parameter saturation {{{
fft.dat[, q.level.num := as.numeric(gsub("q","",q.level))]
sat <- ggplot(fft.dat) + aes(x=true.param.N.true, y=N.cv) +
    geom_bar(stat="identity") + xlab("True value") + ylab("Estimate") + facet_grid(sensor ~ q.level)
plot(sat)
# }}}
