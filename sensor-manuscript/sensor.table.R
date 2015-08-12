# Sensor accuracy and precision table

# Load packages {{{
library(data.table)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(PEcAnRTM)
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
fft.dat[, sensor := factor(sensor, levels=sensor.list)]
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
#}}}

# Prepare xtable {{{
setkey(sensor.table, sensor)
sensor.table <- sensor.table[sensor.list]
#sensor.proper <- c("ASD Field Spec", "AVIRIS NG", "AVIRIS Classic",
                   #"Hyperion", "CHRIS-Proba", "Landsat 5", "Landsat 7",
                   #"Landsat 8", "MODIS", "VIIRS", "AVHRR")
sensor.table[, sensor := sensor.proper[sensor]]
setcolorder(sensor.table, c("sensor", rnames, cnames))
setnames(sensor.table, "sensor", "Sensor")
setnames(sensor.table, cnames, sprintf("$\\alpha(\\mathrm{%s})$", params))
setnames(sensor.table, rnames, sprintf("$\\pi(\\mathrm{%s})$", params))
# }}}

# Write xtable {{{
library(xtable)
cap <- "
Mean accuracy ($\\alpha$) and precision ($\\pi$) by sensor.
"
cap <- gsub("\\n", " ", cap)
out.tab <- xtable(sensor.table, caption=cap, digits=4, label="tab:sensor")
print(out.tab, file="manuscript/tables/tab-sensor.tex", 
      sanitize.text.function=function(x) x,
      include.rownames = FALSE)
# }}}
