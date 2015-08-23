# Sensor accuracy and precision table

# Load packages {{{
library(PEcAnRTM)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
load("../data/simulation.samp.dat.RData")
# }}}

# Parse and format data {{{
simulation.dat[, sensor := factor(sensor, levels=sensor.list)]
simulation.dat <- simulation.dat[!is.na(sensor)]
# }}}

# Get relative uncertainty {{{
simulation.dat[, N.rsd := N.N.sigma / N.mu]
simulation.dat[, Cab.rsd := Cab.Cab.sigma / Cab.mu]
simulation.dat[, Car.rsd := Car.Car.sigma / Car.mu]
simulation.dat[, Cw.rsd := Cw.Cw.sigma / Cw.mu]
simulation.dat[, Cm.rsd := Cm.Cm.sigma / Cm.mu]
simulation.dat[, N.cv := (N.mu - N)/N]
simulation.dat[, Cab.cv := (Cab.mu - Cab)/Cab]
simulation.dat[, Car.cv := (Car.mu - Car)/Car]
simulation.dat[, Cw.cv := (Cw.mu - Cw)/Cw]
simulation.dat[, Cm.cv := (Cm.mu - Cm)/Cm]
#simulation.dat[is.infinite(Cm.cv), Cm.cv := NA]
# }}}

# Isolate and melt relative SD and error {{{
tnames <- params.prospect5
rnames <- sprintf("%s.rsd", params.prospect5)
cnames <- sprintf("%s.cv", params.prospect5)
idvars <- "sensor"
sensor.table <- simulation.dat[, lapply(.SD, mean, na.rm=TRUE), by=sensor, 
                        .SDcols=c(rnames, cnames)]
#}}}

# Prepare xtable {{{
setkey(sensor.table, sensor)
sensor.table <- sensor.table[sensor.list]
sensor.table[, sensor := sensor.proper[sensor]]
setcolorder(sensor.table, c("sensor", rnames, cnames))
setnames(sensor.table, "sensor", "Sensor")
setnames(sensor.table, cnames, sprintf("$\\alpha(\\mathrm{%s})$", params.prospect5))
setnames(sensor.table, rnames, sprintf("$\\pi(\\mathrm{%s})$", params.prospect5))
# }}}

# Write xtable {{{
library(xtable)
cap <- "
Mean inaccuracy ($\\alpha$) and uncertainty ($\\pi$) by sensor.
"
cap <- gsub("\\n", " ", cap)
out.tab <- xtable(sensor.table, caption=cap, digits=4, label="tab:sensor")
out.tab.pre <- print(out.tab, file="", include.rownames=FALSE,
                     sanitize.text.function = function(x) x)
out.tab.post <- out.tab.pre
out.tab.post <- gsub("centering", "centerline{", out.tab.post)
out.tab.post <- gsub("(end\\{tabular\\})", "\\1\n\\}", out.tab.post)
cat(out.tab.post, file="manuscript/tables/tab-sensor.tex")
# }}}
