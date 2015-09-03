#' ---
#' title: Sensor inaccuracy and uncertainty table
#' author: Alexey Shiklomanov
#' output_format: html_document
#' ---

#' # Introduction
#' This script is used to generate the summary table for the sensor simulation 
#' experiment. For each parameter and sensor, the table shows on average how 
#' close the mean inversion estimate was to the true value (inaccuracy) and the 
#' average uncertainty (expressed as relative standard deviation) of the 
#' inversion estimate.

#' # Setup
#' First, we load dependencies and data. `PEcAnRTM` is loaded for parameter 
#' naming and for its implicit use of the `data.table` package; `xtable` allows
#' exporting of R data objects into automatically formatted LaTeX or HTML 
#' tables.  For information on how `simulation.samp.dat.RData` is generated, 
#' see the `load.sim.R` script.

library(PEcAnRTM)
library(xtable)
load("../data/simulation.samp.dat.RData")

#' To facilitate ordering of the table, we convert the representation of the 
#' sensor in the data table (NOTE: `simulation.dat` is a `data.table` object, 
#' NOT a `data.frame`, hence the different syntax) from character to factor. In 
#' addition, we exclude any values that do not have a sensor designation.

simulation.dat[, sensor := factor(sensor, levels=sensor.list)]
simulation.dat <- simulation.dat[!is.na(sensor)]

#' # Calculate statistics
#' This block of code computes the relevant statistics for each sensor. `rsd` 
#' is the relative standard deviation, which we use as a measure of 
#' uncertainty. `cv` is the error in the mean value of the inversion estimate 
#' relative to the true value, which we use as the measure of inaccuracy. We 
#' use the `data.table` column definition syntax (via the `:=` operator) to 
#' perform this efficiently.

simulation.dat[, N.rsd := 100 * N.sigma / N.mu]
simulation.dat[, Cab.rsd := 100 * Cab.sigma / Cab.mu]
simulation.dat[, Car.rsd := 100 * Car.sigma / Car.mu]
simulation.dat[, Cw.rsd := 100 * Cw.sigma / Cw.mu]
simulation.dat[, Cm.rsd := 100 * Cm.sigma / Cm.mu]
simulation.dat[, N.cv := 100 * (N.mu - N)/N]
simulation.dat[, Cab.cv := 100 * (Cab.mu - Cab)/Cab]
simulation.dat[, Car.cv := 100 * (Car.mu - Car)/Car]
simulation.dat[, Cw.cv := 100 * (Cw.mu - Cw)/Cw]
simulation.dat[, Cm.cv := 100 * (Cm.mu - Cm)/Cm]

#' Next, we use `data.table` aggregation syntax to compute the mean values by 
#' sensor of each of the above statistics for each parameter.

rnames <- sprintf("%s.rsd", params.prospect5)
cnames <- sprintf("%s.cv", params.prospect5)
sensor.table <- simulation.dat[, lapply(.SD, mean, na.rm=TRUE), by=sensor, 
                        .SDcols=c(rnames, cnames)]

#' # Table formatting
#' Rather than using ugly R variable names, we add more informative row and 
#' column labels to the table in preparation for export. `sensor.list` and 
#' `sensor.proper` are cross-referenced character vectors defined by PEcAnRTM 
#' containing the names of sensors as used by R (`sensor.list`) and properly 
#' formatted for publication (`sensor.proper`). 

setkey(sensor.table, sensor)
sensor.table <- sensor.table[sensor.list]
sensor.table[, sensor := sensor.proper[sensor]]
setcolorder(sensor.table, c("sensor", rnames, cnames))
setnames(sensor.table, "sensor", "Sensor")
setnames(sensor.table, cnames, sprintf("$\\alpha(\\mathrm{%s})$", params.prospect5))
setnames(sensor.table, rnames, sprintf("$\\pi(\\mathrm{%s})$", params.prospect5))

#' We use xtable to generate a LaTeX version of the table, specifying the 
#' caption, number of significant figures (`digits`), and TeX reference 
#' (`label`). We call the print statement with no file argument to store the 
#' result as a string. The `sanitize.text.function` defines how xtable deals 
#' with escape characters such as `\`; we set this function to identity to use 
#' exactly the text we gave because our column names already have TeX 
#' formatting. Finally, we post-process the resulting raw string to add the 
#' `centerline` option, which centers the table across the entire page rather 
#' than strictly conforming to the left margin as per the default behavior. The 
#' table is written to a `.tex` file in the `manuscript/tables` directory.

cap <- "
Mean uncertainty ($\\pi$) and inaccuracy ($\\alpha$) by sensor, expressed as percentage (x 100).
"
cap <- gsub("\\n", " ", cap)
out.tab <- xtable(sensor.table, caption=cap, digits=3, label="tab:sensor")
out.tab.pre <- print(out.tab, file="", include.rownames=FALSE,
                     sanitize.text.function = function(x) x)
out.tab.post <- out.tab.pre
out.tab.post <- gsub("centering", "centerline{", out.tab.post)
out.tab.post <- gsub("(end\\{tabular\\})", "\\1\n\\}", out.tab.post)
cat(out.tab.post, file="manuscript/tables/tab-sensor.tex")
