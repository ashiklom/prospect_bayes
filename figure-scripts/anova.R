#' Complex analysis of variance for FFT data

if(grepl("figure-scripts", getwd())) setwd("..")

source("preprocess.fft.R")
library(ggplot2)
library(MASS)

form.big <- "%s ~ 1 + plant.type + succession + plant.type:succession +
        plant.type:succession:Label + Height + Height:Age + Site + Site:Plot"
form <- function(x) formula(sprintf(form.big, x))

mod.N1 <- lm(form("N.mu"), data=fft)
N.step <- stepAIC(mod.N1, direction="both")

library(leaps)
N.leaps <- regsubsets(form("N.mu"), data=fft, really.big=T)
