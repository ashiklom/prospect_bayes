#' Arguments: N, Cab, Car, Cw, Cm, sensor, ngibbs, run.name
library(PEcAnRTM)
if(exists("TEST")) {
    arg <- c(77, "identity", 10000, "testrun")
    do.mle <- TRUE
    quiet <- FALSE
} else {
    arg <- commandArgs(trailingOnly=TRUE)
    do.mle <- TRUE
    quiet <- TRUE
}

#' Process arguments and load data
load("simulation.inputs.RData")
parnames <- params.prospect5
index <- as.numeric(arg[1])
true.params <- par.mat[index,]
names(true.params) <- parnames
obs.raw <- obs.mat[,index]
noise <- noise.mat[,index]
sensor <- arg[2]
ngibbs <- as.numeric(arg[3])
run.name <- arg[4]

#' Set up inversion parameters
prior.vals <- prior.defaultvals.prospect(sd.inflate = 3)
start.params <- rlnorm(5, prior.vals$mu, prior.vals$sigma)
names(start.params) <- parnames
start.params["N"] <- start.params["N"] + 1
pm <- c(1,0,0,0,0)
data(sensor.rsr)
model <- function(params, constants) spectral.response(prospect(params, 5)[,1], sensor)
prior <- priorfunc.prospect(prior.vals$mu, prior.vals$sigma)
obs <- spectral.response(obs.raw, sensor)

#' Perform inversion
samples <- invert.slow(observed = obs,
                       inits = start.params,
                       constants = NULL,
                       ngibbs = ngibbs,
                       prior = prior,
                       pm = pm,
                       model = model,
                       do.mle = do.mle,
                       quiet = quiet)

#' Process samples
samples.sub <- burnin.thin(samples)
samples.summary <- summary.mvnorm(samples.sub)

l <- list(noise = noise, obs.raw = obs.raw, sensor = sensor, fname = run.name, samples = samples)
assign(run.name, c(as.list(true.params), as.list(samples.summary), l))

if(exists("TEST")){
    par(mfrow = c(5,2))
    sb <- samples[-5000:0,1:5]
    #sb <- samples[-(floor(ngibbs*0.5)):0,1:5]
    for(i in 1:5){
        plot(samples[,i], type='l')
        abline(h = true.params[i], col="red")
        plot(density(sb[,i]))
        abline(v = true.params[i], col="red")
    }
    pplt <- function() pairs(sb, pch=".")
} else {
    save(list=run.name, file=sprintf("../results-simulation/%s.RData", run.name))
}

