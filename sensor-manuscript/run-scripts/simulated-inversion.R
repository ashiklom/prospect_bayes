#' Arguments: N, Cab, Car, Cw, Cm, sensor, ngibbs, run.name
if(exists("TEST")) {
    arg <- c(1.4, 30, 8.5, 0.01, 0.008, 
             "aviris.classic", 500, "testrun")
} else {
    arg <- commandArgs(trailingOnly=TRUE)
}

#' Extract arguments
parnames <- c("N", "Cab", "Car", "Cw", "Cm")
true.params <- as.numeric(arg[1:5])
names(true.params) <- parnames
sensor <- arg[6]
ngibbs <- as.numeric(arg[7])
run.name <- arg[8]

#' Set up inversion parameters
require(PEcAnRTM)
prior.vals <- prior.defaultvals.prospect(sd.inflate = 3)
start.params <- rlnorm(5, prior.vals$mu, prior.vals$sigma)
names(start.params) <- parnames
start.params["N"] <- start.params["N"] + 1
pm <- c(1,0,0,0,0)
model <- function(params, constants) spectral.response(prospect(params, 5)[,1], sensor)
prior <- priorfunc.prospect(prior.vals$mu, prior.vals$sigma)

#' Simulate data
obs <- model(true.params, NULL)
samples <- invert.slow(observed = obs,
                       inits = start.params,
                       constants = NULL,
                       ngibbs = ngibbs,
                       prior = prior,
                       pm = pm,
                       model = model,
                       do.mle = TRUE,
                       quiet = TRUE)

assign(run.name, list(true.params = true.params, 
                      inits = start.params,
                      sensor = sensor, 
                      samples = samples))

if(exists("TEST")){
    par(mfrow = c(5,2))
    for(i in 1:5){
        plot(samples[,i], type='l')
        abline(h = true.params[i], col="red")
        plot(density(samples[-(floor(ngibbs*0.5)):0, i]))
        abline(v = true.params[i], col="red")
    }
} else {
    save(list=run.name, file=sprintf("../results-simulation/%s.RData", run.name))
}

