library(PEcAnRTM)
# PROSPECT 5
true.params <- c(1.7, 32.8, 8.5, 0.013, 0.008)
start.params <- c(3, 10, 20, 0.1, 0.1)
pmean <- c(0.7, 32.81, 8.51, 0.0129, 0.0077)
psd <- c(0.6, 17.87, 3.2, 0.0073, 0.0035)*1000
pm <- c(1,0,0,0,0)
model <- function(params, constants) spectral.response(prospect(params, 5)[,1], "aviris.classic")
obs <- model(true.params, NULL)
parnames <- c("N", "Cab", "Car", "Cw", "Cm")

# PROSPECT 4
#true.params <- c(1.7, 32.8, 0.013, 0.008)
#start.params <- c(3, 10, 0.1, 0.1)
#pmean <- c(1.7, 32.81, 0.0129, 0.0077)
#psd <- c(0.6, 17.87, 0.0073, 0.0035)*1000
#pm <- c(1,0,0,0)
#obs <- prospect(true.params, 4)[,1]
#model <- function(params, constants) prospect(params, 4)[,1]
#parnames <- c("N", "Cab", "Cw", "Cm")

names(true.params) <- parnames
names(start.params) <- parnames
npar <- length(true.params)
ln.mu <- function(m,s) log(m / sqrt(1 + (s/m)^2))
ln.sd <- function(m,s) sqrt(log(1 + (s/m)^2))
pmu <- ln.mu(pmean, psd)
psigma <- ln.sd(pmean, psd)
prior <- function(x) {
    x[1] <- x[1] - 1
    priors <- mapply(dlnorm, x, pmu, psigma, log=TRUE)
    return(sum(priors))
}

fit <- invert.lsq(obs, start.params, NULL, model, pm)

Rprof()
samples <- invert.slow(observed = obs,
                       inits = start.params,
                       constants = NULL,
                       ngibbs = 10000,
                       prior = prior,
                       pm = pm,
                       model = model,
                       do.mle = TRUE,
                       quiet = FALSE)
Rprof(NULL)

par(mfrow=c(npar,2))
for(i in 1:npar){
    plot(samples[,i], type='l')
    abline(h = true.params[i], col="red")
    plot(density(samples[-8000:0,i]))
    abline(v = true.params[i], col="red")
}
result.params <- colMeans(samples[-8000:0,1:npar])
result.sd <- apply(samples[-8000:0,1:npar], 2, sd)
#plot(prospect(result.params, 5)[,1], type='l')
#lines(prospect(true.params, 5)[,1], col="red")

print("Results")
print(result.params)
print(result.sd)
print(100*result.sd/result.params)
print("True values")
print(true.params)
print("CV (%)")
print(100*(result.params - true.params)/true.params)

summaryRprof()

