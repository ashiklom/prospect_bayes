#'@title Prospect Bayesian Inversion
#'@author Alexey Shiklomanov

source("prospect.R")
source("timer.R")
source("truncnorm.R")

guess.inits <- c(N=1.4, 
              Cab=30,
              Cw=0.017,
              Cm=0.006
              )

source("mle_inversion.R")

## NOTE: obs.spec must be a matrix as follows:
## Column 1 : Wavelengths (400:2500)
## Columns 2-n : Reflectance observations
## Use specdatproc script to generate correct matrices from data.
pinvbayes <- function(obs.spec,
                      ngibbs=100,
                      JumpRSD=0.1,
                      inits='random',
                      notifications=TRUE,
                      ar.step=10,
                      ar.target=0.75,
                      fname = "run_results/loctestrun.dat"
                      )
{
        nwl <- nrow(obs.spec)
        nspec <- ncol(obs.spec)
        JumpSD <- JumpRSD * unlist(guess.inits)
        JumpSD.alpha <- JumpSD


        ### Priors ###
        N.prior <- function(N) dnorm(N - 1, 0, 1.5, log=TRUE) + log(2)
        Cab.prior <- function(Cab) dlnorm(Cab, log(30), 0.9, log=TRUE)
        Cw.prior <- function(Cw) dlnorm(Cw, log(0.017), 0.5, log=TRUE)
        Cm.prior <- function(Cm) dlnorm(Cm, log(0.006), 0.9, log=TRUE)

        # Error
        pwl.p <- c(0.001, 0.001)          # Inverse gamma


        ### Initial conditions ###

        if(inits == "guess"){
                ic <- guess.inits
        } else if(inits == "mle"){
                ic <- p.invert(obs.spec)
        } else {
                ic <- c(abs(rnorm(1, 0, 1.5)) + 1,
                        rlnorm(1, log(30), 0.9),
                        rlnorm(1, log(0.017), 0.5),
                        rlnorm(1, log(0.006), 0.9)
                        )
        }
        N.i <- ic[1]
        Cab.i <- ic[2]
        Cw.i <- ic[3]
        Cm.i <- ic[4]
        print(sprintf("Initial conditions: N %g, Cab %g, Cw %g, Cm %g",
                      N.i, Cab.i, Cw.i, Cm.i))
        sd.i <- 1

        ### Extract indices for random effects ###
        # Matching patterns
        leaf <- "(^.*_)[0-9]{5}.csv"
        plot <- "^[0-9]{4}[A-Za-z]+[0-9]{2}__([A-Za-z]+_[A-Za-z0-9]+_).*"
        species <- "^[0-9]{4}[A-Za-z]+[0-9]{2}__[A-Za-z]+_[A-Za-z0-9]+_([A-Za-z]+_).*"
        
        re.leaf <- unique(gsub(leaf, "\\1", colnames(obs.spec)))
        n.la <- length(re.leaf)
        re.species <- unique(gsub(species, "\\1", colnames(obs.spec)))
        n.sp <- length(re.species)
        la.list <- lapply(re.leaf, grep, colnames(obs.spec))
        names(la.list) <- re.leaf
        sp.list <- lapply(re.species, grep, names(la.list))
        leafsearch <- function(i) {
                j <- seq_along(sp.list)[sapply(sp.list, function(x) i %in% x)]
                k <- which(sp.list[[j]] == i)
                return(c(j,k))
        }
        
        ls.list <- lapply(re.species, function(x)
                lapply(re.leaf[grep(x, re.leaf)], grep, colnames(obs.spec)[grep(x, colnames(obs.spec))]))
        n.ls <- lapply(ls.list, length)
        

        ### Random effects initial conditions
        # Overall leaf
        la.N.i <- rep(0, n.la)
        la.Cab.i <- rep(0, n.la)
        la.Cw.i <- rep(0, n.la)
        la.Cm.i <- rep(0, n.la)
        la.N.sd <- 1
        la.Cab.sd <- 1
        la.Cw.sd  <- 1
        la.Cm.sd <- 1

        # Species
        sp.N.i <- rep(0, n.sp)
        sp.Cab.i <- rep(0, n.sp)
        sp.Cw.i <- rep(0, n.sp)
        sp.Cm.i <- rep(0, n.sp)

        sp.N.sd <- 1
        sp.Cab.sd <- 1
        sp.Cw.sd  <- 1
        sp.Cm.sd <- 1

        # Species-leaf
        ls.N.i <- lapply(n.ls, function(x) rep(0, x))
        ls.Cab.i <- lapply(n.ls, function(x) rep(0, x))
        ls.Cw.i <- lapply(n.ls, function(x) rep(0, x))
        ls.Cm.i <- lapply(n.ls, function(x) rep(0, x))
        ls.N.sd <- rep(1, n.sp)
        ls.Cab.sd <- rep(1, n.sp)
        ls.Cw.sd  <- rep(1, n.sp)
        ls.Cm.sd <- rep(1, n.sp)


        # Random effects Prior
        randeff.s <- c(0.001, 0.001)              # Inverse gamma
                
        
        ### Shortcut functions ###
        prospect <- function(N, Cab, Cw, Cm) prospect4(N, Cab, Cw, Cm, n.a, cab.a, w.a, m.a)

        spec.error <- function(mod.spec, obs.spec){
                if(length(dim(obs.spec))){
                        return(-apply(obs.spec, 2, "-", mod.spec))
                } else {
                        return(mod.spec - obs.spec)
                }
        }

        likelihood <- function(guess.error, sd.i) sum(dnorm(guess.error, 0, sd.i, log=TRUE))

        # Precalculate first model
        prev.spec <- prospect(N.i, Cab.i, Cw.i, Cm.i)
        prev.error <- spec.error(prev.spec, obs.spec)

        ### MCMC storage
                
        prosp.pars <- c("N", "Cab", "Cw", "Cm")

        la.N.vec <- sprintf("la_N_%s", re.leaf)
        la.Cab.vec <- sprintf("la_Cab_%s", re.leaf)
        la.Cw.vec <- sprintf("la_Cw_%s", re.leaf)
        la.Cm.vec <- sprintf("la_Cm_%s", re.leaf)
        la.sd.vec <- sprintf("la_sd_%s", prosp.pars)
        
        sp.N.vec <- sprintf("sp_N_%s", re.species)
        sp.Cab.vec <- sprintf("sp_Cab_%s", re.species)
        sp.Cw.vec <- sprintf("sp_Cw_%s", re.species)
        sp.Cm.vec <- sprintf("sp_Cm_%s", re.species)
        sp.sd.vec <- sprintf("sp_sd_%s", prosp.pars)

        ls.N.vec <- sprintf("ls_N_%s", re.leaf)
        ls.Cab.vec <- sprintf("ls_Cab_%s", re.leaf)
        ls.Cw.vec <- sprintf("ls_Cw_%s", re.leaf)
        ls.Cm.vec <- sprintf("ls_Cm_%s", re.leaf)
        ls.sd.vec <- c(sprintf("ls_sd_N_%s", re.species),
                       sprintf("ls_sd_Cab_%s", re.species),
                       sprintf("ls_sd_Cw_%s", re.species),
                       sprintf("ls_sd_Cm_%s", re.species))

        header <- c(prosp.pars, 
                    la.N.vec, la.Cab.vec, la.Cw.vec, la.Cm.vec, la.sd.vec,
                    sp.N.vec, sp.Cab.vec, sp.Cw.vec, sp.Cm.vec, sp.sd.vec,
                    ls.N.vec, ls.Cab.vec, ls.Cw.vec, ls.Cm.vec, ls.sd.vec,
                    "resid_sd")
        write(header,
              ncolumns=length(header),
              file=fname, 
              sep=",", 
              append=FALSE)

        ## MCMC loop
        tstart <- proc.time()
        ar <- 0
        ar.alpha <- 0
        arp <- 0
        arp.alpha <- 0
        for(g in 1:ngibbs){
                arate <- ar/(4*g)
                if((g == 5) | (g %% (ngibbs/20) == 0) & notifications) laptime(tstart, g, ngibbs)

                if(g %% ar.step == 0){
                        #### TODO ####
                        ## Tweak JumpRSD based on acceptance rate
                        arate <- (ar - arp)/(4*ar.step)
                        arate.alpha <- (ar.alpha - arp.alpha)/(4*(2*n.la + n.sp)*ar.step)
                        JumpSD <- JumpSD * max(arate/ar.target, 0.001)
                        JumpSD.alpha <- JumpSD.alpha * max(arate.alpha/ar.target, 0.001)
                        arp <- ar
                        arp.alpha <- ar.alpha
                }
                ### Sample core PROSPECT parameters ###

                # Sample N
                if(notifications) print("Sampling N")
                guess.N <- rtnorm(1, N.i, JumpSD["N"], Min=1)
                N.spec <- function(i){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.spec <- prospect(guess.N + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        return(guess.error)
                }
                guess.error.re <- lapply(1:n.la, N.spec)
                guess.error <- do.call(cbind, guess.error.re)
                stopifnot(dim(guess.error) == dim(prev.error)) ## Throw error in case of regular expression mismatch
                guess.posterior <- likelihood(guess.error, sd.i) + N.prior(guess.N)              
                prev.posterior <- likelihood(prev.error, sd.i) + N.prior(N.i)
                jnum <- dtnorm(guess.N, N.i, JumpSD["N"], Min=1)
                jden <- dtnorm(N.i, guess.N, JumpSD["N"], Min=1)
                a <- exp((guess.posterior - jnum ) - (prev.posterior - jden ))
                if(is.na(a)) a <- -1
                if(a > runif(1)){
                        N.i <- guess.N
                        prev.error <- guess.error
                        ar <- ar + 1
                }
                # Sample Cab
                if(notifications) print("Sampling Cab")
                guess.Cab <- rtnorm(1, Cab.i, JumpSD["Cab"])
                Cab.spec <- function(i){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               guess.Cab + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        return(guess.error)
                }
                guess.error.re <- lapply(1:n.la, Cab.spec)
                guess.error <- do.call(cbind, guess.error.re)
                guess.posterior <- likelihood(guess.error, sd.i) + Cab.prior(guess.Cab)
                prev.posterior <- likelihood(prev.error, sd.i) + Cab.prior(Cab.i)
                jnum <- dlnorm(guess.Cab, Cab.i, JumpSD["Cab"])
                jden <- dlnorm(Cab.i, guess.Cab, JumpSD["Cab"])
                a <- exp((guess.posterior - jnum) - (prev.posterior - jden))
                if(is.na(a)) a <- -1
                if(a > runif(1)){
                        Cab.i <- guess.Cab
                        prev.error <- guess.error
                        ar <- ar + 1
                }

                # Sample Cw
                if(notifications) print("Sampling Cw")
                guess.Cw <- rtnorm(1, Cw.i, JumpSD["Cw"])
                Cw.spec <- function(i){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               guess.Cw + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        return(guess.error)
                }
                guess.error.re <- lapply(1:n.la, Cw.spec)
                guess.error <- do.call(cbind, guess.error.re)
                guess.posterior <- likelihood(guess.error, sd.i) + Cw.prior(guess.Cw)
                prev.posterior <- likelihood(prev.error, sd.i) + Cw.prior(Cw.i)
                jnum <- dlnorm(guess.Cw, Cw.i, JumpSD["Cw"])
                jden <- dlnorm(Cw.i, guess.Cw, JumpSD["Cw"])
                a <- exp((guess.posterior - jnum) - (prev.posterior - jden))
                if(is.na(a)) a <- -1
                if(a > runif(1)){
                        Cw.i <- guess.Cw
                        prev.error <- guess.error
                        ar <- ar + 1
                }

                # Sample Cm
                if(notifications) print("Sampling Cm")
                guess.Cm <- rtnorm(1, Cm.i, JumpSD["Cm"])
                Cm.spec <- function(i){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               guess.Cm + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        return(guess.error)
                }
                guess.error.re <- lapply(1:n.la, Cm.spec)
                guess.error <- do.call(cbind, guess.error.re)
                guess.posterior <- likelihood(guess.error, sd.i) + Cm.prior(guess.Cm)
                prev.posterior <- likelihood(prev.error, sd.i) + Cm.prior(Cm.i)
                jnum <- dlnorm(guess.Cm, Cm.i, JumpSD["Cm"])
                jden <- dlnorm(Cm.i, guess.Cm, JumpSD["Cm"])
                a <- exp((guess.posterior - jnum) - (prev.posterior - jden))
                if(is.na(a)) a <- -1
                if(a > runif(1)){
                        Cm.i <- guess.Cm
                        prev.error <- guess.error
                        ar <- ar + 1
                }

                #### Sample random effects ####

                ### Overall leaf effects ###
                if(notifications) print("Sampling overall leaf effects")
                ## Sample la.N
                for (i in 1:n.la){
                        guess.la.N <- la.N.i
                        guess.la.N[i] <- rnorm(1, la.N.i[i], JumpSD.alpha["N"])
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.spec <- prospect(N.i + guess.la.N[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.la.N[i], 0, la.N.sd, log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(la.N.i[i], 0, la.N.sd, log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                la.N.i <- guess.la.N
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample la.Cab
                for (i in 1:n.la){
                        guess.la.Cab <- la.Cab.i
                        guess.la.Cab[i] <- rnorm(1, la.Cab.i[i], JumpSD.alpha["Cab"])
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + guess.la.Cab[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.la.Cab[i], 0, la.Cab.sd, log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(la.Cab.i[i], 0, la.Cab.sd, log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                la.Cab.i <- guess.la.Cab
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample la.Cw
                for (i in 1:n.la){
                        guess.la.Cw <- la.Cw.i
                        guess.la.Cw[i] <- rnorm(1, la.Cw.i[i], JumpSD.alpha["Cw"])
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + guess.la.Cw[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.la.Cw[i], 0, la.Cw.sd, log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(la.Cw.i[i], 0, la.Cw.sd, log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                la.Cw.i <- guess.la.Cw
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample la.Cm
                for (i in 1:n.la){
                        guess.la.Cm <- la.Cm.i
                        guess.la.Cm[i] <- rnorm(1, la.Cm.i[i], JumpSD.alpha["Cm"])
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + guess.la.Cm[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.la.Cm[i], 0, la.Cm.sd, log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(la.Cm.i[i], 0, la.Cm.sd, log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                la.Cm.i <- guess.la.Cm
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample la precision
                v1 <- randeff.s[1] + n.la/2

                v2N <- randeff.s[2] + 0.5 * sum(la.N.i^2)
                preN <- rgamma(1, v1, v2N)
                la.N.sd <- 1/sqrt(preN)

                v2Cab <- randeff.s[2] + 0.5 * sum(la.Cab.i^2)
                preCab <- rgamma(1, v1, v2Cab)
                la.Cab.sd <- 1/sqrt(preCab)

                v2Cw <- randeff.s[2] + 0.5 * sum(la.Cw.i^2)
                preCw <- rgamma(1, v1, v2Cw)
                la.Cw.sd <- 1/sqrt(preCw)

                v2Cm <- randeff.s[2] + 0.5 * sum(la.Cm.i^2)
                preCm <- rgamma(1, v1, v2Cm)
                la.Cm.sd <- 1/sqrt(preCm)

                ### Species effects ###
                if(notifications) print("Sampling species effects")
                ## Sample sp.N
                for (i in 1:n.la){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.sp.N <- sp.N.i
                        guess.sp.N[j] <- rnorm(1, sp.N.i[j], JumpSD.alpha["N"])
                        guess.spec <- prospect(N.i + la.N.i[i] + guess.sp.N[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.sp.N[j], 0, sp.N.sd, log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(sp.N.i[j], 0, sp.N.sd, log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                sp.N.i <- guess.sp.N
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample sp.Cab
                for (i in 1:n.la){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.sp.Cab <- sp.Cab.i
                        guess.sp.Cab[j] <- rnorm(1, sp.Cab.i[j], JumpSD.alpha["Cab"])
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + guess.sp.Cab[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.sp.Cab[j], 0, sp.Cab.sd, log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(sp.Cab.i[j], 0, sp.Cab.sd, log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                sp.Cab.i <- guess.sp.Cab
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample sp.Cw
                for (i in 1:n.la){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.sp.Cw <- sp.Cw.i
                        guess.sp.Cw[j] <- rnorm(1, sp.Cw.i[j], JumpSD.alpha["Cw"])
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + guess.sp.Cw[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.sp.Cw[j], 0, sp.Cw.sd, log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(sp.Cw.i[j], 0, sp.Cw.sd, log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                sp.Cw.i <- guess.sp.Cw
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample sp.Cm
                for (i in 1:n.la){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.sp.Cm <- sp.Cm.i
                        guess.sp.Cm[j] <- rnorm(1, sp.Cm.i[j], JumpSD.alpha["Cm"])
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + guess.sp.Cm[j] + ls.Cm.i[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.sp.Cm[j], 0, sp.Cm.sd, log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(sp.Cm.i[j], 0, sp.Cm.sd, log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                sp.Cm.i <- guess.sp.Cm
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample sp precision
                v1 <- randeff.s[1] + n.sp/2

                v2N <- randeff.s[2] + 0.5 * sum(sp.N.i^2)
                preN <- rgamma(1, v1, v2N)
                sp.N.sd <- 1/sqrt(preN)

                v2Cab <- randeff.s[2] + 0.5 * sum(sp.Cab.i^2)
                preCab <- rgamma(1, v1, v2Cab)
                sp.Cab.sd <- 1/sqrt(preCab)

                v2Cw <- randeff.s[2] + 0.5 * sum(sp.Cw.i^2)
                preCw <- rgamma(1, v1, v2Cw)
                sp.Cw.sd <- 1/sqrt(preCw)

                v2Cm <- randeff.s[2] + 0.5 * sum(sp.Cm.i^2)
                preCm <- rgamma(1, v1, v2Cm)
                sp.Cm.sd <- 1/sqrt(preCm)


                ### Species-leaf effects ###
                if(notifications) print("Sampling species-leaf effects")
                ## Sample ls.N
                for (i in 1:n.la){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.ls.N <- ls.N.i
                        guess.ls.N[[j]][k] <- rnorm(1, ls.N.i[[j]][k], JumpSD.alpha["N"])
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + guess.ls.N[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.ls.N[[j]][k], 0, ls.N.sd[j], log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(ls.N.i[[j]][k], 0, ls.N.sd[j], log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                ls.N.i <- guess.ls.N
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample ls.Cab
                for (i in 1:n.la){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.ls.Cab <- ls.Cab.i
                        guess.ls.Cab[[j]][k] <- rnorm(1, ls.Cab.i[[j]][k], JumpSD.alpha["Cw"])
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + guess.ls.Cab[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.ls.Cab[[j]][k], 0, ls.Cab.sd[j], log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(ls.Cab.i[[j]][k], 0, ls.Cab.sd[j], log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                ls.Cab.i <- guess.ls.Cab
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample ls.Cw
                for (i in 1:n.la){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.ls.Cw <- ls.Cw.i
                        guess.ls.Cw[[j]][k] <- rnorm(1, ls.Cw.i[[j]][k], JumpSD.alpha["Cw"])
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + guess.ls.Cw[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + ls.Cm.i[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.ls.Cw[[j]][k], 0, ls.Cw.sd[j], log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(ls.Cw.i[[j]][k], 0, ls.Cw.sd[j], log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                ls.Cw.i <- guess.ls.Cw
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample ls.Cm
                for (i in 1:n.la){
                        jk <- leafsearch(i)
                        j <- jk[1]
                        k <- jk[2]
                        guess.ls.Cm <- ls.Cm.i
                        guess.ls.Cm[[j]][k] <- rnorm(1, ls.Cm.i[[j]][k], JumpSD.alpha["Cw"])
                        guess.spec <- prospect(N.i + la.N.i[i] + sp.N.i[j] + ls.N.i[[j]][k],
                                               Cab.i + la.Cab.i[i] + sp.Cab.i[j] + ls.Cab.i[[j]][k],
                                               Cw.i + la.Cw.i[i] + sp.Cw.i[j] + ls.Cw.i[[j]][k],
                                               Cm.i + la.Cm.i[i] + sp.Cm.i[j] + guess.ls.Cm[[j]][k])
                        guess.error <- prev.error
                        guess.error[,la.list[[i]]] <- spec.error(guess.spec, obs.spec[, la.list[[i]]])
                        guess.posterior <- likelihood(guess.error, sd.i) + dnorm(guess.ls.Cm[[j]][k], 0, ls.Cm.sd[j], log=TRUE)
                        prev.posterior <- likelihood(prev.error, sd.i) + dnorm(ls.Cm.i[[j]][k], 0, ls.Cm.sd[j], log=TRUE)
                        a <- exp(guess.posterior - prev.posterior)
                        if(is.na(a)) a <- -1
                        if(a > runif(1)){
                                ls.Cm.i <- guess.ls.Cm
                                prev.error <- guess.error
                                ar.alpha <- ar.alpha + 1
                        }
                }

                ## Sample ls precision
                
                v1 <- sapply(1:n.sp, function(j) randeff.s[1] + n.ls[[j]]/2)
                v2N <- sapply(1:n.sp, function(j) randeff.s[2] + (n.ls[[j]] - 1) * var(ls.N.i[[j]]))
                preN <- rgamma(n.sp, v1, v2N)
                ls.N.sd <- 1/sqrt(preN)

                v2Cab <- sapply(1:n.sp, function(j) randeff.s[2] + (n.ls[[j]] - 1) * var(ls.Cab.i[[j]]))
                preCab <- rgamma(n.sp, v1, v2Cab)
                ls.Cab.sd <- 1/sqrt(preCab)

                v2Cw <- sapply(1:n.sp, function(j) randeff.s[2] + (n.ls[[j]] - 1) * var(ls.Cw.i[[j]]))
                preCw <- rgamma(n.sp, v1, v2Cw)
                ls.Cw.sd <- 1/sqrt(preCw)

                v2Cm <- sapply(1:n.sp, function(j) randeff.s[2] + (n.ls[[j]] - 1) * var(ls.Cm.i[[j]]))
                preCm <- rgamma(n.sp, v1, v2Cm)
                ls.Cm.sd <- 1/sqrt(preCm)


                ### Sample residual precision ### 
                nprec <- 1
                u1p <- nspec*nwl/2
                u2p <- (nspec*nwl - 1) * var(c(prev.error))
                u1 <- pwl.p[1] + u1p
                u2 <- pwl.p[2] + u2p
                pwl.i <- rgamma(nprec, u1, u2)
                sd.i <- 1/sqrt(pwl.i)

                # Store values 
                write(c(N.i, Cab.i, Cw.i, Cm.i,
                        la.N.i, la.Cab.i, la.Cw.i, la.Cm.i,
                        la.N.sd, la.Cab.sd, la.Cw.sd, la.Cm.sd,
                        sp.N.i, sp.Cab.i, sp.Cw.i, sp.Cm.i,
                        sp.N.sd, sp.Cab.sd, sp.Cw.sd, sp.Cm.sd,
                        unlist(ls.N.i), unlist(ls.Cab.i), unlist(ls.Cw.i), unlist(ls.Cm.i), 
                        ls.N.sd, ls.Cab.sd, ls.Cw.sd, ls.Cm.sd,
                        sd.i), 
                      ncolumns=length(header),
                      sep=",",
                      file=fname,
                      append=TRUE)
        }
}

