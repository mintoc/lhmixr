#' @title Fit finite mixture von Bertalanffy growth model.
#'
#' @description \code{vb_growth_mix} fits sex-specific growth models where some of the animals are of unknown sex. Optimization is via the Expectation-Maximisation algorithm. Equality constraints across sexes can be implemented for any combination of parameters using the \code{binding} argument. Assumes a normal distribution currently.
#' @param start.list A list with a list called par containing starting values for: "mixprop", "growth.par" (see Examples).
#' @param data A data.frame with columns: "age", "length" and "obs.sex". "obs.sex" must have values "female", "male", "immature".
#' @param binding A (4x2) parameter index matrix with rows named (in order): "lnlinf", "lnk", "lnnt0", "lnsigma" and the left column for the female parameter index and right column for mal parameter index. Used to impose arbitrary equality constraints across the sexes (see Examples).  
#' @param maxiter.em Integer for maximum number of EM iterations (1e3 default).
#' @param abstol Tolerance for EM observed data log likelihood convergence (1e-8 default).
#' @param plot.fit Logical, if TRUE fit plotted per iteration. Red and blue circles are used for known females and males, respectively. Immature / unsexed animals are plotted as triangle with the colour indicating the expected probability of being female or male (FALSE default).
#' @param verbose Logical, if TRUE iteration and observed data log-likelihood printed.
#' @param optim.method Character, complete data optimisation method to use in \code{optim}.
#' @param estimate.mixprop Logical, if TRUE the mixing proportion is estimated, otherwise fixed at the starting value.
#' @param distribution Character with options: "normal" or "lognormal".
#' @return List containing the components:
#' \item{logLik.vec}{Observed data log-likelihood at each iteration.}
#' \item{logLik}{Observed data log-likelihood on the last EM iteration.}
#' \item{complete_data}{Data frame of the data (re-ordered) with component probabilities (tau).}
#' \item{coefficients}{Parameter estimates (on the real line) and associated standard errors on the real line.}
#' \item{vcov}{Estimated variance covariance matrix of the parameters estimated on the real line. Can be used to obtain parameter standard errors on the natural scale.}
#' \item{convergence}{Binary with a "0" denoting convergence of the EM algorithm.}
#' @examples
#' set.seed(1010)
#' sim.dat <- sim_vb_data(nfemale = 50, nmale = 50, mean_ageF = 4, mean_ageM = 4,
#'                       growth_parF = c(linf = 30, k = 0.5, t0 = -1, sigma = 0.1),
#'                       growth_parM = c(linf = 25, k = 0.5, t0 = -1, sigma = 0.1),
#'                       mat_parF = c(A50 = 5, MR = 2), mat_parM = c(A50 = 3, MR = 2),
#'                       distribution = "lognormal")
#' 
#' ## Model fit with contrained Brody's growth coefficient
#' ## Set up the constraint
#' binding <- matrix(c(1:2, rep(3, 2), 4:7), ncol = 2, byrow = TRUE)
#' rownames(binding) <- c("lnlinf", "lnk", "lnnt0", "lnsigma")
#' colnames(binding) <- c("female", "male")
#' ## starting values 
#' start.par <- c(c(log(30), log(25)), rep(log(0.3), 1), rep(log(1), 2), rep(log(.1), 2))
#' start.list <- list(par = list(mixprop = 0.5, growth.par = start.par))
#' vb.bind.fit <- vb_growth_mix(data = sim.dat, start.list = start.list,
#'                              binding = binding, distribution = "lognormal",
#'                              abstol = 1e-6)
#' options(device.ask.default = TRUE)
#'

vb_growth_mix <- function(start.list, data, binding, maxiter.em = 1e3, abstol = 1e-8, plot.fit = FALSE, verbose = TRUE, optim.method = "BFGS", estimate.mixprop = TRUE, distribution){
  ## check mixprop starting values
  if(!"mixprop" %in% names(start.list[["par"]])){
    stop("No starting value for mixing proportion provided, specify 'mixprop = value' in start.list list")
  }
  ## check length of the starting parameters
  if(max(binding) != length(start.list[["par"]][["growth.par"]])){
    stop("Mis-match in the length of growth.par and that specified by binding.")
  }
  ## observed log-likelihood container
  ollike <- rep(NA, maxiter.em)
  ## if plotting set up some variables
  if(plot.fit){
    blue2red<-colorRampPalette(c("blue","red"))
    col.vec<-blue2red(100)
    breaks<-seq(-0.1,1.1,length=100)
    data$jitter.age <- jitter(data$age)
  }
  ## split the data 
  classified.data <- data[data$obs.sex %in% c("female", "male"), ]
  unclassified.data <- data[data$obs.sex == "immature", ]
  ## define growth functions
  female_growth_fit <- function(x){linfF * (1 - exp(-kF * (x - t0F)))}
  male_growth_fit <- function(x){linfM * (1 - exp(-kM * (x - t0M)))}
  ## EM ITERATIONS 
  for(i in 1:maxiter.em){
    if(i==1){
      par <- start.list[["par"]]
    }
    ##--------
    ## E-STEP
    ##--------
    growth.par <- par[["growth.par"]]
    linfF <- exp(growth.par[binding["lnlinf", "female"]])
    linfM <- exp(growth.par[binding["lnlinf", "male"]])
    kF <- exp(growth.par[binding["lnk", "female"]])
    kM <- exp(growth.par[binding["lnk", "male"]])
    t0F <- - exp(growth.par[binding["lnnt0", "female"]])
    t0M <- - exp(growth.par[binding["lnnt0", "male"]])
    sigmaF <- exp(growth.par[binding["lnsigma", "female"]])
    sigmaM <- exp(growth.par[binding["lnsigma", "male"]])    
    ## MIXING PROPORTION
    mixprop <- par[["mixprop"]]
    ## unclassified means
    muF.unclass <- female_growth_fit(unclassified.data$age)
    muM.unclass <- male_growth_fit(unclassified.data$age)
    ## classified means
    muF.class <- female_growth_fit(classified.data$age)
    muM.class <- male_growth_fit(classified.data$age)
    ## EXPECTED VALUE OF CLASSIFICATION
    ## classified data (known)
    classified.data$tau <- ifelse(classified.data$obs.sex == "female", 1, ifelse(classified.data$obs.sex == "male", 0, NA))
    ## classification for unclassified data (missing)
    unclassified.data$tau <- get_growth_post_prob(mixprop = mixprop, muF = muF.unclass, muM = muM.unclass, sigmaF = sigmaF, sigmaM = sigmaM, data = unclassified.data, distribution = distribution)
    ## make the complete data
    complete.data <- rbind(classified.data, unclassified.data)
    ##-----------------------------
    ## FILL IN OBSERVED LIKELIHOOD 
    ##-----------------------------
    ## NOTE: IF ONLY INCLUDING IMMATURE IN MIXPROP CALCULATION OMIT MIXPROP FROM CLASSIFIED
    if(distribution == "normal"){
      ll.F.class <- sum(classified.data$obs.sex == "female") * log(mixprop) +
        sum(dnorm(classified.data$length, mean = muF.class, sd = sigmaF, log=TRUE)[classified.data$obs.sex == "female"])
      ##
      ll.M.class <- sum(classified.data$obs.sex == "male") * log(1 - mixprop) +
        sum(dnorm(classified.data$length, mean = muM.class, sd = sigmaM, log=TRUE)[classified.data$obs.sex == "male"])
      ## unclassified component - finite mixture density
      ll.miss <- sum(log(
                       mixprop * dnorm(unclassified.data$length, mean = muF.unclass, sd = sigmaF) +
                       (1-mixprop) * dnorm(unclassified.data$length, mean = muM.unclass, sd = sigmaM)))
    }
    if(distribution == "lognormal"){
      ## female classified
      ll.F.class <- sum(classified.data$obs.sex == "female") * log(mixprop) +
        sum(dlnorm(classified.data$length, meanlog = log(muF.class) - sigmaF^2 / 2, sdlog = sigmaF, log=TRUE)[classified.data$obs.sex == "female"])
      ## male classified
      ll.M.class <- sum(classified.data$obs.sex == "male") * log(1 - mixprop) +
        sum(dlnorm(classified.data$length, meanlog = log(muM.class) - sigmaM^2 / 2, sdlog = sigmaM, log=TRUE)[classified.data$obs.sex == "male"])
      ## unclassified component - finite mixture density
      ll.miss <- sum(log(
                       mixprop * dlnorm(unclassified.data$length, meanlog = log(muF.unclass) - sigmaF^2 / 2, sdlog = sigmaF) +
                       (1-mixprop) * dlnorm(unclassified.data$length, meanlog = log(muM.unclass) - sigmaM^2 / 2, sdlog = sigmaM)))
    }
    ##
    ollike[i] <- ll.F.class + ll.M.class + ll.miss
    ##------
    ## PLOT
    ##------
    if(plot.fit){
      ##par(ask  = FALSE) ## so example runs through
      tau.col <- col.vec[cut(complete.data$tau, breaks)]
      par(mfrow=c(1, 1), mar = c(2, 2, 1, 1), oma = c(2, 2, 1, 1))
      age.pred <- seq(min(complete.data$jitter.age), max(complete.data$jitter.age), length=50)
      plot(complete.data$jitter.age, complete.data$length,
           pch=ifelse(complete.data$obs.sex=="immature",17, 19),
           col=paste(tau.col,40, sep=""),
           ylim=c(0, max(complete.data$length)),
           xlim=c(0, max(complete.data$jitter.age)),
           xlab="", ylab="")
      ##
      lines(age.pred, female_growth_fit(age.pred), col = "red")
      lines(age.pred, male_growth_fit(age.pred), col = "blue")
      mtext(side = 2, line = 2.5, text = "Length")
    }
    ##--------
    ## M-STEP (UPDATE PARAMETERS)
    ##--------
    ## Mixing proportion
    if(estimate.mixprop){
      mixprop <- sum(complete.data$tau)/length(complete.data$tau)
    }
    par[["mixprop"]] <- mixprop
    ## GROWTH MODEL - USE NLMINB AS BFGS FLIPPING OUT OCCASSIONALLY, SEEMS FASTER ALSO
    complete.data$weights <- complete.data$tau
    vb_fit <- optim(vb_bind_nll, par = growth.par, gr = vb_bind_gr, binding = binding, data = complete.data, method = optim.method, distribution = distribution)
    ##vb_fit <- nlminb(objective = vb_bind_nll, start = growth.par, gradient = vb_bind_gr, binding = binding, data = complete.data, distribution = distribution)
    par[["growth.par"]] <- vb_fit$par
    ##--------
    ## OUTPUT
    ##--------
    if(verbose){
      cat(paste("EM iteration:", i, "|", "Observed data log-likelihood: ", ollike[i], "\n"))
    }
    ## CONVERGENCE
    if(i>2){
      if(abs(ollike[i] - ollike[i-1]) <  abstol | i == maxiter.em){
        ## STANDARD ERRORS
        ## one fit of observed data log-likelihood
        oll <- function(theta, estimate.mixprop, distribution){
          linfF <- exp(theta[binding["lnlinf", "female"]])
          linfM <- exp(theta[binding["lnlinf", "male"]])
          kF <- exp(theta[binding["lnk", "female"]])
          kM <- exp(theta[binding["lnk", "male"]])
          t0F <- - exp(theta[binding["lnnt0", "female"]])
          t0M <- - exp(theta[binding["lnnt0", "male"]])
          sigmaF <- exp(theta[binding["lnsigma", "female"]])
          sigmaM <- exp(theta[binding["lnsigma", "male"]])
          if(estimate.mixprop){
            mixprop <- plogis(theta[max(binding) + 1])
          }else{
            mixprop <- mixprop
          }
          ## predicted means
          ## unclassified
          muF.unclass <- linfF * (1 - exp(-kF * (unclassified.data$age - t0F))) 
          muM.unclass <- linfM * (1 - exp(-kM * (unclassified.data$age - t0M))) 
          ## classified 
          muF.class <- linfF * (1 - exp(-kF * (classified.data$age - t0F))) 
          muM.class <- linfM * (1 - exp(-kM * (classified.data$age - t0M))) 
          if(distribution == "normal"){
            ## female classified
            ll.F.class <- sum(classified.data$obs.sex == "female") * log(mixprop) +
              sum(dnorm(classified.data$length, mean = muF.class, sd = sigmaF, log=TRUE)[classified.data$obs.sex == "female"])
            ## male classified
            ll.M.class <- sum(classified.data$obs.sex == "male") * log(1 - mixprop) +
              sum(dnorm(classified.data$length, mean = muM.class, sd = sigmaM, log=TRUE)[classified.data$obs.sex == "male"])
            ## unclassified component - finite mixture density
            ll.miss <- sum(log(
                             mixprop * dnorm(unclassified.data$length, mean = muF.unclass, sd = sigmaF) +
                             (1-mixprop) * dnorm(unclassified.data$length, mean = muM.unclass, sd = sigmaM)))
          }
          if(distribution == "lognormal"){
            ## female classified
            ll.F.class <- sum(classified.data$obs.sex == "female") * log(mixprop) +
              sum(dlnorm(classified.data$length, meanlog = log(muF.class) - sigmaF^2 / 2, sdlog = sigmaF, log=TRUE)[classified.data$obs.sex == "female"])
            ## male classified
            ll.M.class <- sum(classified.data$obs.sex == "male") * log(1 - mixprop) +
              sum(dlnorm(classified.data$length, meanlog = log(muM.class) - sigmaM^2 / 2, sdlog = sigmaM, log=TRUE)[classified.data$obs.sex == "male"])
            ## unclassified component - finite mixture density
            ll.miss <- sum(log(
                             mixprop * dlnorm(unclassified.data$length, meanlog = log(muF.unclass) - sigmaF^2 / 2, sdlog = sigmaF) +
                             (1-mixprop) * dlnorm(unclassified.data$length, meanlog = log(muM.unclass) - sigmaM^2 / 2, sdlog = sigmaM)))
          }
          ##
          oll <- ll.F.class + ll.M.class + ll.miss
          return(-oll)
        }
        ## INCLUDE GRADIENTS HERE ALSO
        if(estimate.mixprop){
          oll.fit <- optim(fn = oll, par = c(par[["growth.par"]], qlogis(mixprop)), hessian = TRUE, control = list(maxit = 1e4),  estimate.mixprop = TRUE, distribution = distribution, method = optim.method)
        }else{
          oll.fit <- optim(fn = oll, par = c(par[["growth.par"]]), hessian = TRUE, control = list(maxit = 1e4),  estimate.mixprop = FALSE, distribution = distribution, method = optim.method)
        }
        ##print(oll.fit$par)
        ## check to make sure final optim fit close to EM
        if(!(round(-oll.fit$value / ollike[i], 4) == 1)){
          warning(paste("EM solution and optim solution differ by ", -oll.fit$value - ollike[i], ", final parameter values may differ from final EM values.", sep = ""))
        }
        ## collate the final estimates
        theta <- oll.fit$par
        par.vcov <- solve(oll.fit$hessian)
        par.se <- sqrt(diag(par.vcov))
        female.pars <- c(oll.fit$par[binding[, "female"]], oll.fit$par[max(binding) + 1])
        female.se <- c(par.se[binding[, "female"]], par.se[max(binding) + 1])
        male.pars <- c(oll.fit$par[binding[, "male"]], - oll.fit$par[max(binding) + 1])
        male.se <- c(par.se[binding[, "male"]], par.se[max(binding) + 1])
        theta.df <- data.frame(Parameter = c(rownames(binding), "logitpi"),
                               Female = female.pars,
                               Female.Std.Error = female.se,
                               Male = male.pars,
                               Male.Std.Error = male.se)
        ## RESULTS OUTPUT
        res <- list()
        res$logLik.vec <- ollike[1:i]
        res$logLik <- ollike[i]
        res$complete.data <- complete.data
        res$coefficients <- theta.df
        res$vcov <- par.vcov
        res$convergence <- ifelse(i == maxiter.em, 1, 0)
        return(res)
      }
    }
    ## clean-up within iteration
    rm(list = ls()[!ls()%in%c("classified.data", "unclassified.data","maxiter.em","par", "ollike","abstol","plot.fit", "col.vec", "breaks", "verbose", "vb_bind_nll", "binding", "optim.method", "estimate.mixprop", "distribution", "female_growth_fit", "male_growth_fit")])
  }
}
