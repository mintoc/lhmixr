#' @title Fit finite mixture growth models.
#'
#' @description \code{growth_mix} fits sex-specific growth models where some of the animals are of unknown sex. Optimization is via the Expectation-Maximisation algorithm. Mean form for growth is flexible (tested on "nls", "gam" in package mgcv, "segmented" in package "segmented"). Assumes a normal distribution for growth currently.
#' @param start.fit A list with starting models for (optionally): "female_growth_fit", "male_growth_fit", "mixprop".
#' @param data A data.frame with columns: "age", "length" (for growth) and "obs.sex". "obs.sex" must have values "female", "male", "immature". 
#' @param maxiter.em Integer for maximum number of EM iterations (1e3 default).
#' @param abstol Tolerance for EM observed data log likelihood convergence (1e-9 default).
#' @param plot.fit Logical, if TRUE fit plotted per iteration.
#' @param verbose Logical, if TRUE iteration and observed data log-likelihood printed.
#' @return List containing the components:
#' \item{ollike}{Observed data log-likelihood at each iteration.}
#' \item{complete_data}{Data frame of the data (re-ordered) with component probabilities (tau).}
#' \item{mixprop}{Estimate of the female mixing proportion.}
#' \item{par}{A list containing final complete data model estimates.}
#' \item{convergence}{Binary with a "0" denoting convergence.}
#' @examples
#' set.seed(1010)
#' sim.dat <- sim_vb_data(nfemale = 50, nmale = 50, mean_ageF = 4, mean_ageM = 4,
#'                       growth_parF = c(linf = 30, k = 0.5, t0 = -1, sigma = 1),
#'                       growth_parM = c(linf = 25, k = 0.5, t0 = -1, sigma = 1),
#'                       mat_parF = c(A50 = 5, MR = 2), mat_parM = c(A50 = 3, MR = 2))
#' 
#' ## set weights to one initially
#' sim.dat$weights <- 1
#' 
#' ## Additive growth model fit
#' start.fit<-list(
#'             female_growth_fit = mgcv::gam(length ~ s(age, k=5),
#'                                           data=subset(sim.dat, obs.sex=="female")),
#'             male_growth_fit = mgcv::gam(length ~ s(age, k=5),
#'                                         data=subset(sim.dat, obs.sex=="male")),
#'             mixprop = 0.5
#'             )
#'
#' add.fit <- growth_mix(data = sim.dat, start.fit = start.fit)
#'

growth_mix <- function(start.fit, data, maxiter.em = 1e3, abstol = 1e-9, plot.fit = TRUE, verbose = TRUE){
  ## check mixprop starting values
  if(!"mixprop" %in% names(start.fit)){
    stop("No starting value for mixing proportion provided, specify 'mixprop = value' in start.fit list")
  }
  ## observed log likelihood container
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
  ## EM ITERATIONS 
  for(i in 1:maxiter.em){
    if(i==1){
      par <- start.fit
    }
    ##--------
    ## E-STEP
    ##--------
    ## MIXING PROPORTION
    mixprop <- par[["mixprop"]]
    female_growth_fit <- par[["female_growth_fit"]]
    male_growth_fit <- par[["male_growth_fit"]]
    ## unclassified means
    muF.unclass <- predict(female_growth_fit, newdata = unclassified.data)
    muM.unclass <- predict(male_growth_fit, newdata = unclassified.data)
    ## classified means
    muF.class <- predict(female_growth_fit, newdata = classified.data)
    muM.class <- predict(male_growth_fit, newdata = classified.data)
    ## maximum likelihood estimates of weighted sigma
    sigmaF <- sqrt(sum(female_growth_fit$weights * residuals(female_growth_fit, "response")^2) / sum(female_growth_fit$weights))
    sigmaM <- sqrt(sum(male_growth_fit$weights * residuals(male_growth_fit, "response")^2) / sum(male_growth_fit$weights))
    ## EXPECTED VALUE OF CLASSIFICATION
    ## classified data (known)
    classified.data$tau <- ifelse(classified.data$obs.sex == "female", 1, ifelse(classified.data$obs.sex == "male", 0, NA))
    ## classification for unclassified data (missing)
    unclassified.data$tau <- get_growth_post_prob(mixprop = mixprop, muF = muF.unclass, muM = muM.unclass, sigmaF = sigmaF, sigmaM = sigmaM, data = unclassified.data)
    ## make the complete data
    complete.data <- rbind(classified.data, unclassified.data)
    ##-----------------------------
    ## FILL IN OBSERVED LIKELIHOOD 
    ##-----------------------------
    ## NOTE: IF ONLY INCLUDING IMMATURE IN MIXPROP CALCULATION OMIT MIXPROP FROM CLASSIFIED
    ll.F.class <- sum(classified.data$obs.sex == "female") * log(mixprop) +
      sum(dnorm(classified.data$length, mean = muF.class, sd = sigmaF, log=TRUE)[classified.data$obs.sex == "female"])
    ##
    ll.M.class <- sum(classified.data$obs.sex == "male") * log(1 - mixprop) +
      sum(dnorm(classified.data$length, mean = muM.class, sd = sigmaM, log=TRUE)[classified.data$obs.sex == "male"])
    ## unclassified component - finite mixture density
    ll.miss <- sum(log(
                     mixprop * dnorm(unclassified.data$length, mean = muF.unclass, sd = sigmaF) +
                     (1-mixprop) * dnorm(unclassified.data$length, mean = muM.unclass, sd = sigmaM)))
    ##
    ollike[i] <- ll.F.class + ll.M.class + ll.miss
    ##------
    ## PLOT
    ##------
    if(plot.fit){
      tau.col <- col.vec[cut(complete.data$tau, breaks)]
      par(mfrow=c(1, 1), mar = c(2, 2, 1, 1), oma = c(2, 2, 1, 1))
      age.pred <- seq(min(complete.data$jitter.age), max(complete.data$jitter.age), length=50)
      plot(complete.data$jitter.age, complete.data$length,
           pch=ifelse(complete.data$obs.sex=="immature",17, 19),
           col=paste(tau.col,40, sep=""), ylim=c(0, max(complete.data$length)),
           xlab="", ylab="")
      ##
      lines(age.pred, predict(female_growth_fit, newdata=data.frame(age=age.pred)), col = "red")
      lines(age.pred, predict(male_growth_fit, newdata=data.frame(age=age.pred)), col = "blue")
      mtext(side = 2, line = 2.5, text = "Length")
    }
    ##--------
    ## M-STEP (UPDATE PARAMETERS) 
    ##--------
    ## Mixing proportion
    mixprop <- sum(complete.data$tau)/length(complete.data$tau)
    par[["mixprop"]] <- mixprop
    ## for segmented update also needs to happen on the inner linear model
    if(class(par[["female_growth_fit"]])[1] == "segmented"){
      female.lm <- lm(length ~ age, data = complete.data, weights = tau)
      male.lm <- lm(length ~ age, data = complete.data, weights = 1 - tau)
    }
    ## NOTE: SHOULD CHECK FOR CONVERGENCE HERE BUT VARIOUS METHODS STORE IT DIFFERENTLY
    ## GROWTH MODELS
    par[["female_growth_fit"]] <- update(par[["female_growth_fit"]], data=complete.data, weights = tau)
    par[["male_growth_fit"]] <- update(par[["male_growth_fit"]], data=complete.data, weights = 1 - tau)
    ##--------
    ## OUTPUT
    ##--------
    if(verbose){
      cat(paste("EM iteration:", i, "|", "Observed data log-likelihood: ", ollike[i], "\n"))
    }
    ## CONVERGENCE
    if(i>2){
      if(abs(ollike[i] - ollike[i-1]) <  abstol | i == maxiter.em){
        res <- list()
        res$ollike <- ollike[1:i]
        res$complete.data <- complete.data
        res$mixprop <- mixprop
        res$par <- par
        res$convergence <- ifelse(i == maxiter.em, 1, 0)
        return(res)
      }
    }
    ## clean-up within iteration
    rm(list = ls()[!ls()%in%c("classified.data", "unclassified.data","maxiter.em","par", "ollike","abstol","plot.fit", "col.vec", "breaks", "verbose")])
  }
}
