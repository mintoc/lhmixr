#' @title Fit finite mixture maturity model.
#'
#' @description \code{maturity_mix} fits sex-specific maturity models where some of the animals are of unknown sex. Optimization is via the Expectation-Maximisation algorithm. 
#' @param start.list A list with a list called par containing starting values for: "mixprop", "betaF" (female parameters) and "betaM" (male parameters) (see Examples).
#' @param data A data.frame with columns: "maturity" (binary), "x" (single explanatory covariate currently), and "obs.sex". "obs.sex" must have values "female", "male", "unclassified".
#' @param binding A (kx2) parameter index matrix (not implemented yet)
#' @param maxiter.em Integer for maximum number of EM iterations (1e3 default).
#' @param reltol Relative tolerance for EM observed data log likelihood convergence (1e-8 default).
#' @param plot.fit Logical, if TRUE fit plotted per iteration (binary maturity jittered for visualisation). Red and blue circles are used for known females and males, respectively. Unclassified animals are plotted as triangle with the colour indicating the expected probability of being female or male (FALSE default).
#' @param verbose Logical, if TRUE iteration and observed data log-likelihood printed.
#' @param estimate.mixprop Logical, if TRUE the mixing proportion is estimated, otherwise fixed at the starting value.
#' @return List containing the components:
#' \item{logLik.vec}{Observed data log-likelihood at each iteration.}
#' \item{logLik}{Observed data log-likelihood on the last EM iteration.}
#' \item{complete_data}{Data frame of the data (re-ordered) with component probabilities (tau).}
#' \item{coefficients}{Parameter estimates (on logit scale) by sex and associated standard errors (not currently implemented).}
#' \item{vcov}{Estimated variance covariance matrix of the parameters estimated on the real line. Can be used to obtain parameter standard errors on the natural scale (not currently implemented).}
#' \item{convergence}{Binary with a "0" denoting convergence of the EM algorithm.}
#' @examples
#' set.seed(1011)
#' sim.dat<-sim_mat_data(nfemale = 30, nmale = 30, shapes = c(3, 3),
#'                      mat_parF = c(x50 = 0.5, xr = 0.1),
#'                      mat_parM = c(x50 = 0.4, xr = 0.1),
#'                      xcutoff = c(female = 0.4, male = 0.3))
#' 
#' start.list <- list(par = list(
#'                        mixprop = 0.5,
#'                        betaF = matrix(c(-5, 10)),
#'                        betaM = matrix(c(-5, 10))
#'                    ))
#' 
#' mat.fit <- maturity_mix(start.list = start.list, data = sim.dat)

maturity_mix <- function(start.list, data, binding= NULL, maxiter.em = 1e3, reltol = 1e-8, plot.fit = FALSE, verbose = TRUE, estimate.mixprop = TRUE){
    ## check mixprop starting values
    if(!"mixprop" %in% names(start.list[["par"]])){
        stop("No starting value for mixing proportion provided, specify 'mixprop = value' in start.list list")
    }
    ollike <- rep(NA, maxiter.em)
    ## split the data
    data$jmaturity <- jitter(data$maturity, factor = 0.1)
    classified.data <- data[data$obs.sex %in% c("female", "male"), ]
    unclassified.data <- data[data$obs.sex == "unclassified", ]
    ## plot setup
    if(plot.fit){
        blue2red<-colorRampPalette(c("blue", "red"))
        col.vec<-blue2red(100)
        breaks<-seq(-0.1, 1.1, length=100)
    }
    ## EM ITERATIONS 
    for(i in 1:maxiter.em){
        if(i==1){
            par <- start.list[["par"]]
        }
        ##--------
        ## E-STEP
        ##--------
        ## MIXING PROPORTION
        mixprop <- par[["mixprop"]]
        ## PARAMETERS
        betaF <- matrix(par[["betaF"]])
        betaM <- matrix(par[["betaM"]])
        ## unclassified probs
        pF.unclass <- plogis(cbind(1, unclassified.data$x) %*% betaF)
        pM.unclass <- plogis(cbind(1, unclassified.data$x) %*% betaM)
        ## classified means
        pF.class <- plogis(cbind(1, classified.data$x) %*% betaF)
        pM.class <- plogis(cbind(1, classified.data$x) %*% betaM)
        ## EXPECTED VALUE OF CLASSIFICATION
        ## classified data (known)
        classified.data$tau <- ifelse(classified.data$obs.sex == "female", 1, ifelse(classified.data$obs.sex == "male", 0, NA))
        ## classification for unclassified data (missing)
        unclassified.data$tau <- get_maturity_post_prob(mixprop = mixprop, pF = pF.unclass, pM = pM.unclass, data = unclassified.data)
        ## make the complete data
        complete.data <- rbind(classified.data, unclassified.data)
        ##-----------------------------
        ## FILL IN OBSERVED LIKELIHOOD 
        ##-----------------------------
        ## NOTE: IF ONLY INCLUDING IMMATURE IN MIXPROP CALCULATION OMIT MIXPROP FROM CLASSIFIED
        ll.F.class <- sum(classified.data$obs.sex == "female") * log(mixprop) +
            sum(dbinom(classified.data$maturity, size = 1, prob = pF.class, log=TRUE)[classified.data$obs.sex == "female"])
        ##
        ll.M.class <- sum(classified.data$obs.sex == "male") * log(1 - mixprop) +
            sum(dbinom(classified.data$maturity, size = 1, prob = pM.class, log=TRUE)[classified.data$obs.sex == "male"])
        ## unclassified component - finite mixture density
        ll.miss <- sum(log(
            mixprop * dbinom(unclassified.data$maturity, size = 1, prob = pF.unclass) +
            (1-mixprop) * dbinom(unclassified.data$maturity, size = 1, prob = pM.unclass)))
        ##
        ollike[i] <- ll.F.class + ll.M.class + ll.miss
        ##------
        ## PLOT
        ##------
        if(plot.fit){
            ##par(ask  = FALSE) ## so example runs through
            tau.col <- col.vec[cut(complete.data$tau, breaks)]
            par(mfrow=c(1, 1), mar = c(2, 2, 1, 1), oma = c(2, 2, 1, 1))
            x.pred <- seq(min(complete.data$x), max(complete.data$x), length = 50)
            ##
            predF <- c(plogis(cbind(1, x.pred) %*% betaF))
            predM <- c(plogis(cbind(1, x.pred) %*% betaM))
            plot(complete.data$x, complete.data$jmaturity,
                 pch=ifelse(complete.data$obs.sex=="unclassified",17, 19),
                 col=paste(tau.col,40, sep=""),
                 ylim=c(0, 1),
                 xlab="", ylab="")
            ##
            lines(x.pred, predF, col = "red")
            lines(x.pred, predM, col = "blue")
        }
        ##--------
        ## M-STEP (UPDATE PARAMETERS)
        ##--------
        ## Mixing proportion
        if(estimate.mixprop){
            mixprop <- sum(complete.data$tau)/length(complete.data$tau)
        }
        par[["mixprop"]] <- mixprop
        ## MATURITY MODEL
        complete.data$weights <- complete.data$tau
        fitF <- glm(cbind(maturity, 1 - maturity) ~ x, weights = weights, family = binomial, data = complete.data)
        fitM <- glm(cbind(maturity, 1 - maturity) ~ x, weights = 1 - weights, family = binomial, data = complete.data)
        betaF <- matrix(coef(fitF))
        betaM <- matrix(coef(fitM))
        par[["betaF"]] <- betaF
        par[["betaM"]] <- betaM
        ##--------
        ## OUTPUT
        ##--------
        if(verbose){
            cat(paste("EM iteration:", i, "|", "Observed data log-likelihood: ", ollike[i], "\n"))
        }
        ## CONVERGENCE
        if(i>2){
            ##if(abs(ollike[i] - ollike[i-1]) <  abstol | i == maxiter.em){
            if(abs(ollike[i] - ollike[i-1]) <  abs(ollike[i-1] * reltol) | i == maxiter.em){
                ## RESULTS OUTPUT
                res <- list()
                res$logLik.vec <- ollike[1:i]
                res$logLik <- ollike[i]
                res$complete.data <- complete.data
                res$coefficients <- list(mixprop = par[["mixprop"]],
                                         betaF = par[["betaF"]],
                                         betaM = par[["betaM"]])
                ##res$vcov <- par.vcov
                res$convergence <- ifelse(i == maxiter.em, 1, 0)
                return(res)
            }
        }
        ## clean-up within iteration
        rm(list = ls()[!ls()%in%c("classified.data", "unclassified.data","maxiter.em","par", "ollike","reltol","plot.fit", "col.vec", "breaks", "verbose", "estimate.mixprop")])
    }
}
