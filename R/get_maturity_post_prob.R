#' @title Posterior probability of sex being female based on maturity alone
#'
#' @description \code{get_growth_post_prob} returns the probability of
#' the observation(s) arising from the female component given a set of parameters
#' and an assumed distribution (binomial). The component 
#' probability is given by Bayes' theorem. Used internally.
#' @param mixprop Numeric scalar of mixing proportion (overall sex ratio)
#' @param pF Numeric vector with predicted female maturity probability.
#' @param pM Numeric vector with predicted male maturity probability.
#' @param data A data.frame with binary column "maturity". 
#' @return Numeric vector of the posterior probability of being female.
#' @examples
#' get_maturity_post_prob(mixprop = 0.5, pF = 0.5, pM = 0.5, sigmaF = 1,
#'                           data = data.frame(maturity = 1))

get_maturity_post_prob <- function(mixprop, pF, pM, data){
    tau <- mixprop * dbinom(data$maturity, size = 1, prob = pF) /
        (mixprop * dbinom(data$maturity, size = 1, prob = pF) + (1-mixprop) * dbinom(data$maturity, size = 1, prob = pM))
    return(tau)
}
