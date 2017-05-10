#' @title Posterior probability of sex being female based on growth
#'
#' @description \code{get_growth_post_prob} returns the probability of
#' the observation(s) arising from the female component given a set of growth parameters
#' and an assumed distribution (normal or lognormal). The component probability is
#' given by Bayes' theorem. Used internally.
#' @param mixprop Numeric scalar of mixing proportion (overall sex ratio)
#' @param muF Numeric vector with predicted female lengths
#' @param muM Numeric vector with predicted male lengths
#' @param sigmaF Numeric scalar for female residual standard deviation
#' @param sigmaM Numeric scalar for male residual standard deviation
#' @param data A data.frame with column "length". Note predicted means "muF" and "muM" must come from corresponding ages.
#' @param distribution Character with options: "normal" or "lognormal".
#' @return Numeric vector of the posterior probability of being female.
#' @source Minto, C., Hinde, J. and Coelho, R. (2017). Including unsexed individuals in sex-specific growth models.
#' \emph{Canadian Journal of Fisheries and Aquatic Sciences}. DOI: 10.1139/cjfas-2016-0450.
#' @examples
#' get_growth_post_prob(mixprop = 0.5, muF = 4, muM = 6, sigmaF = 1,
#'                           sigmaM = 1, data = data.frame(length = 4.5), distribution = "normal")
#' @export

get_growth_post_prob <- function(mixprop, muF, muM, sigmaF, sigmaM, data, distribution){
  if(distribution == "normal"){
    tau <- mixprop * dnorm(data$length, mean = muF, sd = sigmaF)/
      (mixprop * dnorm(data$length, mean = muF, sd = sigmaF)+(1-mixprop) * dnorm(data$length, mean = muM, sd = sigmaM))
  }
  if(distribution == "lognormal"){
    tau <- mixprop * dlnorm(data$length, meanlog = log(muF) - sigmaF^2 /2, sdlog = sigmaF)/
        (mixprop * dlnorm(data$length, meanlog = log(muF)  - sigmaF^2 /2, sdlog = sigmaF)+(1-mixprop) * dlnorm(data$length, meanlog = log(muM) - sigmaM^2 /2, sdlog = sigmaM))
    if(any(is.nan(tau))){
        warning("Very small probabilities of either group membership at given values. \n See function 'get_growth_post_prob'")
        ## should use underflow methods for very small numbers here
    }
  }  
  return(tau)
}
