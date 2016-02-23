#' @title Negative log-likelihood for potentially constrained von Bertalanffy growth model (typically used internally).
#'
#' @description \code{vb_bind_nll} returns the negative log-likelihood for the von Bertalanffy model. Equality constraints across sexes can be implemented for any combination of parameters using the \code{binding} argument. Assumes a normal distribution currently.
#' @param theta A parameter vector of the same length as the maximum of \code{binding}. Unconstrained parameters take the order: lnlinfF, lnlinfM, lnkF, lnkM, lnnt0F, lnnt0M, lnsigmaF, lnsigmaM.
#' @param binding A (4x2) parameter index matrix with rows named (in order): "lnlinf", "lnk", "lnnt0", "lnsigma" and the left column for the female parameter index and right column for mal parameter index. Used to impose arbitrary equality constraints across the sexes (see Examples).
#' @param data A data.frame with columns: "age", "length" and "weights". "weights" are set to 1 or 0 for known females or males, respectively; proportions otherwise.
#' @return Negative log-likelihood:
#' set.seed(1010)
#' sim.dat <- sim_vb_data(nfemale = 50, nmale = 50, mean_ageF = 4, mean_ageM = 4,
#'                       growth_parF = c(linf = 30, k = 0.5, t0 = -1, sigma = 1),
#'                       growth_parM = c(linf = 25, k = 0.5, t0 = -1, sigma = 1),
#'                       mat_parF = c(A50 = 5, MR = 2), mat_parM = c(A50 = 3, MR = 2))
#' 
#' ## Unconstrained model 
#' binding <- matrix(c(1:8), ncol = 2, byrow = TRUE)
#' rownames(binding) <- c("lnlinf", "lnk", "lnnt0", "lnsigma")
#' colnames(binding) <- c("female", "male")
#' ## starting values 
#' start.par <- c(rep(log(25), 2), rep(log(0.2), 2), rep(log(3), 2), rep(log(1), 2))
#' vb_bind_nll(theta = start.par, binding = binding, data = data.frame(age = rep(1, 2), length = rep(10, 2), weights = c(1, 0)))

vb_bind_nll <- function(theta, binding, data) {
  linfF <- exp(theta[binding["lnlinf", "female"]])
  linfM <- exp(theta[binding["lnlinf", "male"]])
  kF <- exp(theta[binding["lnk", "female"]])
  kM <- exp(theta[binding["lnk", "male"]])
  t0F <- -exp(theta[binding["lnnt0", "female"]])
  t0M <- -exp(theta[binding["lnnt0", "male"]])
  sigmaF <- exp(theta[binding["lnsigma", "female"]])
  sigmaM <- exp(theta[binding["lnsigma", "male"]])
  muF <- linfF * (1 - exp(-kF * (data$age - t0F)))
  llF <- sum(data$weights * dnorm(data$length, mean = muF, 
                                  sd = sigmaF, log = TRUE))
  muM <- linfM * (1 - exp(-kM * (data$age - t0M)))
  llM <- sum((1 - data$weights) * dnorm(data$length, mean = muM, 
                                        sd = sigmaM, log = TRUE))
  return(-(llF + llM))
}
