#' @title Gradient of the negative log-likelihood for potentially constrained von Bertalanffy growth model (typically used internally).
#'
#' @description \code{vb_bind_gr} returns the parameter gradients of negative log-likelihood for the von Bertalanffy model. Equality constraints across sexes can be implemented for any combination of parameters using the \code{binding} argument.
#' @param theta A parameter vector of the same length as the maximum of \code{binding}. Unconstrained parameters take the order: lnlinfF, lnlinfM, lnkF, lnkM, lnnt0F, lnnt0M, lnsigmaF, lnsigmaM.
#' @param binding A (4x2) parameter index matrix with rows named (in order): "lnlinf", "lnk", "lnnt0", "lnsigma" and the left column for the female parameter index and right column for mal parameter index. Used to impose arbitrary equality constraints across the sexes (see Examples).
#' @param data A data.frame with columns: "age", "length" and "weights". "weights" are set to 1 or 0 for known females or males, respectively; proportions otherwise.
#' @param distribution Character with options: "normal" or "lognormal".
#' @return Vector of parameter gradients:
#' 
#' ## Unconstrained model 
#' binding <- matrix(c(1:8), ncol = 2, byrow = TRUE)
#' rownames(binding) <- c("lnlinf", "lnk", "lnnt0", "lnsigma")
#' colnames(binding) <- c("female", "male")
#' ## starting values 
#' start.par <- c(rep(log(25), 2), rep(log(0.2), 2), rep(log(3), 2), rep(log(1), 2))
#' vb_bind_gr(theta = start.par, binding = binding, data = data.frame(age = rep(1, 2), length = rep(10, 2), weights = c(1, 0)))

vb_bind_gr <- function(theta, binding, data, distribution){
  lnlinfF <- theta[binding["lnlinf", "female"]]
  lnlinfM <- theta[binding["lnlinf", "male"]]
  lnkF <- theta[binding["lnk", "female"]]
  lnkM <- theta[binding["lnk", "male"]]
  lnnt0F <- theta[binding["lnnt0", "female"]]
  lnnt0M <- theta[binding["lnnt0", "male"]]
  lnsigmaF <- theta[binding["lnsigma", "female"]]
  lnsigmaM <- theta[binding["lnsigma", "male"]]  
  ##
  age <- data$age
  length <- data$length
  weights <- data$weights
  if(distribution == "normal"){
    all.gradients <- c(
                       ## d/dlnlinfF
                       sum(-exp(lnlinfF-2*lnsigmaF)*weights*(1-exp(-exp(lnkF)*(age+exp(lnnt0F))))*(length-exp(lnlinfF)*(1-exp(-exp(lnkF)*(age+exp(lnnt0F)))))),
                       ## d/dlnlinfM
                       sum(-exp(lnlinfM-2*lnsigmaM)*(1-weights)*(1-exp(-exp(lnkM)*(age+exp(lnnt0M))))*(length-exp(lnlinfM)*(1-exp(-exp(lnkM)*(age+exp(lnnt0M)))))),  
                       ## d/dkF
                       sum(-exp(-2*lnsigmaF)*weights*(age+exp(lnnt0F))*exp(-exp(lnkF)*(age+exp(lnnt0F))+lnkF+lnlinfF)*(length-exp(lnlinfF)*(1-exp(-exp(lnkF)*(age+exp(lnnt0F)))))),
                       ## d/dkM
                       sum(-exp(-2*lnsigmaM)*(1-weights)*(age+exp(lnnt0M))*exp(-exp(lnkM)*(age+exp(lnnt0M))+lnkM+lnlinfM)*(length-exp(lnlinfM)*(1-exp(-exp(lnkM)*(age+exp(lnnt0M)))))),  
                       ## d/dlnnt0F
                       sum(-weights*exp(-exp(lnkF)*(age+exp(lnnt0F))-2*lnsigmaF+lnnt0F+lnkF+lnlinfF)*(length-exp(lnlinfF)*(1-exp(-exp(lnkF)*(age+exp(lnnt0F)))))),  
                       ## d/dlnnt0M
                       sum(-(1-weights)*exp(-exp(lnkM)*(age+exp(lnnt0M))-2*lnsigmaM+lnnt0M+lnkM+lnlinfM)*(length-exp(lnlinfM)*(1-exp(-exp(lnkM)*(age+exp(lnnt0M)))))),
                       ## d/dlnsigmaF
                       sum(weights*(1-exp(-2*lnsigmaF)*(length-exp(lnlinfF)*(1-exp(-exp(lnkF)*(age+exp(lnnt0F)))))^2)),  
                       ## d/dlnsigmaM
                       sum((1-weights)*(1-exp(-2*lnsigmaM)*(length-exp(lnlinfM)*(1-exp(-exp(lnkM)*(age+exp(lnnt0M)))))^2))
                       )
  }
  if(distribution == "lognormal"){    
    all.gradients <- c(                     
                       ## d/dlnlinfF
                       sum(-exp(-2*lnsigmaF)*weights*(log(length)-log(exp(lnlinfF)*(1-exp(-exp(lnkF)*(age+exp(lnnt0F)))))+exp(2*lnsigmaF)/2)),
                       ## d/dlnlinfM
                       sum(-exp(-2*lnsigmaM)*(1-weights)*(log(length)-log(exp(lnlinfM)*(1-exp(-exp(lnkM)*(age+exp(lnnt0M)))))+exp(2*lnsigmaM)/2)),
                       ## d/dkF
                       sum(-exp(-2*lnsigmaF)*weights*(age+exp(lnnt0F))*exp(lnkF-exp(lnkF)*(age+exp(lnnt0F)))*(log(length)-log(exp(lnlinfF)*(1-exp(-exp(lnkF)*(age+exp(lnnt0F)))))+exp(2*lnsigmaF)/2)/(1-exp(-exp(lnkF)*(age+exp(lnnt0F))))),
                       ## d/dkM
                       sum(-exp(-2*lnsigmaM)*(1-weights)*(age+exp(lnnt0M))*exp(lnkM-exp(lnkM)*(age+exp(lnnt0M)))*(log(length)-log(exp(lnlinfM)*(1-exp(-exp(lnkM)*(age+exp(lnnt0M)))))+exp(2*lnsigmaM)/2)/(1-exp(-exp(lnkM)*(age+exp(lnnt0M))))),
                       ## d/dlnnt0F
                       sum(-exp(-2*lnsigmaF)*weights*exp(-exp(lnkF)*(age+exp(lnnt0F))+lnnt0F+lnkF)*(log(length)-log(exp(lnlinfF)*(1-exp(-exp(lnkF)*(age+exp(lnnt0F)))))+exp(2*lnsigmaF)/2)/(1-exp(-exp(lnkF)*(age+exp(lnnt0F))))),
                       ## d/dlnnt0M
                       sum(-exp(-2*lnsigmaM)*(1-weights)*exp(-exp(lnkM)*(age+exp(lnnt0M))+lnnt0M+lnkM)*(log(length)-log(exp(lnlinfM)*(1-exp(-exp(lnkM)*(age+exp(lnnt0M)))))+exp(2*lnsigmaM)/2)/(1-exp(-exp(lnkM)*(age+exp(lnnt0M))))),
                       ## d/dlnsigmaF
                       sum(-weights*(exp(-2*lnsigmaF)*(log(length)-log(exp(lnlinfF)*(1-exp(-exp(lnkF)*(age+exp(lnnt0F)))))+exp(2*lnsigmaF)/2)^2-log(length)+log(exp(lnlinfF)*(1-exp(-exp(lnkF)*(age+exp(lnnt0F)))))-exp(2*lnsigmaF)/2-1)),
                       ## d/dlnsigmaM
                       sum(-(1-weights)*(exp(-2*lnsigmaM)*(log(length)-log(exp(lnlinfM)*(1-exp(-exp(lnkM)*(age+exp(lnnt0M)))))+exp(2*lnsigmaM)/2)^2-log(length)+log(exp(lnlinfM)*(1-exp(-exp(lnkM)*(age+exp(lnnt0M)))))-exp(2*lnsigmaM)/2-1))
                       )
  }
  ## gradients should be summed across bound parameters
  ##return(all.gradients[sort(unique(c(binding)))])
  return(tapply(all.gradients, c(t(binding)), sum))
}
