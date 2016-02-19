#' @title von Bertalanffy growth function.
#'
#' @description \code{vb_lengths} returns the predicted length-at-age for given named set of parameters for the von Bertalanffy growth function:
#' \deqn{L=L_{\infty}(1-e^{-k(A-t_0)})}
#'
#' @param theta A numeric vector with named values "linf", "k", "t0".
#' @param age A numeric vector of ages.
#' @return Predicted length-at-age.
#' @examples
#' vb_lengths(theta=c("linf"=30,"k"=0.2,"t0"=-1), age=0:10)
vb_lengths<- function(theta, age) {
 pred_length<-theta["linf"]*(1-exp(-theta["k"]*(age-theta["t0"])))
 return(pred_length)
}
