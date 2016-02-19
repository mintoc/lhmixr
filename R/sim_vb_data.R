#' @title Simulate sex-specific von Bertalanffy data with missing classifications. 
#'
#' @description \code{sim_vb_data} simulates sex-specific growth data according to
#' the von Bertalanffy growth model and a logistic model governing maturity.
#' @param nfemale Numeric scalar for number of female observations.
#' @param nmale Numeric scalar for number of male observations.
#' @param mean_ageF Numeric scalar for female mean age - used to generate ages from rpois(,lambda=mean_ageF)
#' @param mean_ageM Numeric scalar for male mean age - used to generate ages from rpois(,lambda=mean_ageM)
#' @param growth_parF Named ("linf", "k", "t0", "sigma") numeric vector with female growth parameters
#' @param growth_parM Named ("linf", "k", "t0", "sigma") numeric vector with male growth parameters
#' @param mat_parF Named ("A50", "MR") numeric vector with female maturation parameters
#' A50 is the age at 50\% maturity, MR is age range between 25\% and 75\% maturity.  
#' @param mat_parM Named ("A50", "MR") numeric vector with male maturation parameters.
#' @return data.frame with columns "age", "length", "true.sex", "obs.sex" (observed sex assuming immature animals unclassified), "maturiy" (binary: 1 if mature; 0 if immature).
#' @examples
#' sim.dat<-sim_vb_data(nfemale = 30, nmale = 30, mean_ageF = 3, mean_ageM = 3,
#'                      growth_parF = c(linf = 30, k = 0.2, t0 = -1, sigma = 1),
#'                      growth_parM = c(linf = 25, k = 0.2, t0 = -1, sigma = 1),
#'                      mat_parF = c(A50 = 3, MR = 1), mat_parM = c(A50 = 2, MR = 1))
#' 
#' plot(jitter(sim.dat$age), sim.dat$length,
#'      xlim=c(0, max(sim.dat$age)), ylim=c(0, max(sim.dat$length)),
#'      col = c("red", "blue", "grey")[match(sim.dat$obs.sex,c("female", "male", "immature"))],
#'      pch = 19, xlab="age", ylab="Length")

sim_vb_data<-function(nfemale, nmale, mean_ageF, mean_ageM, growth_parF, growth_parM, mat_parF, mat_parM){
  ## female ages
  age.vecF<-rpois(nfemale, lambda = mean_ageF)
  ## male ages
  age.vecM<-rpois(nmale, lambda = mean_ageM)
  age.vec<-c(age.vecF, age.vecM)
  ## female lengths
  length.vecF<-rnorm(nfemale, mean = vb_lengths(theta = growth_parF, age = age.vecF), sd = growth_parF["sigma"])
  ## male lengths
  length.vecM<-rnorm(nmale, mean = vb_lengths(theta = growth_parM, age = age.vecM), sd = growth_parM["sigma"])
  length.vec<-c(length.vecF, length.vecM)
  ## female maturity 
  matF<-rbinom(nfemale, size = 1, prob = plogis(2*log(3)/mat_parF["MR"]*(age.vecF-mat_parF["A50"])))
  ## male maturity
  matM<-rbinom(nmale, size = 1, prob = plogis(2*log(3)/mat_parM["MR"]*(age.vecM-mat_parM["A50"])))
  mat<-c(matF, matM)
  ## true and observed sexes
  true.sex<-rep(c("female","male"), times = c(nfemale, nmale))
  obs.sex<-ifelse(mat == 1, true.sex, "immature")
  ##
  df<-data.frame(age = age.vec, length = length.vec, true.sex = true.sex, obs.sex = obs.sex, maturity = mat)
  return(df)
}
