#' @title Simulate sex-specific maturity data with missing classifications. 
#'
#' @description \code{sim_mat_data} simulates sex-specific maturity data according to
#' a logistic maturity model.
#' @param nfemale Numeric scalar for number of female observations.
#' @param nmale Numeric scalar for number of male observations.

#' @param shapes Numeric vector for beta-distribution shape parameters - used to generate x values
#' @param mat_parF Named ("x50", "xr") numeric vector with female maturity parameters ("x50" is the x value at which 50\% are mature; "xr" is the x range between 25\% and 75\% mature)
#' @param mat_parM Named ("x50", "xr") numeric vector with male maturity parameters
#' @param xcutoff Named ("female", "male") numeric vector with x cutoffs below which immature animals are unclassified
#' @return data.frame with columns "x", "true.sex", "obs.sex" (observed sex including unclassified), "maturity" (binary: 1 if mature; 0 if immature).
#' @examples
#' sim.dat<-sim_mat_data(nfemale = 30, nmale = 30, shapes = c(3, 3),
#'                      mat_parF = c(x50 = 0.5, xr = 0.1),
#'                      mat_parM = c(x50 = 0.4, xr = 0.1),
#'                      xcutoff = c(female = 0.5, male = 0.4))

sim_mat_data<-function(nfemale, nmale, shapes, mat_parF, mat_parM, xcutoff){
    ## female x
    x.vecF <- rbeta(nfemale, shape1 = shapes[1], shape2 = shapes[2])
    ## male x
    x.vecM <- rbeta(nmale, shape1 = shapes[1], shape2 = shapes[2])
    x.vec <- c(x.vecF, x.vecM)
    ## female maturity 
    matF <- rbinom(nfemale, size = 1, prob = plogis( 2 * log(3) / mat_parF["xr"]*(x.vecF - mat_parF["x50"])))
    ## male maturity
    matM <- rbinom(nfemale, size = 1, prob = plogis( 2 * log(3) / mat_parM["xr"]*(x.vecM - mat_parM["x50"])))
    mat <- c(matF, matM)
    ## true sex
    female.sex <- rep("female", times = nfemale)
    male.sex <- rep("male", times = nmale)
    true.sex <- rep(c("female","male"), times = c(nfemale, nmale))
    ## observed sex
    female.obs.sex <- ifelse(matF == 0 & x.vecF < xcutoff["female"], "unclassified", female.sex)
    male.obs.sex <- ifelse(matM == 0 & x.vecM < xcutoff["male"], "unclassified", male.sex)
    obs.sex <- c(female.obs.sex, male.obs.sex)
    ##
    df <- data.frame(x = x.vec, true.sex = true.sex, obs.sex = obs.sex, maturity = mat)
    return(df)
}
