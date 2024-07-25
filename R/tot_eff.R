# Total effect on odds-ratio scale ####

#' Total effect on odds-ratio scale
#'
#' @param eta,zeta Linear predictors of M and Y models, respectively.
#' @param a1,b1,b2 Coefficients of X in M, M in Y and X in Y, respectively.
#' @param s1,s2 SD of linear predictor at x+1 in Y and M, respectively
#' @param s3,s4 SD of linear predictor at x in Y and M, respectively
#'
#' @return The total effect of X on Y on the odds-ratio scale.
#' @export
#'
tot_eff = function(eta, zeta, a1, b1, b2, s1, s2, s3, s4){
  (psi(zeta + b1 + b2, s1) * psi(eta + a1, s2) + psi(zeta + b2, s1) * (1 - psi(eta + a1, s2))) * (1 - psi(zeta, s3) + psi(eta, s4) * (psi(zeta, s3) - psi(zeta + b1, s3))) / (1 - psi(zeta + b2, s1) + psi(eta + a1, s2) * (psi(zeta + b2, s1) - psi(zeta + b1 + b2, s2))) / (psi(zeta + b1, s3) * psi(eta, s3) + psi(zeta, s3) * (1 - psi(eta, s4)))
}
