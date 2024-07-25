

#' SD of REs evaluated at given vector of predictors
#'
#' @param pred_vec Vector of predictor values.
#' @param Sigma Covariance matrix of random effects.
#'
#' @return The SD of the random effects evaluated at the given vector of predictors.
#' @export
get_gamma <- function(pred_vec, Sigma){
  gamma = sqrt(pred_vec %*% Sigma %*% pred_vec)
  return(gamma)
}


#' Expected nested counterfactual. Binary Y, binary M.
#'
#' @param x Level of exposure variable, \eqn{X}.
#' @param x_m Level of exposure for mediator. I.e. Mediator set to \eqn{M_{x\_m}}.
#' @param w Level of covariates, \eqn{W}.
#' @param b_Y,b_M Coefficient vectors for \eqn{Y}-model and \eqn{M}-model, respectively.
#' @param theta_Y,theta_M Covariance parameters of random effects in \eqn{Y}-model and \eqn{M}-model, respectively. See details.
#'
#' @return The expected nested counterfactual of \eqn{Y} with \eqn{X = x} and \eqn{M = M_{x\_m}}.
#' @export
#'
#' @details
#' Contents of \code{b_Y} are \code{(b_Y_0, b_Y_X, b_Y_M, B_Y_W)}. Contents of \code{b_M} are \code{(b_M_0, b_M_X, B_M_W)}.
#'
#' Contents of \code{theta_Y} are \code{(s_Y_0, cor_Y_0X, cor_Y_0M, s_Y_X, cor_Y_XM, s_Y_M)}. Contents of \code{theta_M} are \code{(s_M_0, cor_M_0X, s_M_X)}.
#'
#'
ENC <- function(x, x_m, w, b_Y, theta_Y, b_M, theta_M){
  mu_Y = b_Y[1] + x * b_Y[2] + w %*% b_Y[4]
  mu_M = b_M[1] + x_m * b_M[2] + w %*% b_M[3]

  gamma_Y_2 = get_gamma(c(1, x, 1), theta_Y)
  gamma_Y_1 = get_gamma(c(1, x, 0), theta_Y)
  gamma_M = get_gamma(c(1,x_m), theta_M)

  psi_Y_2 = psi(mu_Y + b_Y[3], gamma_Y_2)
  psi_Y_1 = psi(mu_Y, gamma_Y_1)
  psi_M_2 = psi(mu_M, gamma_M)
  psi_M_1 = psi(-mu_M, gamma_M)

  return(psi_Y_2 * psi_M_2 + psi_Y_1 * psi_M_1)
}


