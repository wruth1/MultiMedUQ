

#' SD of REs evaluated at given vector of predictors
#'
#' @param pred_vec Vector of predictor values.
#' @param Sigma Covariance matrix of random effects.
#'
#' @name get_gamma
#'
#' @return The SD of the random effects evaluated at the given vector of predictors.
#' @export
Sigma2gamma <- function(pred_vec, Sigma){
  gamma = sqrt(pred_vec %*% Sigma %*% pred_vec)
  return(gamma)
}

#' @param theta Covariance parameters of random effects. See details.
#'
#' @rdname get_gamma
#'
#' @export
theta2gamma <- function(pred_vec, theta){
  Sigma = theta2Sigma(theta)
  return(Sigma2gamma(pred_vec, Sigma))
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
  mu_Y = b_Y[1] + x * b_Y[2] + w %*% b_Y[4:length(b_Y)]
  mu_M = b_M[1] + x_m * b_M[2] + w %*% b_M[3:length(b_M)]

  gamma_Y_1 = get_gamma(c(1, x, 1), theta_Y)  # SD of REs in Y model when M=1
  gamma_Y_0 = get_gamma(c(1, x, 0), theta_Y)  # SD of REs in Y model when M=0
  gamma_M = get_gamma(c(1,x_m), theta_M)      # SD of REs in M model

  psi_Y_1 = psi(mu_Y + b_Y[3], gamma_Y_1) # P(Y=1 | M=1)
  psi_Y_0 = psi(mu_Y, gamma_Y_0)          # P(Y=1 | M=0)
  psi_M_1 = psi(mu_M, gamma_M)            # P(M=1)
  psi_M_0 = psi(-mu_M, gamma_M)           # P(M=0)

  return(psi_Y_1 * psi_M_1 + psi_Y_0 * psi_M_0)
}




# Gradient of ENC ####

## First, gradients of the arguments to psi. ####

grad_mu_Y <- function(x, x_m, w, b_Y, theta_Y, b_M, theta_M){

  # Fixed effects for Y
  d_b_Y_0 = 1
  d_b_Y_X = x
  d_b_Y_M = 0
  d_B_Y_W = w
  d_b_Y = c(d_b_Y_0, d_b_Y_X, d_b_Y_M, d_B_Y_W)

  # Random effect parameters for Y
  d_theta_Y = rep(0, times = length(theta_Y))

  # Fixed effects for M
  d_b_M = rep(0, times = length(b_M))

  # Random effect parameters for M
  d_theta_M = rep(0, times = length(theta_M))

  return(c(d_b_Y, d_theta_Y, d_b_M, d_theta_M))
}


grad_mu_M <- function(x, x_m, w, b_Y, theta_Y, b_M, theta_M){

  # Fixed effects for Y
  d_b_Y = rep(0, times = length(b_Y))

  # Random effect parameters for Y
  d_theta_Y = rep(0, times = length(theta_Y))

  # Fixed effects for M
  d_b_M_0 = 1
  d_b_M_X = x_m
  d_B_M_W = w
  d_b_M = c(d_b_M_0, d_b_M_X, d_B_M_W)

  # Random effect parameters for M
  d_theta_M = rep(0, times = length(theta_M))

  return(c(d_b_Y, d_theta_Y, d_b_M, d_theta_M))
}

grad_b_Y_M <- function(x, x_m, w, b_Y, theta_Y, b_M, theta_M){

  # Fixed effects for Y
  d_b_Y_0 = 0
  d_b_Y_X = 0
  d_b_Y_M = 1
  d_B_Y_W = rep(0, times=length(w))
  d_b_Y = c(d_b_Y_0, d_b_Y_X, d_b_Y_M, d_B_Y_W)

  # Random effect parameters for Y
  d_theta_Y = rep(0, times = length(theta_Y))

  # Fixed effects for M
  d_b_M = rep(0, times = length(b_M))

  # Random effect parameters for M
  d_theta_M = rep(0, times = length(theta_M))

  return(c(d_b_Y, d_theta_Y, d_b_M, d_theta_M))
}


# Note the additional argument, m, at the front.
grad_gamma_Y <- function(m, x, x_m, w, b_Y, theta_Y, b_M, theta_M){

  output = 2 * get_gamma(c(1, x, m), theta_Y)

  # Fixed effects for Y
  d_b_Y = rep(0, times = length(b_Y))


  # Random effect parameters for Y
  d_s_Y_0 = 1
  d_cor_Y_0X = 2*x
  d_cor_Y_0M = 2*m
  d_s_Y_X = x^2
  d_cor_Y_XM = 2*x*m
  d_s_Y_M = m^2
  d_theta_Y = c(d_s_Y_0, d_cor_Y_0X, d_cor_Y_0M, d_s_Y_X, d_cor_Y_XM, d_s_Y_M)


  # Fixed effects for M
  d_b_M = rep(0, times = length(b_M))


  # Random effect parameters for M
  d_theta_M = rep(0, times = length(theta_M))


  return(c(d_b_Y, d_theta_Y, d_b_M, d_theta_M))

}

grad_gamma_M <- function(x, x_m, w, b_Y, theta_Y, b_M, theta_M){

  output = 2 * get_gamma(c(1, x_m), theta_M)

  # Fixed effects for Y
  d_b_Y = rep(0, times = length(b_Y))

  # Random effect parameters for Y
  d_theta_Y = rep(0, times = length(theta_Y))

  # Fixed effects for M
  d_b_M = rep(0, times = length(b_M))

  # Random effect parameters for M
  d_s_M_0 = 1
  d_cor_M_0X = 2*x_m
  d_s_M_X = x_m^2
  d_theta_M = c(d_s_M_0, d_cor_M_0X, d_s_M_X)

  return(c(d_b_Y, d_theta_Y, d_b_M, d_theta_M))

}

grad_psi_Y <- function(m, x, x_m, w, b_Y, theta_Y, b_M, theta_M){

  mu_Y = b_Y[1] + x * b_Y[2] + w %*% b_Y[4:length(b_Y)]
  gamma_Y = theta2gamma(c(1, x, m), theta_Y)

  grad_mu = grad_mu_Y(x, x_m, w, b_Y, theta_Y, b_M, theta_M)
  grad_gamma = grad_gamma_Y(m, x, x_m, w, b_Y, theta_Y, b_M, theta_M)

  return(d1_psi(mu_Y, gamma_Y) * grad_mu + d2_psi(mu_Y, gamma_Y) * grad_gamma)

}

grad_psi_M <- function(x, x_m, w, b_Y, theta_Y, b_M, theta_M){

  mu_M = b_M[1] + x_m * b_M[2] + w %*% b_M[3:length(b_M)]
  gamma_M = theta2gamma(c(1, x_m), theta_M)

  grad_mu = grad_mu_M(x, x_m, w, b_Y, theta_Y, b_M, theta_M)
  grad_gamma = grad_gamma_M(x, x_m, w, b_Y, theta_Y, b_M, theta_M)

  return(d1_psi(mu_M, gamma_M) * grad_mu + d2_psi(mu_M, gamma_M) * grad_gamma)

}



############## THE ABOVE FUNCTIONS NEED TO BE TESTED!!!!!!


#
#
# # Fixed effects for Y
# d_b_Y = rep(0, times = length(b_Y))
# d_b_Y_0 = 0
# d_b_Y_X = 0
# d_b_Y_M = 0
# d_B_Y_W = rep(0, times=length(w))
# d_b_Y = c(d_b_Y_0, d_b_Y_X, d_b_Y_M, d_B_Y_W)
#
# # Random effect parameters for Y
# d_theta_Y = rep(0, times = length(theta_Y))
# d_s_Y_0 = 0
# d_cor_Y_0X = 0
# d_cor_Y_0M = 0
# d_s_Y_X = 0
# d_cor_Y_XM = 0
# d_s_Y_M = 0
# d_theta_Y = c(d_s_Y_0, d_cor_Y_0X, d_cor_Y_0M, d_s_Y_X, d_cor_Y_XM, d_s_Y_M)
#
#
# # Fixed effects for M
# d_b_M = rep(0, times = length(b_M))
# d_b_M_0 = 0
# d_b_M_X = 0
# d_B_M_W = rep(0, times=length(w))
# d_b_M = c(d_b_M_0, d_b_M_X, d_B_M_W)
#
# # Random effect parameters for M
# d_theta_M = rep(0, times = length(theta_M))
# d_s_M_0 = 0
# d_cor_M_0X = 0
# d_s_M_X = 0
# d_theta_M = c(d_s_M_0, d_cor_M_0X, d_s_M_X)
#
# return(c(d_b_Y, d_theta_Y, d_b_M, d_theta_M))
