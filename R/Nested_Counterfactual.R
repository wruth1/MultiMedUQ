

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
  return(as.numeric(gamma))
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

#' Vectors used to compute gamma
#'
#' @param x,m Values of X and M at which gamma is evaluated.
#' @param which_REs Which random effects to include in the calculation. Default is all. Shorthands are available. See details.
#'
#' @name gamma_vecs
#'
#' @return A vector to be used in the quadratic form defining gamma.
#'
#' @details
#' The following shorthands for random effects are available:
#' \itemize{
#' \item "all": All REs
#' \item "Y.all": All REs for Y
#' \item "M.all": All REs for M
#' }
#' Additionally, individual REs can be specified:
#' \itemize{
#' \item "Y.Int": Intercept for Y
#' \item "Y.X": Slope for X in Y
#' \item "Y.M": Slope for M in Y
#' \item "M.Int": Intercept for M
#' \item "M.X": Slope for M
#' }
#'
Y_vec_gamma <- function(x, m, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){
  RE_names = expand_REs(which_REs)

  Y_vec = c()
  if("Y.Int" %in% RE_names){
    Y_vec = c(Y_vec, 1)
  }
  if("Y.X" %in% RE_names){
    Y_vec = c(Y_vec, x)
  }
  if("Y.M" %in% RE_names){
    Y_vec = c(Y_vec, m)
  }

  return(Y_vec)
}

#' @rdname gamma_vecs
M_vec_gamma <- function(x, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){
  RE_names = expand_REs(which_REs)

  M_vec = c()
  if("M.Int" %in% RE_names){
    M_vec = c(M_vec, 1)
  }
  if("M.X" %in% RE_names){
    M_vec = c(M_vec, x)
  }

  return(M_vec)
}


#' Expected nested counterfactual. Binary Y, binary M.
#'
#' @param x Level of exposure variable, \eqn{X}.
#' @param x_m Level of exposure for mediator. I.e. Mediator set to \eqn{M_{x\_m}}.
#' @param w Level of covariates, \eqn{W}.
#' @param b_Y,b_M Coefficient vectors for \eqn{Y}-model and \eqn{M}-model, respectively.
#' @param theta_Y,theta_M Covariance parameters of random effects in \eqn{Y}-model and \eqn{M}-model, respectively. See details.
#' @param which_REs Which random effects to include in the calculation. Default is all. Shorthands are available. See details.
#'
#' @return The expected nested counterfactual of \eqn{Y} with \eqn{X = x} and \eqn{M = M_{x\_m}}.
#' @export
#'
#' @details
#' Contents of \code{b_Y} are \code{(b_Y_0, b_Y_X, b_Y_M, B_Y_W)}. Contents of \code{b_M} are \code{(b_M_0, b_M_X, B_M_W)}.
#'
#' Contents of \code{theta_Y} are \code{(s_Y_0, cor_Y_0X, cor_Y_0M, s_Y_X, cor_Y_XM, s_Y_M)}. Contents of \code{theta_M} are \code{(s_M_0, cor_M_0X, s_M_X)}.
#'
#' The following shorthands for random effects are available:
#' \itemize{
#' \item "all": All REs
#' \item "Y.all": All REs for Y
#' \item "M.all": All REs for M
#' }
#' Additionally, individual REs can be specified:
#' \itemize{
#' \item "Y.Int": Intercept for Y
#' \item "Y.X": Slope for X in Y
#' \item "Y.M": Slope for M in Y
#' \item "M.Int": Intercept for M
#' \item "M.X": Slope for M
#' }
#'
#'
ENC <- function(x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){

  RE_names = expand_REs(which_REs)


  mu_Y = as.numeric(b_Y[1] + x * b_Y[2] + w %*% b_Y[4:length(b_Y)])
  mu_M = as.numeric(b_M[1] + x_m * b_M[2] + w %*% b_M[3:length(b_M)])

  # SD of REs in Y model when M=1
  Y_vec_1 = Y_vec_gamma(x, 1, which_REs)
  gamma_Y_1 = theta2gamma(Y_vec_1, theta_Y)

  # SD of REs in Y model when M=0
  Y_vec_0 = Y_vec_gamma(x, 0, which_REs)
  gamma_Y_0 = theta2gamma(Y_vec_0, theta_Y)

  # SD of REs in M model
  M_vec = M_vec_gamma(x_m, which_REs)
  gamma_M = theta2gamma(M_vec, theta_M)

  psi_Y_1 = psi(mu_Y + b_Y[3], gamma_Y_1) # P(Y=1 | M=1)
  psi_Y_0 = psi(mu_Y, gamma_Y_0)          # P(Y=1 | M=0)
  psi_M_1 = psi(mu_M, gamma_M)            # P(M=1)
  psi_M_0 = psi(-mu_M, gamma_M)           # P(M=0)

  return(psi_Y_1 * psi_M_1 + psi_Y_0 * psi_M_0)
}


#' Expected nested counterfactuals for all configurations of \eqn{X} and \eqn{X_M}.
#'
#' @param w Level of covariates, \eqn{W}.
#' @param b_Y,b_M Coefficient vectors for \eqn{Y}-model and \eqn{M}-model, respectively.
#' @param theta_Y,theta_M Covariance parameters of random effects in \eqn{Y}-model and \eqn{M}-model, respectively. See details.
#' @param which_REs Which random effects to include in the calculation. Default is all. Shorthands are available. See details.
#'
#' @details
#' The following shorthands for random effects are available:
#' \itemize{
#' \item "all": All REs
#' \item "Y.all": All REs for Y
#' \item "M.all": All REs for M
#' }
#' Additionally, individual REs can be specified:
#' \itemize{
#' \item "Y.Int": Intercept for Y
#' \item "Y.X": Slope for X in Y
#' \item "Y.M": Slope for M in Y
#' \item "M.Int": Intercept for M
#' \item "M.X": Slope for M
#' }
#'
#'
#' @return A vector of expected nested counterfactuals for all configurations of \eqn{X} and \eqn{X_M}. Order of output is \code{ENC(1,1), ENC(1,0), ENC(0,1), ENC(0,0)}.
#' @export
#'
all_ENCs <- function(w, b_Y, theta_Y, b_M, theta_M, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){
  ENC_11 = ENC(1, 1, w, b_Y, theta_Y, b_M, theta_M, which_REs)
  ENC_10 = ENC(1, 0, w, b_Y, theta_Y, b_M, theta_M, which_REs)
  ENC_01 = ENC(0, 1, w, b_Y, theta_Y, b_M, theta_M, which_REs)
  ENC_00 = ENC(0, 0, w, b_Y, theta_Y, b_M, theta_M, which_REs)

  return(c(ENC_11, ENC_10, ENC_01, ENC_00))
}



# Gradient of ENC ####

## First, gradients of the arguments to psi. ####

### Means

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


### SDs

# Note the additional argument, m, at the front.
grad_gamma_Y <- function(m, x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){

  Y_vec = Y_vec_gamma(x, m, which_REs)
  divisor = 2 * theta2gamma(Y_vec, theta_Y)

  # Fixed effects for Y
  d_b_Y = rep(0, times = length(b_Y))


  # Random effect parameters for Y
  ## Check which random effects are included. Exclude any that aren't from future calculations
  s_Y_0 = ifelse("Y.Int" %in% which_REs, theta_Y[1], 0)                     # SD of intercept
  cor_Y_0X = ifelse(all(c("Y.Int", "Y.X") %in% which_REs), theta_Y[2], 0)   # Corr between intercept and X coeff
  cor_Y_0M = ifelse(all(c("Y.Int", "Y.M") %in% which_REs), theta_Y[3], 0)   # Corr between intercept and M coeff
  s_Y_X = ifelse("Y.X" %in% which_REs, theta_Y[4], 0)                       # SD of X coeff
  cor_Y_XM = ifelse(all(c("Y.X", "Y.M") %in% which_REs), theta_Y[5], 0)     # Corr between X and M coeffs
  s_Y_M = ifelse("Y.M" %in% which_REs, theta_Y[6], 0)                       # SD of M coeff

  d_s_Y_0 = ifelse("Y.Int" %in% which_REs, 2 * s_Y_0 + 2* x * s_Y_X * cor_Y_0X + 2*m * s_Y_M * cor_Y_0M, NULL)
  d_cor_Y_0X = ifelse(all(c("Y.Int", "Y.X") %in% which_REs), 2*x * s_Y_0 * s_Y_X, NULL)
  d_cor_Y_0M = ifelse(all(c("Y.Int", "Y.M") %in% which_REs), 2*m * s_Y_0 * s_Y_M, NULL)
  d_s_Y_X = ifelse("Y.X" %in% which_REs, 2 * x^2 * s_Y_X + 2 * x * s_Y_0 * cor_Y_0X + 2 * x * m * s_Y_M * cor_Y_XM, NULL)
  d_cor_Y_XM = ifelse(all(c("Y.X", "Y.M") %in% which_REs), 2*x*m * s_Y_X * s_Y_M, NULL)
  d_s_Y_M = ifelse("Y.M" %in% which_REs, 2 * m^2 * s_Y_M + 2 * m * s_Y_0 * cor_Y_0M + 2 * x * m * s_Y_X * cor_Y_XM, NULL)

  d_theta_Y = c(d_s_Y_0, d_cor_Y_0X, d_cor_Y_0M, d_s_Y_X, d_cor_Y_XM, d_s_Y_M)

  if(length(d_theta_Y) > length(theta_Y)){
    stop("Gradients requested for random effects not present in theta_Y.")
  } else if(length(d_theta_Y) < length(theta_Y)){
    stop("Random effects present in theta_Y whose gradients are not requested.")
  }


  # Fixed effects for M
  d_b_M = rep(0, times = length(b_M))


  # Random effect parameters for M
  d_theta_M = rep(0, times = length(theta_M))


  return(c(d_b_Y, d_theta_Y, d_b_M, d_theta_M) / divisor)

}

grad_gamma_M <- function(x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){

  M_vec = M_vec_gamma(x_m, which_REs)
  divisor = 2 * theta2gamma(M_vec, theta_M)

  # Fixed effects for Y
  d_b_Y = rep(0, times = length(b_Y))

  # Random effect parameters for Y
  d_theta_Y = rep(0, times = length(theta_Y))

  # Fixed effects for M
  d_b_M = rep(0, times = length(b_M))

  # Random effect parameters for M
  s_M_0 = ifelse("M.Int" %in% which_REs, theta_M[1], 0)
  cor_M_0X = ifelse(all(c("M.Int", "M.X") %in% which_REs), theta_M[2], 0)
  s_M_X = ifelse("M.X" %in% which_REs, theta_M[3], 0)

  d_s_M_0 = ifelse("M.Int" %in% which_REs, 2 * s_M_0 + 2*x_m * s_M_X * cor_M_0X, NULL)
  d_cor_M_0X = ifelse(all(c("M.Int", "M.X") %in% which_REs), 2*x_m * s_M_0 * s_M_X, NULL)
  d_s_M_X = ifelse("M.X" %in% which_REs, 2 * x_m^2 * s_M_X + 2 * x_m * s_M_0 * cor_M_0X, NULL)

  d_theta_M = c(d_s_M_0, d_cor_M_0X, d_s_M_X)

  if(length(d_theta_M) > length(theta_M)){
    stop("Gradients requested for random effects not present in theta_M.")
  } else if(length(d_theta_M) < length(theta_M)){
    stop("Random effects present in theta_M whose gradients are not requested.")
  }

  return(c(d_b_Y, d_theta_Y, d_b_M, d_theta_M) / divisor)

}




## Now, the gradient of psi ####

### For a Y-probability
#' Gradient of psi used to compute a probability for Y or M
#'
#' @param m Value of M
#' @param x,x_m values of X in the models for Y and M respectively
#' @param w Vector of covariates
#' @param b_Y,b_M Fixed effects for Y and M models respectively
#' @param theta_Y,theta_M SDs and correlations for Y and M models respectively
#' @param which_REs Which random effects to include in the calculation. Default is all. Shorthands are available. See details.
#' @name grad_psi
#'
#' @details
#' The following shorthands for random effects are available:
#' \itemize{
#' \item "all": All REs
#' \item "Y.all": All REs for Y
#' \item "M.all": All REs for M
#' }
#' Additionally, individual REs can be specified:
#' \itemize{
#' \item "Y.Int": Intercept for Y
#' \item "Y.X": Slope for X in Y
#' \item "Y.M": Slope for M in Y
#' \item "M.Int": Intercept for M
#' \item "M.X": Slope for M
#' }
#'
#' @return Gradient of the function \code{psi(mu, sigma)} with respect to \code{b_Y}, \code{theta_Y}, \code{b_M}, and \code{theta_M}.
#' @export
grad_psi_Y <- function(m, x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){

  mu_Y = as.numeric(b_Y[1] + x * b_Y[2] + w %*% b_Y[4:length(b_Y)])
  Y_vec = Y_vec_gamma(x, m, which_REs)
  gamma_Y = theta2gamma(Y_vec, theta_Y)

  grad_mu = grad_mu_Y(x, x_m, w, b_Y, theta_Y, b_M, theta_M)
  grad_M_contrib = m * grad_b_Y_M(x, x_m, w, b_Y, theta_Y, b_M, theta_M)
  grad_gamma = grad_gamma_Y(m, x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs)

  return(d1_psi(mu_Y + m*b_Y[3], gamma_Y) * (grad_mu +  grad_M_contrib) + d2_psi(mu_Y + m*b_Y[3], gamma_Y) * grad_gamma)

}

### For an M-probability

#' @export
#' @rdname grad_psi
grad_psi_M <- function(m, x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){

  mu_M = as.numeric(b_M[1] + x_m * b_M[2] + w %*% b_M[3:length(b_M)])
  M_vec = M_vec_gamma(x_m, which_REs)
  gamma_M = theta2gamma(M_vec, theta_M)

  grad_mu = (2*m - 1) * grad_mu_M(x, x_m, w, b_Y, theta_Y, b_M, theta_M)
  grad_gamma = grad_gamma_M(x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs)

  return(d1_psi(mu_M, gamma_M) * grad_mu + d2_psi(mu_M, gamma_M) * grad_gamma)

}


#' Gradient of expected nested counterfactual. Binary Y, binary M.
#'
#' @param x Level of exposure variable, \eqn{X}.
#' @param x_m Level of exposure for mediator. I.e. Mediator set to \eqn{M_{x\_m}}.
#' @param w Level of covariates, \eqn{W}.
#' @param b_Y,b_M Coefficient vectors for \eqn{Y}-model and \eqn{M}-model, respectively.
#' @param theta_Y,theta_M Covariance parameters of random effects in \eqn{Y}-model and \eqn{M}-model, respectively. See details.
#' @param which_REs Which random effects to include in the calculation. Default is all. Shorthands are available. See details.
#'
#' @return Gradient of the expected nested counterfactual of \eqn{Y} with \eqn{X = x} and \eqn{M = M_{x\_m}}.
#' @export
#'
#' @details
#' Contents of \code{b_Y} are \code{(b_Y_0, b_Y_X, b_Y_M, B_Y_W)}. Contents of \code{b_M} are \code{(b_M_0, b_M_X, B_M_W)}.
#'
#' Contents of \code{theta_Y} are \code{(s_Y_0, cor_Y_0X, cor_Y_0M, s_Y_X, cor_Y_XM, s_Y_M)}. Contents of \code{theta_M} are \code{(s_M_0, cor_M_0X, s_M_X)}.
#'
#'
#' The following shorthands for random effects are available:
#' \itemize{
#' \item "all": All REs
#' \item "Y.all": All REs for Y
#' \item "M.all": All REs for M
#' }
#' Additionally, individual REs can be specified:
#' \itemize{
#' \item "Y.Int": Intercept for Y
#' \item "Y.X": Slope for X in Y
#' \item "Y.M": Slope for M in Y
#' \item "M.Int": Intercept for M
#' \item "M.X": Slope for M
#' }
#'
#'
grad_ENC <- function(x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){

  mu_Y = as.numeric(b_Y[1] + x * b_Y[2] + w %*% b_Y[4:length(b_Y)])
  mu_M = as.numeric(b_M[1] + x_m * b_M[2] + w %*% b_M[3:length(b_M)])

  # SD of REs in Y model when M=1
  Y_vec_1 = Y_vec_gamma(x, 1, which_REs)
  gamma_Y_1 = theta2gamma(Y_vec_1, theta_Y)

  # SD of REs in Y model when M=0
  Y_vec_0 = Y_vec_gamma(x, 0, which_REs)
  gamma_Y_0 = theta2gamma(Y_vec_0, theta_Y)

  # SD of REs in M model
  M_vec = M_vec_gamma(x_m, which_REs)
  gamma_M = theta2gamma(M_vec, theta_M)


  psi_Y_1 = psi(mu_Y + b_Y[3], gamma_Y_1) # P(Y=1 | M=1)
  psi_Y_0 = psi(mu_Y, gamma_Y_0)          # P(Y=1 | M=0)
  psi_M_1 = psi(mu_M, gamma_M)            # P(M=1)
  psi_M_0 = psi(-mu_M, gamma_M)           # P(M=0)


  grad_psi_Y_1 = grad_psi_Y(1, x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs)
  grad_psi_Y_0 = grad_psi_Y(0, x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs)

  grad_psi_M_1 = grad_psi_M(1, x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs)
  grad_psi_M_0 = grad_psi_M(0, x, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs)


  output = psi_Y_1 * grad_psi_M_1 + grad_psi_Y_1 * psi_M_1 + psi_Y_0 * grad_psi_M_0 + grad_psi_Y_0 * psi_M_0
  return(output)
}





#### Covariance matrix of all ENC estimates ####
## I.e. All combinations of levels of X and X_M (there are only two of each)
## Organization is (X,X_M) = (1,1), (1,0), (0,1), (0,0)


#' "Jacobian" of expected nested counterfactuals. I.e. Gradients for every combination of \eqn{X, X_M}.
#'
#' @param w Level of covariates, \eqn{W}.
#' @param b_Y,b_M Coefficient vectors for \eqn{Y}-model and \eqn{M}-model, respectively.
#' @param theta_Y,theta_M Covariance parameters of random effects in \eqn{Y}-model and \eqn{M}-model, respectively. See details.
#' @param which_REs Which random effects to include in the calculation. Default is all. Shorthands are available. See details.
#'
#' @name Jacob_ENC
#'
#' @details
#' The following shorthands for random effects are available:
#' \itemize{
#' \item "all": All REs
#' \item "Y.all": All REs for Y
#' \item "M.all": All REs for M
#' }
#' Additionally, individual REs can be specified:
#' \itemize{
#' \item "Y.Int": Intercept for Y
#' \item "Y.X": Slope for X in Y
#' \item "Y.M": Slope for M in Y
#' \item "M.Int": Intercept for M
#' \item "M.X": Slope for M
#' }
#'
#'
#' @return Gradients for every combination of \eqn{X} and \eqn{X_M}, organized as a 4-by-many matrix. Order of \eqn{(X, X_M)} levels is (1,1), (1,0), (0,1), (0,0).
#' @export
#'
Jacob_ENC_pars <- function(w, b_Y, theta_Y, b_M, theta_M, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){
  grad_ENC_11 = grad_ENC(1, 1, w, b_Y, theta_Y, b_M, theta_M, which_REs)
  grad_ENC_10 = grad_ENC(1, 0, w, b_Y, theta_Y, b_M, theta_M, which_REs)
  grad_ENC_01 = grad_ENC(0, 1, w, b_Y, theta_Y, b_M, theta_M, which_REs)
  grad_ENC_00 = grad_ENC(0, 0, w, b_Y, theta_Y, b_M, theta_M, which_REs)

  return(rbind(grad_ENC_11, grad_ENC_10, grad_ENC_01, grad_ENC_00))
}

#' @rdname Jacob_ENC
#' @param fit_Y,fit_M Fitted models for Y and M.
#'
#' @export
Jacob_ENC_models <- function(w, fit_Y, fit_M, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){
  info_Y = get_model_pars(fit_Y)
  info_M = get_model_pars(fit_M)

  b_Y = info_Y[["b"]]
  theta_Y = info_Y[["theta"]]
  b_M = info_M[["b"]]
  theta_M = info_M[["theta"]]

  return(Jacob_ENC_pars(w, b_Y, theta_Y, b_M, theta_M))
}



#
#' Covariance matrix of all ENC estimates
#'
#' @param w Level of covariates, \eqn{W}.
#' @param fit_Y,fit_M Fitted models for Y and M.
#' @param which_REs Which random effects to include in the calculation. Default is all. Shorthands are available. See details.
#'
#' @details
#' Note: Uses the \eqn{K}-adjusted covariance matrix, not the asymptotic covariance matrix.
#'
#' The following shorthands for random effects are available:
#' \itemize{
#' \item "all": All REs
#' \item "Y.all": All REs for Y
#' \item "M.all": All REs for M
#' }
#' Additionally, individual REs can be specified:
#' \itemize{
#' \item "Y.Int": Intercept for Y
#' \item "Y.X": Slope for X in Y
#' \item "Y.M": Slope for M in Y
#' \item "M.Int": Intercept for M
#' \item "M.X": Slope for M
#' }
#'
#' @return The covariance matrix of all ENC estimates. Order of \eqn{(X, X_M)} levels is (1,1), (1,0), (0,1), (0,0).
#' @export
all_covs_ENC <- function(w, fit_Y, fit_M, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){
  Sigma = all_pars_cov_mat(fit_Y, fit_M)
  Jacob = Jacob_ENC_models(w, fit_Y, fit_M, which_REs)
  return(Jacob %*% Sigma %*% t(Jacob))
}

# Run the first few lines of "test-Reg_Par_Covs.R" to get the arguments for the following function.
# Q = all_covs_ENC(w, fit_Y, fit_M)
