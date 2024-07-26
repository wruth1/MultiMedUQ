
#' Logistic and Inverse-Logistic (expit) transformations
#'
#' These functions are useful for transforming probabilities to the real line and back. logit(x) = log(x/(1-x)) and expit(x) = exp(x)/(1+exp(x)).
#'
#' @param x A numeric vector
#' @name logit_expit
#'
#' @return A numeric vector of the same length as x.
#' @export
#'
logit <- function(x){
  log(x/(1-x))
}

#' @rdname logit_expit
#' @export
#'
#' @examples
#' logit(c(0.1, 0.2, 0.3))
#' expit(c(0.1, 0.2, 0.3))
#'
#' logit(expit(c(0.1, 0.2, 0.3)))
#' expit(logit(c(0.1, 0.2, 0.3)))
expit <- function(x){
  exp(x)/(1+exp(x))
}



#' Convert a covariance matrix to a vector of SDs and correlations
#'
#' @param Sigma A covariance matrix
#'
#' @return theta, a vector of SDs and correlations. Order matches that of merDeriv.
#' @export
Sigma2theta <- function(Sigma){
  # Extract lower triangle of Sigma
  # Note: This orders parameters the same way as merDeriv
  theta = Sigma[lower.tri(Sigma, diag=T)]
  return(theta)
}

#' Convert a vector, theta, to its corresponding covariance matrix
#'
#' @param theta A vector of SDs and correlations. Order matches that of merDeriv.
#'
#' @return Sigma, a covariance matrix
#' @export
theta2Sigma <- function(theta){
  # Find the size of the covariance matrix which generated theta
  l = length(theta)
  n = (sqrt(1 + 8*l) - 1)/2
  if(!(n %% 1 == 0)){
    stop("theta is not the correct length for a covariance matrix")
  }

  # First, create matrix of SDs and correlations
  ## Create a matrix with zeros
  Sigma_raw = matrix(0, nrow = n, ncol = n)
  ## Fill in the lower triangle
  Sigma_raw[lower.tri(Sigma_raw, diag=T)] = theta
  ## Fill in the upper triangle
  Sigma_raw[upper.tri(Sigma_raw)] = t(Sigma_raw)[upper.tri(Sigma_raw)]

  # Next, convert to covariance matrix
  Sigma = Sigma_raw
  ## Fix the off-diagonals
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      Sigma[i,j] = Sigma[i,j] * (Sigma[i,i] * Sigma[j,j])
      Sigma[j,i] = Sigma[i,j]
    }
  }
  ## Fix the diagonals
  for(i in 1:n){
    Sigma[i,i] = Sigma[i,i]^2
  }

  return(Sigma)
}

