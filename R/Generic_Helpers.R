
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
#' @param theta A vector of SDs and correlations. Order matches that of merDeriv (i.e. SD, cor, SD for a 2x2).
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

  if(n == 1) return(matrix(theta^2))

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



#' Extract fitted parameters of interest from an lme4 model.
#'
#' @param fit An lme4 model
#' @param format The format of the output. Can be "list" or "vector".
#'
#' @return A list with elements b (fixed effects) and theta (RE SDs and correlations), or a vector with b followed by theta.
#' @export
#'
get_model_pars <- function(fit, format="list"){
  b = lme4::fixef(fit)

  # Estimated raneff covariance matrix
  # warning("Confirm that order of theta is correct.")
  info_cov = as.data.frame(lme4::VarCorr(fit))
  info_cov_sort = info_cov[order(info_cov$var1),]
  theta = info_cov_sort$sdcor

  if(format == "list"){
    return(list(b=b, theta=theta))
  } else if(format == "vector"){
    return(c(b, theta))
  } else{
    stop("Invalid format")
  }
}





#' Expand shorthand notation for random effects in models for Y and M
#'
#' @param RE_input A character vector of REs to include. May contain shorthands. See details.
#'
#'
#' @return A character vector containing names of the REs to include in the models.
#' @export
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
expand_REs <- function(RE_input){
  all_REs = c()

  if(identical(RE_input, "all")){ # All REs
    all_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")
  } else{
    if("Y.all" %in% RE_input){ # All Y REs
      all_REs = c(all_REs, "Y.Int", "Y.X", "Y.M")
    } else{ # Some Y REs
      if("Y.Int" %in% RE_input){
        all_REs = c(all_REs, "Y.Int")
      }
      if("Y.X" %in% RE_input){
        all_REs = c(all_REs, "Y.X")
      }
      if("Y.M" %in% RE_input){
      all_REs = c(all_REs, "Y.M")
      }
    }

    if("M.all" %in% RE_input){ # All M REs
      all_REs = c(all_REs, "M.Int", "M.X")
    } else{ # Some M REs
      if("M.Int" %in% RE_input){
        all_REs = c(all_REs, "M.Int")
      }
      if("M.X" %in% RE_input){
        all_REs = c(all_REs, "M.X")
      }
    }
  }

  return(unique(all_REs))
}
