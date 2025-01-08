
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



#' Convert an SD-corr matrix to a vector of SDs and correlations
#'
#' @param Sigma An SD-corr matrix
#'
#' @return theta, a vector of SDs and correlations. Order matches that of merDeriv (e.g. SD, cor, SD for a 2x2).
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


# Get the number of REs used to construct theta
theta2num_REs <- function(theta){
  l = length(theta)
  n = (sqrt(1 + 8*l) - 1)/2
  return(n)
}


# Re-arrange a vector of SDs followed by correlations to my theta order
## I.e. SD, corr, corr, ..., SD, corr, corr, ..., SD
## Output: a vector of indices, inds, s.t. sd_corr[inds] == theta
SD_corr2theta_indices <- function(num_pars = NULL, num_vars = NULL){
  if(!is.null(num_pars) & is.null(num_vars)){
    num_vars = (sqrt(1 + 8*num_pars) - 1 ) / 2
    if(num_vars %% 1 != 0){
      stop("Incorrect number of parameters")
    }
  } else if(is.null(num_pars) & !is.null(num_vars)){
    num_pars = num_vars * (num_vars + 1) / 2
  } else{
    stop("Must specify exactly one of num_pars or num_vars")
  }

  # Container to store indices
  output = rep(0, times = num_pars)


  # Destination indices of SDs (formula derived by hand)
  ## Also useful as anchors for correlation destination indices
  inds_SD_dest = ((1:num_vars) - 1) * (1 + num_vars - (1:num_vars)/2) + 1

  # Fill-in indices
  for(i in seq_len(num_vars)){
    ### SDs
    this_ind_SD_dest = inds_SD_dest[i]
    output[this_ind_SD_dest] = i

    ### Correlations
    this_num_corrs = num_vars - i    # Number of correlations for this variable
    num_corrs_so_far = (i-1)*(num_vars - i/2)    # Number of correlations already accounted for (this helps us index d_corr)
    for(j in seq_len(this_num_corrs)){
        this_ind_corr_dest = this_ind_SD_dest + j   # Index in my theta of current correlation

        this_ind_corr_origin = num_corrs_so_far + j    # Index in d_corr of current correlation

        output[this_ind_corr_dest] = num_vars + this_ind_corr_origin
    }
  }

  return(output)
}







#' Extract fitted parameters of interest from an lme4 model.
#'
#' @param fit An lme4 model
#' @param format The format of the output. Can be "list" or "vector".
#'
#' @details
#' Note: It is unnecessary to specify which random effects are required, since we extract exactly what is contained in the fitted \code{lme4} object.
#'
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


sds_corrs2theta <- function(sd_corr_vec){
  num_vars = (sqrt(1 + 8*length(sd_corr_vec)) - 1 ) / 2

  sds = sd_corr_vec[1:num_vars]
  corrs = sd_corr_vec[(num_vars+1):length(sd_corr_vec)]
  len_corrs = length(corrs)


  # Build SD-Correlation matrix
  sd_corr_mat = matrix(0, nrow = num_vars, ncol = num_vars)

  for(i in 1:num_vars){
    for(j in 1:num_vars){
      if(i == j) {
        sd_corr_mat[i,j] = sds[i]
      } else {
         sd_corr_mat[i,j] = corrs[(i - 1) + (j - 1)]
      }
    }
  }

  # Convert to theta
  theta = Sigma2theta(sd_corr_mat)
  return(theta)
}


#' Extract the parameters (in our order) from a glmmTMB object
#'
#' @param fit A GLMM fit using glmmTMB
#' @param format The format of the output. Can be "list" or "vector".
#'
#' @return A vector, theta, of model parameters. Order is fixed effects, then covariance parameters. The latter is organized as, e.g., SD, corr, corr, SD, corr, SD for a 3x3 covariance matrix.
#' @export
#'
get_model_pars_TMB <- function(fit, format="list"){
  b = glmmTMB::fixef(fit)[[1]]


  theta = broom.mixed::tidy(fit) %>% dplyr::filter(effect == "ran_pars") %>% dplyr::pull(estimate) %>% sds_corrs2theta
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
#' \item "All": All REs
#' \item "Y.All": All REs for Y
#' \item "M.All": All REs for M
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



  if(identical(RE_input, "All")){ # All REs
    all_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")
  } else{
    if("Y.All" %in% RE_input){ # All Y REs
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

    if("M.All" %in% RE_input){ # All M REs
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



# Number of parameters required to encode the joint distribution of the specified list of REs
## Note: This is just the sum of the first r positive integers, where r is the number of REs.
REs2theta_length = function(REs){
  num_REs = length(REs)

  num_REs * (num_REs + 1) / 2
}

# Similar to the previous function, but this one takes the number of random effects.
num_REs2theta_length = function(num_REs){
  num_REs * (num_REs + 1) / 2
}



# Return the number of REs for Y
num_Y_REs = function(which_REs){
  RE_names = expand_REs(which_REs)

  # Extract REs starting with "Y"
  # Note: The following will fail if there are no Y REs
  Y_REs = RE_names[grepl("^Y", RE_names)]

  return(length(Y_REs))
}


# Return the number of REs for M
num_M_REs = function(which_REs){
  RE_names = expand_REs(which_REs)

  # Extract REs starting with "M"
  # Note: The following will fail if there are no M REs
  M_REs = RE_names[grepl("^M", RE_names)]

  return(length(M_REs))
}
