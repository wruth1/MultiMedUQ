
#' Joint (asymptotic) covariance matrix of fixed and random effects' parameters from both models
#'
#' @param fit_Y Model for \code{Y}
#' @param fit_M Model for \code{M}
#'
#' @return A matrix of the joint covariance matrix of fixed and random effects' parameters from both models
#' @export
#'
#' @details
#'
#' Note: We assume that parameter estimators from the two models are independent. This may not be true in practice. See, e.g., Bauer, Preacher, & Gil (2006) for a method to incorporate inter-model dependence.
#'
#' Note: The actual covariance matrix of the parameter estimators requires further division by K
#'
all_pars_asymp_cov_mat = function(fit_Y, fit_M){
  Y_cov = merDeriv::vcov.glmerMod(fit_Y, full=TRUE, ranpar = "sd")
  M_cov = merDeriv::vcov.glmerMod(fit_M, full=TRUE, ranpar = "sd")

  Y_length = nrow(Y_cov)
  M_length = nrow(M_cov)
  joint_cov = matrix(0, nrow = M_length + Y_length, ncol = M_length + Y_length)
  joint_cov[1:Y_length, 1:Y_length] = Y_cov
  joint_cov[(Y_length + 1):(Y_length + M_length), (Y_length + 1):(Y_length + M_length)] = M_cov

  # Make output matrix symmetric
  joint_cov = (joint_cov + t(joint_cov))/2
  return(joint_cov)
}
