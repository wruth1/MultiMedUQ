
#' Joint covariance matrix of fixed and random effects' parameters from both models
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
#' Note: I thought that the actual covariance matrix of the parameter estimators requires further division by K. However, numerical experimentation showed that this is not the case. The output of this function scales (approximately) linearly with the number of groups.
#'
#' Note: It is unnecessary to specify which random effects are included because we extract all present REs from the fitted \code{lme4} model.
#'
all_pars_cov_mat = function(fit_Y, fit_M){
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


#? The below function was meant to scale an asymptotic covariance matrix obtained from merDeriv::vcov.glmerMod() to match the observed number of groups. However, it turns out that vcov() already does this and I was double-adjusting. Oopie.
#? I'm going to keep the function here as a reminder to not do that again.
# #' Joint (sample-size adjusted) covariance matrix of fixed and random effects' parameters from both models
# #'
# #' @param fit_Y Model for \code{Y}
# #' @param fit_M Model for \code{M}
# #'
# #' @return A matrix of the joint covariance matrix of fixed and random effects' parameters from both models. Obtained from asymptotic theory, then re-scaled to match the observed number of groups.
# #' @export
# #'
# #' @details
# #'
# #' Note: We assume that parameter estimators from the two models are independent. This may not be true in practice. See, e.g., Bauer, Preacher, & Gil (2006) for a method to incorporate inter-model dependence.
# #'
# #' Note: It is unnecessary to specify which random effects are included because we extract all present REs from the fitted \code{lme4} model.
# #'
# all_pars_cov_mat = function(fit_Y, fit_M){
#   asymp_cov = all_pars_asymp_cov_mat(fit_Y, fit_M)

#   K_Y = lme4::ngrps(fit_Y)
#   K_M = lme4::ngrps(fit_M)
#   if(K_Y != K_M){
#     stop("Number of groups in models for Y and M do not match.")
#   }
#   K = K_Y

#   return(asymp_cov/K)
# }
