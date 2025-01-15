

# --------------- Extract fitted GLMM parameters (SD and corr) --------------- #

#' Extract fitted internal parameters from a glmmTMB object
#'
#' Uses the \code{glmmTMB} parameterization of log-SDs and Cholesky factors from the scaled correlation matrix (see the [TMB documentation](https://kaskr.github.io/adcomp/classdensity_1_1UNSTRUCTURED__CORR__t.html) for more details).
#'
#' @param fit_TMB A model fit using \code{glmmTMB}
#'
#'
#'
#' @returns A list containing fixed effects (\code{TMB_FE_pars}), log-standard deviations (\code{TMB_SD_pars}), and Cholesky factors from the (scaled) correlation matrix (\code{TMB_corr_pars}).
#' @export
TMB_pars_list <- function(fit_TMB) {
    TMB_pars = fit_TMB$fit$par
    TMB_FE_pars = TMB_pars[names(TMB_pars)=="beta"]
    TMB_RE_pars = TMB_pars[names(TMB_pars)=="theta"]


    num_RE_pars = length(TMB_RE_pars)
    num_RE_vars = (sqrt(1 + 8*num_RE_pars) - 1 ) / 2

    TMB_SD_pars = TMB_RE_pars[1:num_RE_vars]
    TMB_corr_pars = TMB_RE_pars[(num_RE_vars+1):length(TMB_RE_pars)]

    output = list(TMB_FE_pars = TMB_FE_pars, TMB_SD_pars = TMB_SD_pars, TMB_corr_pars = TMB_corr_pars)
    return(output)
}

#' Extract GLMM parameters from a glmmTMB object
#'
#' Uses our parameterization of SDs and correlations. Mostly serves to translate from glmmTMB parameterization.
#'
#' @inheritParams TMB_pars_list
#'
#' @returns A list containing fixed effects (\code{TMB_FEs}), standard deviations (\code{TMB_SDs}), and correlations (\code{TMB_corrs}).
#' @export
TMB_2_GLMM_pars_list <- function(fit_TMB) {
    TMB_pars = TMB_pars_list(fit_TMB)

    TMB_FEs = TMB_pars$TMB_FE_pars
    TMB_SD_pars = TMB_pars$TMB_SD_pars
    TMB_corr_pars = TMB_pars$TMB_corr_pars

    TMB_SDs = exp(TMB_SD_pars)

    ## Compute correlation matrix
    TMB_corrs = glmmTMB::get_cor(TMB_corr_pars, return_val = "vec")

    output = list(TMB_FEs = TMB_FEs, TMB_SDs = TMB_SDs, TMB_corrs = TMB_corrs)
    return(output)
}

#' Extract GLMM parameters from a glmmTMB object
#'
#' Uses our parameterization of SDs and correlations. Output is a vector organized according to our usual parameterization (see, e.g., [Sigma2theta()])
#'
#' @inheritParams TMB_pars_list
#'
#' @returns Theta, a vector of model parameters. Order is FEs for Y, RE parameters for Y, FEs for M, RE parameters for M. See, e.g., [Sigma2theta()], for details on the order of RE parameters.
#' @export
TMB_2_GLMM_pars <- function(fit_TMB) {
    pars_list = TMB_2_GLMM_pars_list(fit_TMB)

    FEs = pars_list$TMB_FEs
    SDs = pars_list$TMB_SDs
    corrs = pars_list$TMB_corrs

    TMB_theta = c(SDs, corrs)

    my_theta = TMB_theta[SD_corr2theta_indices(num_pars = length(TMB_theta))]

    output = c(FEs, my_theta)
    return(output)
}



# ---------------------------------------------------------------------------- #
#               Gradient of transformation from TMB to GLMM pars               #
# ---------------------------------------------------------------------------- #

#! Orientation: rows index outputs and columns index inputs

#' Gradient of GLMM parameters under our parameterization WRT glmmTMB parameterization.
#'
#' @inheritParams TMB_pars_list
#'
#'
#' @returns A matrix of gradients. Rows index outputs and columns index inputs.
#' @export
#'
#' @name grad_TMB_2_GLMM
#'
#' @details
#'
#'
#' `TMB_2_GLMM_grad_FE()` Returns the gradients of the fixed effects parameters (WRT all glmmTMB parameters) \
#' `TMB_2_GLMM_grad_SD()` Returns the gradients of the standard deviations (WRT all glmmTMB parameters) \
#' `TMB_2_GLMM_grad_FE()` Returns the gradients of the correlations (WRT all glmmTMB parameters)
#'
#' `TMB_2_GLMM_grad()` Combines the three sets of gradients and permutes them to match our parameterization.
#'
#'
TMB_2_GLMM_grad_FE <- function(fit_TMB) {
    TMB_pars = TMB_pars_list(fit_TMB)
    num_FE_pars = length(TMB_pars$TMB_FE_pars)
    num_RE_pars = length(TMB_pars$TMB_SD_pars) + length(TMB_pars$TMB_corr_pars)

    d_FE_FE = diag(num_FE_pars)
    d_FE_RE = matrix(0, nrow = num_FE_pars, ncol = num_RE_pars)
    d_FE = cbind(d_FE_FE, d_FE_RE)

    return(d_FE)
}


#'
#' @rdname grad_TMB_2_GLMM
#'
TMB_2_GLMM_grad_SD <- function(fit_TMB) {
    TMB_pars = TMB_pars_list(fit_TMB)
    num_FE_pars = length(TMB_pars$TMB_FE_pars)
    num_SD_pars = length(TMB_pars$TMB_SD_pars)
    num_corr_pars = length(TMB_pars$TMB_corr_pars)

    TMB_SD_pars = TMB_pars$TMB_SD_pars

    d_SD_FE = matrix(0, nrow = num_SD_pars, ncol = num_FE_pars)
    d_SD_SD = diag(exp(TMB_SD_pars))
    d_SD_corr = matrix(0, nrow = num_SD_pars, ncol = num_corr_pars)
    d_SD = cbind(d_SD_FE, d_SD_SD, d_SD_corr)

    return(d_SD)
}

#'
#' @rdname grad_TMB_2_GLMM
#'
TMB_2_GLMM_grad_corr <- function(fit_TMB) {
    TMB_pars = TMB_pars_list(fit_TMB)
    num_FE_pars = length(TMB_pars$TMB_FE_pars)
    num_SD_pars = length(TMB_pars$TMB_SD_pars)
    num_corr_pars = length(TMB_pars$TMB_corr_pars)

    TMB_corr_pars = TMB_pars$TMB_corr_pars

    d_corr_FE = matrix(0, nrow = num_corr_pars, ncol = num_FE_pars)
    d_corr_SD = matrix(0, nrow = num_corr_pars, ncol = num_SD_pars)
    d_corr_corr = t(numDeriv::jacobian(glmmTMB::get_cor, TMB_corr_pars, return_val = "vec")) # Jacobian output dimension is range-by-domain
    d_corr = cbind(d_corr_FE, d_corr_SD, d_corr_corr)

    return(d_corr)
}

#'
#' @rdname grad_TMB_2_GLMM
#'
TMB_2_GLMM_grad <- function(fit_TMB) {
    d_FE = TMB_2_GLMM_grad_FE(fit_TMB)
    d_SD = TMB_2_GLMM_grad_SD(fit_TMB)
    d_corr = TMB_2_GLMM_grad_corr(fit_TMB)

    num_FE_pars = nrow(d_FE)
    num_RE_pars = nrow(d_SD) + nrow(d_corr)

    inds_RE_unshuffle = SD_corr2theta_indices(num_pars = num_RE_pars)
    inds_all_unshuffle = c(1:num_FE_pars, num_FE_pars + inds_RE_unshuffle)

    grad_shuffled = rbind(d_FE, d_SD, d_corr)
    grad_unshuffled = grad_shuffled[inds_all_unshuffle, ]

    return(grad_unshuffled)
}



# ---------------------------------------------------------------------------- #
#                         Standard errors from glmmTMB                         #
# ---------------------------------------------------------------------------- #

#' Covariance matrix for GLMM parameters from `glmmTMB`
#'
#' @inheritParams TMB_pars_list
#'
#' @name TMB_SE
#'
#' @returns Estimated covariance matrix for the sampling distribution of GLMM parameters. `glmmTMB_SE_mat()` uses `glmmTMB`'s internal parameterization, while `TMB_2_GLMM_SE()` uses our parameterization.
#' @export
#'
glmmTMB_SE_mat <- function(fit_TMB) {
    estimates = fit_TMB$fit$par
    grad_fun = fit_TMB$obj$gr
    SE_mat = solve(numDeriv::jacobian(grad_fun, estimates))

    return(SE_mat)
}

#'
#' @rdname TMB_SE
#'
TMB_2_GLMM_SE <- function(fit_TMB) {
    TMB_SE = glmmTMB_SE_mat(fit_TMB)

    Jacob = TMB_2_GLMM_grad(fit_TMB)

    GLMM_SE = Jacob %*% TMB_SE %*% t(Jacob)

    return(GLMM_SE)
}


# Estimated covariance matrix for all parameters (i.e. from both models)
#' Estimated sampling covariance matrix for all GLMM parameters (i.e. from both models)
#'
#' @param fit_Y A GLMM fit using \code{glmmTMB} for the outcome model
#' @param fit_M A GLMM fit using \code{glmmTMB} for the mediator model
#'
#' @returns Estimated sampling covariance matrix for all parameters from both models following our parameterization.
#' @export
#'
all_pars_cov_mat_TMB <- function(fit_Y, fit_M) {
    Y_cov = TMB_2_GLMM_SE(fit_Y)
    M_cov = TMB_2_GLMM_SE(fit_M)

    Y_length = nrow(Y_cov)
    M_length = nrow(M_cov)
    joint_cov = matrix(0, nrow = M_length + Y_length, ncol = M_length + Y_length)
    joint_cov[1:Y_length, 1:Y_length] = Y_cov
    joint_cov[(Y_length + 1):(Y_length + M_length), (Y_length + 1):(Y_length + M_length)] = M_cov

    # Make output matrix symmetric
    joint_cov = (joint_cov + t(joint_cov))/2
    return(joint_cov)
}


# fit_Y = fit_Y_TMB
# fit_M = fit_M_TMB
