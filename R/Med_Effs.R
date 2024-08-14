
#' Compute mediation effects on different scales
#'
#' @param ENC1,ENC2 Two expected nested counterfactuals
#' @name ME_scales
#'
#' @return The mediation effect on the given scale.
#' @export
ME_diff = function(ENC1, ENC2){
  return(ENC1 - ENC2)
}

#' @rdname ME_scales
#' @export
ME_rat = function(ENC1, ENC2){
  return(ENC1 / ENC2)
}

#' @rdname ME_scales
#' @export
ME_OR = function(ENC1, ENC2){
  return(ENC1 / (1 - ENC1) / (ENC2 / (1 - ENC2)))
}


#' Compute a mediation effect on the given scale(s)
#'
#' @param ENC1,ENC2 Two expected nested counterfactuals
#' @param scale The scale(s) of the mediation effect. Can be "diff", "rat" or "OR".
#'
#' @return The mediation effect on the given scale.
#' @export
get_ME = function(ENC1, ENC2, scale = c("diff", "rat", "OR")){
  output = c()
  if("diff" %in% scale){
    output = c(output, ME_diff(ENC1, ENC2))
  }
  if("rat" %in% scale){
    output = c(output, ME_rat(ENC1, ENC2))
  }
  if("OR" %in% scale){
    output = c(output, ME_OR(ENC1, ENC2))
  }

  if(length(output) == 0){
    stop("Unknown scale")
  } else{
    return(output)
  }
}



# Different flavours of mediation effect

#' Mediation effects of \eqn{X} on \eqn{Y} mediated by \eqn{M}
#'
#' @param scale The scale(s) of the mediation effect. Can be "diff", "rat" or "OR".
#' @param w Level of covariates, \eqn{W}.
#' @param b_Y,b_M Coefficient vectors for \eqn{Y}-model and \eqn{M}-model, respectively.
#' @param theta_Y,theta_M Covariance parameters of random effects in \eqn{Y}-model and \eqn{M}-model, respectively. See details.
#'
#' @name Med_Effs
#'
#' @return The specified mediation effect of \eqn{X} on \eqn{Y} mediated by \eqn{M}.
#' @export
#'
#' @details
#' Contents of \code{b_Y} are \code{(b_Y_0, b_Y_X, b_Y_M, B_Y_W)}. Contents of \code{b_M} are \code{(b_M_0, b_M_X, B_M_W)}.
#'
#' Contents of \code{theta_Y} are \code{(s_Y_0, cor_Y_0X, cor_Y_0M, s_Y_X, cor_Y_XM, s_Y_M)}. Contents of \code{theta_M} are \code{(s_M_0, cor_M_0X, s_M_X)}.
#'
total_effect <- function(scale = c("diff", "rat", "OR"), w, b_Y, theta_Y, b_M, theta_M){

  ENC1 = ENC(1, 1, w, b_Y, theta_Y, b_M, theta_M)
  ENC2 = ENC(0, 0, w, b_Y, theta_Y, b_M, theta_M)

  return(get_ME(ENC1, ENC2, scale))
}


#' @rdname Med_Effs
#'
#' @param x_m_ref Reference value of \eqn{X} for determining the mediator level.
#'
#' @export
direct_effect <- function(scale = c("diff", "rat", "OR"), w, b_Y, theta_Y, b_M, theta_M, x_m_ref = 0){

  ENC1 = ENC(1, x_m_ref, w, b_Y, theta_Y, b_M, theta_M)
  ENC2 = ENC(0, x_m_ref, w, b_Y, theta_Y, b_M, theta_M)

  return(get_ME(ENC1, ENC2, scale))
}


#' @rdname Med_Effs
#'
#' @param x_ref Reference value of \eqn{X}.
#'
#' @export
indirect_effect <- function(scale = c("diff", "rat", "OR"), w, b_Y, theta_Y, b_M, theta_M, x_ref = 1){

  ENC1 = ENC(x_ref, 1, w, b_Y, theta_Y, b_M, theta_M)
  ENC2 = ENC(x_ref, 0, w, b_Y, theta_Y, b_M, theta_M)

  return(get_ME(ENC1, ENC2, scale))
}



#' Total, direct and indirect mediation effects on the specified scale(s)
#'
#' @param scale The scale(s) of the mediation effect. Can be "diff", "rat" or "OR".
#' @param w Level of covariates, \eqn{W}.
#' @param b_Y,b_M Coefficient vectors for \eqn{Y}-model and \eqn{M}-model, respectively.
#' @param theta_Y,theta_M Covariance parameters of random effects in \eqn{Y}-model and \eqn{M}-model, respectively. See details.
#' @param x_ref,x_m_ref Reference values of \eqn{X}, respectively for \eqn{X} and for determining the value of \eqn{M}.
#'
#' @return All mediation effects (total, direct and indirect), on the specified scale(s). Order is total, direct, indirect. Within each effect, order is as specified in \code{scale}. Default is difference, ratio, odds-ratio.
#' @export
#'
all_MEs <- function(scale = c("diff", "rat", "OR"), w, b_Y, theta_Y, b_M, theta_M, x_ref = 1, x_m_ref = 0){
  MEs = c(total_effect(scale, w, b_Y, theta_Y, b_M, theta_M),
    direct_effect(scale, w, b_Y, theta_Y, b_M, theta_M, x_m_ref),
    indirect_effect(scale, w, b_Y, theta_Y, b_M, theta_M, x_ref))

  ME_names = as.vector(t(outer(c("total", "direct", "indirect"), scale, paste, sep = "_")))
  names(MEs) = ME_names

  return(MEs)
}





# Gradient of mediation effects as functions of ENCs (expected nested counterfactuals)
## Recall that the order of \eqn{(X, X_M)} levels in ENC is (1,1), (1,0), (0,1), (0,0).

#' Title
#'
#' @param ENC_11,ENC_10,ENC_01,ENC_00 Expected nested counterfactuals at specified levels of x and x_m respectively.
#'
#' @name grad_med_effs
#'
#' @return Gradient of the mediation effect of specified flavour and scale.
#' @export
grad_TE_diff <- function(ENC_11, ENC_10, ENC_01, ENC_00){
  return(c(1, 0, 0, -1))
}

#' @rdname grad_med_effs
#' @export
grad_TE_rat <- function(ENC_11, ENC_10, ENC_01, ENC_00){
  return(c(1/ENC_00, 0, 0, -ENC_11/ENC_00^2))
}

#' @rdname grad_med_effs
#' @export
grad_TE_or <- function(ENC_11, ENC_10, ENC_01, ENC_00){
  d1 = (1 / (1 - ENC_11)^2) / (ENC_00 / (1 - ENC_00))
  d2 = 0
  d3 = 0
  d4 = (ENC_11 / (1 - ENC_11)) / (ENC_00^2)

  return(c(d1, d2, d3, d4))
}


#' @rdname grad_med_effs
#' @export
grad_DE_diff <- function(ENC_11, ENC_10, ENC_01, ENC_00){
  return(c(0, 1, 0, -1))
}

#' @rdname grad_med_effs
#' @export
grad_DE_rat <- function(ENC_11, ENC_10, ENC_01, ENC_00){
  return(c(0, 1/ENC_00, 0, -ENC_10/ENC_00^2))
}

#' @rdname grad_med_effs
#' @export
grad_DE_or <- function(ENC_11, ENC_10, ENC_01, ENC_00){
  d1 = 0
  d2 = (1 / (1 - ENC_10)^2) / (ENC_00 / (1 - ENC_00))
  d3 = 0
  d4 = (ENC_10 / (1 - ENC_10)) / (ENC_00^2)

  return(c(d1, d2, d3, d4))
}

#' @rdname grad_med_effs
#' @export
grad_IE_diff <- function(ENC_11, ENC_10, ENC_01, ENC_00){
  return(c(1, 0, -1, 0))
}

#' @rdname grad_med_effs
#' @export
grad_IE_rat <- function(ENC_11, ENC_10, ENC_01, ENC_00){
  return(c(1/ENC_00, 0, -ENC_11/ENC_01^2, 0))
}

#' @rdname grad_med_effs
#' @export
grad_IE_or <- function(ENC_11, ENC_10, ENC_01, ENC_00){
  d1 = (1 / (1 - ENC_11)^2) / (ENC_01 / (1 - ENC_01))
  d2 = 0
  d3 = (ENC_11 / (1 - ENC_11)) / (ENC_01^2)
  d4 = 0

  return(c(d1, d2, d3, d4))
}



#' Jacobian of all mediation effects wrt all ENCs
#'
#' @param scale The scale(s) of the mediation effect. Can be "diff", "rat" or "OR".
#' @param w Level of covariates, \eqn{W}.
#' @param b_Y,b_M Coefficient vectors for \eqn{Y}-model and \eqn{M}-model, respectively.
#' @param theta_Y,theta_M Covariance parameters of random effects in \eqn{Y}-model and \eqn{M}-model, respectively. See details.
#' @name all_grad_MEs
#'
#' @return A matrix of partial derivatives. Dimension is (3 * length(scale))-by-4.
#' @export
#'
all_grad_MEs_pars <- function(scale = c("diff", "rat", "OR"), w, b_Y, theta_Y, b_M, theta_M){

  ENCs = all_ENCs(w, b_Y, theta_Y, b_M, theta_M)

  all_TE_grads = c()
  all_DE_grads = c()
  all_IE_grads = c()

  if("diff" %in% scale){
    all_TE_grads = rbind(all_TE_grads, grad_TE_diff(ENCs[1], ENCs[2], ENCs[3], ENCs[4]))
    all_DE_grads = rbind(all_DE_grads, grad_DE_diff(ENCs[1], ENCs[2], ENCs[3], ENCs[4]))
    all_IE_grads = rbind(all_IE_grads, grad_IE_diff(ENCs[1], ENCs[2], ENCs[3], ENCs[4]))
  }
  if("rat" %in% scale){
    all_TE_grads = rbind(all_TE_grads, grad_TE_rat(ENCs[1], ENCs[2], ENCs[3], ENCs[4]))
    all_DE_grads = rbind(all_DE_grads, grad_DE_rat(ENCs[1], ENCs[2], ENCs[3], ENCs[4]))
    all_IE_grads = rbind(all_IE_grads, grad_IE_rat(ENCs[1], ENCs[2], ENCs[3], ENCs[4]))
  }
  if("OR" %in% scale){
    all_TE_grads = rbind(all_TE_grads, grad_TE_or(ENCs[1], ENCs[2], ENCs[3], ENCs[4]))
    all_DE_grads = rbind(all_DE_grads, grad_DE_or(ENCs[1], ENCs[2], ENCs[3], ENCs[4]))
    all_IE_grads = rbind(all_IE_grads, grad_IE_or(ENCs[1], ENCs[2], ENCs[3], ENCs[4]))
  }

  return(rbind(all_TE_grads, all_DE_grads, all_IE_grads))

}

#' @rdname all_grad_MEs
#' @export
all_grad_MEs_models <- function(scale, w, fit_Y, fit_M){
  info_Y = get_model_pars(fit_Y)
  info_M = get_model_pars(fit_M)

  b_Y = info_Y[["b"]]
  theta_Y = info_Y[["theta"]]
  b_M = info_M[["b"]]
  theta_M = info_M[["theta"]]

  return(all_grad_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M))
}



#' Covariance matrix of all mediation effects on various scales.
#'
#' @param scale The scale(s) of the mediation effect. Can be "diff", "rat" or "OR".
#' @param w Level of covariates, \eqn{W}.
#' @param fit_Y,fit_M Fitted models for Y and M.
#'
#' @return A covariance matrix for all mediation effects (total, direct and indirect) on the specified scale(s).
#' @export
all_cov_MEs <- function(scale = c("diff", "rat", "OR"), w, fit_Y, fit_M){
  cov_ENCs = all_covs_ENC(w, fit_Y, fit_M)

  grad_MEs = all_grad_MEs_models(scale, w, fit_Y, fit_M)

  cov_MEs = grad_MEs %*% cov_ENCs %*% t(grad_MEs)

  return(cov_MEs)
}
