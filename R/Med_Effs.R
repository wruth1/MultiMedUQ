
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
#' @return All mediation effects (total, direct and indirect), on the specified scale(s)
#' @export
#'
all_MEs <- function(scale = c("diff", "rat", "OR"), w, b_Y, theta_Y, b_M, theta_M, x_ref = 1, x_m_ref = 0){
  c(total_effect(scale, w, b_Y, theta_Y, b_M, theta_M),
    direct_effect(scale, w, b_Y, theta_Y, b_M, theta_M, x_m_ref),
    indirect_effect(scale, w, b_Y, theta_Y, b_M, theta_M, x_ref))
}





