# Intermediate integral for computing mediation effects ####

#' Intermediate integral for computing mediation effects
#'
#' This function computes the integral of the function \code{dnorm(x)/(1 + exp(-mu - sigma*x))} from -Inf to Inf. I.e. The expected inverse-logit of a RV with mean \code{mu} and SD \code{sigma}.
#' This is useful for computing various mediation effects.
#'
#' @param mu Mean of the latent Gaussian.
#' @param sigma SD of the latent Gaussian.
#'
#' @export
#'
#' @return The integral of the function \code{dnorm(x)/(1 + exp(-mu - sigma*x))} from -Inf to Inf.
#'
psi <- function(mu, sigma){
  integrand <- function(x){
    A = 1 + exp(-mu - sigma*x)
    B = stats::dnorm(x)

    return(B/A)
  }

  integral = stats::integrate(integrand, -Inf, Inf)
  # print(integral$message)
  # print(integral$abs.error)

  return(integral$value)
}


## Gradient of psi ####

#' Gradient of the intermediate integral function, psi
#'
#' @param mu Mean of the latent Gaussian.
#' @param sigma SD of the latent Gaussian.
#'
#' @name grad_psi
#'
#' @export
#'
#' @return The gradient of the function \code{psi(mu, sigma)} with respect to \code{mu} and \code{sigma}.
#'
d1_psi = function(mu, sigma){
  integrand = function(x){
    A = stats::dnorm(x)

    r = mu + sigma*x
    B = exp(-r/2) + exp(r/2)

    # return(C*B/A^2)
    return(A / (B^2))
  }

  return(stats::integrate(integrand, -Inf, Inf)$value)
}


#' @rdname grad_psi
#' @export
d2_psi = function(mu, sigma){
  integrand = function(x){
    A = stats::dnorm(x)

    r = mu + sigma*x
    B = exp(-r/2) + exp(r/2)

    return(x * A / B^2)
  }

  return(stats::integrate(integrand, -Inf, Inf)$value)
}

#'
#' #' Gradient of psi, chain rule version
#' #'
#' #' This function computes the gradient of the function \code{psi(mu, sigma)}, where \code{mu} and \code{sigma} are themselves functions of the variable with respect to which we are differentiating.
#' #'
#' #' @param mu,sigma Arguments of psi
#' #' @param grad_mu,grad_sigma Gradients of the arguments, \code{mu} and \code{sigma}.
#' #'
#' #' @return Total gradient of psi.
#' #' @export
#' grad_psi_chain <- function(mu, sigma, grad_mu, grad_sigma){
#'   d1_psi(mu, sigma) * grad_mu + d2_psi(mu, sigma) * grad_sigma
#' }
#'
