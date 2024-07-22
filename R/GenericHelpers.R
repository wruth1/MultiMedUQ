
#' Title
#'
#' @param x A numeric vector
#' @name logit_expit
#'
#' @return A numeric vector of the same length as x
#' @export
#'
#' @examples
#' logit(c(0.1, 0.2, 0.3))
#' expit(c(0.1, 0.2, 0.3))
#'
#' logit(expit(c(0.1, 0.2, 0.3))
#' expit(logit(c(0.1, 0.2, 0.3))
logit = function(x){
  log(x/(1-x))
}

#' @rdname logit_expit
#' @export
expit = function(x){
  exp(x)/(1+exp(x))
}
