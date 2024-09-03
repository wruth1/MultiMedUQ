# General ####

logit = function(x) log(x/(1-x))

expit = function(x) exp(x)/(1+exp(x))



# Cont Resp, Bin Med, Fixed-Effs ####
## Some notation:
  ## eta: a_0 + a_1 x + A_2 W = a_0 + a_x * x + A_2' w. Linear predictor for M
  ## zeta: b_0 + b_1 * m + b_2 * x + B_3' w = b_0 + b_m * m + b_x * x + B_3' w. Linear predictor for Y

  ## delta: E(M | X=x+1) - E(M | X=x)
  ## gamma: E(Y | X=x+1) - E(Y | X=x). Total effect of X on Y


get_delta = function(eta, a_x){
  Q1 = exp(-eta)
  Q2 = (1 - exp(-a_x))
  Q3 = (1 + exp(-eta - a_x))
  Q4 = (1 + exp(-eta))

  return(Q1 * Q2 / (Q3 * Q4))
}


get_gamma = function(delta, b_m, b_x){
  return(delta * b_m + b_x)
}


## Gradient of total effect wrt the a's and b's

### First, some intermediate quantities
d_delta_d_eta = function(eta, a_x){
  delta = get_delta(eta, a_x)

  Q1 = 1 - exp(-2*eta - a_x)
  Q2 = 1 + exp(-eta)
  Q3 = 1 + exp(-eta - a_x)

  return( (-1) * delta * Q1 / (Q2 * Q3))
}

d_delta_d_a_x = function(eta, a_x){
  Q1 = exp(-2*eta - a_x)
  Q2 = exp(-eta - a_x)
  Q3 = 1 + exp(-eta)
  Q4 = 1 + exp(-eta - a_x)

  return((Q1 + Q2)/(Q3 * Q4^2))
}


### Now, the gradient
d_gamma_d_a_0 <- function(dd_de, b_m) {
  return(dd_de * b_m)
}

d_gamma_d_a_1 <- function(dd_de, dd_da_x, x, b_m) {
  return(b_m * x * dd_de + b_m * dd_da_x)
}

d_gamma_d_A_2 <- function(dd_de, W, b_m) {
  return(c(dd_de * b_m) * W)
}




d_gamma_d_b_0 <- function() {
  return(0)
}

d_gamma_d_b_1 <- function(eta, a_x) {
  return(get_delta(eta, a_x))
}

d_gamma_d_b_2 <- function() {
  return(1)
}

d_gamma_d_B_3 <- function(W) {
  return(rep(0, times = length(W)))
}

d_gamma_d_theta <- function(eta, x, W, a, b) {
  a_0 <- a[1]
  a_x <- a[2]
  A_2 <- a[3:length(a)]

  b_0 <- b[1]
  b_m <- b[2]
  b_x <- b[3]
  B_3 <- b[4:length(b)]

  dd_de <- d_delta_d_eta(eta, a_x)
  dd_da_x <- d_delta_d_a_x(eta, a_x)

  dg_da_0 <- d_gamma_d_a_0(dd_de, b_m)
  dg_da_1 <- d_gamma_d_a_1(dd_de, dd_da_x, x, b_m)
  dg_dA_2 <- d_gamma_d_A_2(dd_de, W, b_m)

  dg_db_0 <- d_gamma_d_b_0()
  dg_db_1 <- d_gamma_d_b_1(eta, a_x)
  dg_db_2 <- d_gamma_d_b_2()
  dg_dB_3 <- d_gamma_d_B_3(W)

  return(c(dg_da_0, dg_da_1, dg_dA_2, dg_db_0, dg_db_1, dg_db_2, dg_dB_3))
}







######################################################
#### Binary Response, Binary Mediator, Fixed-Effs ####
######################################################


# Compute mediation effect

## Probability of Y=1 given inputs
PY1 <- function(eta, zeta, beta) {
  Q1 <- 1 + exp(-zeta)
  Q2 <- 1 + exp(-eta)
  Q3 <- 1 + exp(-zeta - beta)

  num <- Q1 + (Q2 - 1) * Q3
  den <- Q1 * Q2 * Q3

  return(num / den)
}

## Probability of Y=0 given inputs
PY0 <- function(eta, zeta, beta) {
  Q1 <- exp(-zeta)
  Q2 <- exp(-eta)
  Q3 <- exp(-beta)

  num <- Q1 * (Q2 + Q3 + Q1 * Q3 + Q1 * Q2 * Q3)
  den <- (1 + Q1) * (1 + Q2) * (1 + Q1 * Q3)

  return(num / den)
}


get_odds_ref <- function(eta, zeta, beta) {
  num <- PY1(eta, zeta, beta)
  den <- PY0(eta, zeta, beta)

  return(num / den)
}

get_odds <- function(eta, zeta, beta) {
  Q1 <- exp(-zeta)
  Q1_inv <- 1 / Q1
  Q2 <- exp(-eta)
  Q3 <- exp(-beta)

  num <- 1 + Q1_inv + Q2 * (Q3 + Q1_inv)
  den <- Q2 + Q3 + Q1 * Q3 + Q1 * Q2 * Q3

  return(num / den)
}

get_odds_ratio <- function(eta, a_x, zeta, b_x, b_m) {
  odds_1 <- get_odds(eta + a_x, zeta + b_x, b_m)
  odds_2 <- get_odds(eta, zeta, b_m)
  OR_hat <- odds_1 / odds_2
  return(OR_hat)
}






# Gradients of OR, obtained from Maple
# I verified all of these against a simple finite difference calculation

d_OR_d_eta <- function(eta, a_x, zeta, b_x, b_m) {
  -exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m)) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) - (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) ^ 2 / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-eta - a_x) - exp(-zeta - b_x - eta - a_x - b_m)) + (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) ^ 2 * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * exp(-eta) * (1 + exp(-zeta - b_m)) + (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (-exp(-eta) - exp(-zeta - eta - b_m))
}

d_OR_d_a_x <- function(eta, a_x, zeta, b_x, b_m) {
  -exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m)) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) - (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) ^ 2 / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-eta - a_x) - exp(-zeta - b_x - eta - a_x - b_m))
}

d_OR_d_zeta <- function(eta, a_x, zeta, b_x, b_m) {
  (-exp(-zeta - b_x) - exp(-eta - a_x) * exp(-zeta - b_x - b_m)) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) - (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) ^ 2 / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-zeta - b_x - b_m) - exp(-zeta - b_x - eta - a_x - b_m)) - (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) ^ 2 * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-zeta) - exp(-eta) * exp(-zeta - b_m)) + (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (-exp(-zeta - b_m) - exp(-zeta - eta - b_m))
}



d_OR_d_b_x <- function(eta, a_x, zeta, b_x, b_m) {
  (-exp(-zeta - b_x) - exp(-eta - a_x) * exp(-zeta - b_x - b_m)) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) +
    (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) -
    (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) ^ 2 / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-zeta - b_x - b_m) - exp(-zeta - b_x - eta - a_x - b_m))
}

d_OR_d_b_m <- function(eta, a_x, zeta, b_x, b_m) {
  -exp(-eta - a_x) * exp(-zeta - b_x - b_m) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) -
    (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) ^ 2 / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * (-exp(-b_m) - exp(-zeta - b_x - b_m) - exp(-zeta - b_x - eta - a_x - b_m)) +
    (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) ^ 2 * exp(-zeta) * (exp(-eta) + exp(-b_m) + exp(-zeta - b_m) + exp(-zeta - eta - b_m)) * exp(-eta) * exp(-zeta - b_m) +
    (1 + exp(-zeta - b_x) + exp(-eta - a_x) * (1 + exp(-zeta - b_x - b_m))) / exp(-zeta - b_x) / (exp(-eta - a_x) + exp(-b_m) + exp(-zeta - b_x - b_m) + exp(-zeta - b_x - eta - a_x - b_m)) / (1 + exp(-zeta) + exp(-eta) * (1 + exp(-zeta - b_m))) * exp(-zeta) * (-exp(-b_m) - exp(-zeta - b_m) - exp(-zeta - eta - b_m))
}





# Build gradient wrt reg pars
d_OR_d_a0 <- function(eta, a_x, zeta, b_x, b_m) {
  return(d_OR_d_eta(eta, a_x, zeta, b_x, b_m))
}

d_OR_d_a1 <- function(eta, a_x, zeta, b_x, b_m, x_ref) {
  return(x_ref * d_OR_d_eta(eta, a_x, zeta, b_x, b_m) + d_OR_d_a_x(eta, a_x, zeta, b_x, b_m))
}

d_OR_d_A2 <- function(eta, a_x, zeta, b_x, b_m, W_ref) {
  return(W_ref * c(d_OR_d_eta(eta, a_x, zeta, b_x, b_m)))
}

d_OR_d_b0 <- function(eta, a_x, zeta, b_x, b_m) {
  return(d_OR_d_zeta(eta, a_x, zeta, b_x, b_m))
}

d_OR_d_b1 <- function(eta, a_x, zeta, b_x, b_m) {
  return(d_OR_d_b_m(eta, a_x, zeta, b_x, b_m))
}

d_OR_d_b2 <- function(eta, a_x, zeta, b_x, b_m, x_ref) {
  return(x_ref * d_OR_d_zeta(eta, a_x, zeta, b_x, b_m) + d_OR_d_b_x(eta, a_x, zeta, b_x, b_m))
}

d_OR_d_B3 <- function(eta, a_x, zeta, b_x, b_m, W_ref) {
  return(W_ref * c(d_OR_d_zeta(eta, a_x, zeta, b_x, b_m)))
}

d_OR_d_theta <- function(eta, a_x, zeta, b_x, b_m, x_ref, W_ref) {
  return(c(d_OR_d_a0(eta, a_x, zeta, b_x, b_m),
            d_OR_d_a1(eta, a_x, zeta, b_x, b_m, x_ref),
            d_OR_d_A2(eta, a_x, zeta, b_x, b_m, W_ref),
            d_OR_d_b0(eta, a_x, zeta, b_x, b_m),
            d_OR_d_b1(eta, a_x, zeta, b_x, b_m),
            d_OR_d_b2(eta, a_x, zeta, b_x, b_m, x_ref),
            d_OR_d_B3(eta, a_x, zeta, b_x, b_m, W_ref)))
}









######################################################
#### Binary Response, Binary Mediator, Mixed-Effs ####
######################################################


phi = function(mu, sigma){
  integrand = function(x){
    A = 1 + exp(-mu - sigma*x)
    B = dnorm(x)

    return(B/A)
  }

  integral = integrate(integrand, -Inf, Inf)
  # print(integral$message)
  # print(integral$abs.error)

  return(integral$value)
}

# Partial derivatives of phi
# I validated these against a simple finite difference approximation

## d phi / d mu
d1_phi = function(mu, sigma){
  integrand = function(x){
    A = dnorm(x)

    r = mu + sigma*x
    B = exp(-r/2) + exp(r/2)

    # return(C*B/A^2)
    return(A / (B^2))
  }

  return(integrate(integrand, -Inf, Inf)$value)
}


## d phi / d sigma
d2_phi = function(mu, sigma){
  integrand = function(x){
    A = dnorm(x)

    r = mu + sigma*x
    B = exp(-r/2) + exp(r/2)

    return(x * A / B^2)
  }

  return(integrate(integrand, -Inf, Inf)$value)
}




  a1 = b_M[2]
  b1 = b_Y[3]
  b2 = b_Y[2]

  s1 = sigma_fun(1, theta_Y[1], theta_Y[3], theta_Y[2])
  s2 = sigma_fun(1, theta_M[1], theta_M[3], theta_M[2])

  x = 1
  x_m = 1
  w = c(0,0)

  zeta = b_Y[1]  + b_Y[4] * w[1] + b_Y[5] * w[2]  
  eta = b_M[1]  + b_M[3] * w[1] + b_M[4] * w[2]
  
  Y1_X1_M1 = phi(zeta + b1 + b2, s1)  # P(Y=1 | X=1, M=1)
  M1_X1 = phi(eta + a1, s2)           # P(M=1 | X=1)
  Y1_X1_M0 = phi(zeta + b2, s1)       # P(Y=1 | X=1, M=0)
  M0_X1 = (1 - phi(eta + a1, s2))     # P(M=0 | X=1)

  # These should match, but don't    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Y1_X1_M1 * M1_X1 + Y1_X1_M0 * M0_X1
  ENC(1, 1, w, b_Y, theta_Y, b_M, theta_M, which_REs)
  ENC(1, 0, w, b_Y, theta_Y, b_M, theta_M, which_REs)




# Total mediation effect (on odds-ratio scale)
## In my notation: Phi(eta, zeta, a_x, b_m, b_x, sigma_U(x+1), sigma_V(x+1), sigma_U(x), sigma_V(x))
## Note: eta, a and U come from the M-model, while zeta, b and V come from the Y-model
## eta and zeta are linear predictors
##
## Phi = (P(Y=1 | X=1) / P(Y=0 | X=1)) / (P(Y=1 | X=0) / P(Y=0 | X=0))
Phi = function(eta, zeta, a1, b1, b2, s1, s2, s3, s4){
  

  p_Y1_X1 = phi(zeta + b1 + b2, s1) * phi(eta + a1, s2) + phi(zeta + b2, s1) * (1 - phi(eta + a1, s2))
  p_Y0_X1 = 1 - phi(zeta + b2, s1) + phi(eta + a1, s2) * (phi(zeta + b2, s1) - phi(zeta + b1 + b2, s1))
  p_Y1_X0 = (phi(zeta + b1, s3) * phi(eta, s4) + phi(zeta, s3) * (1 - phi(eta, s4)))
  p_Y0_X0 = 1 - phi(zeta, s3) + phi(eta, s4) * (phi(zeta, s3) - phi(zeta + b1, s3))

  num = p_Y1_X1 / p_Y0_X1
  denom = p_Y1_X0 / p_Y0_X0
  return(num/denom)
}

# [1 - phi(zeta + b2, s1) + phi(eta + a1, s2) * (phi(zeta + b2, s1) - phi(zeta + b1 + b2, s2))] * [(phi(zeta + b1, s3) * phi(eta, s3) + phi(zeta, s3) * (1 - phi(eta, s4)))]
  
  
#   (phi(zeta + b1 + b2, s1) * phi(eta + a1, s2) + phi(zeta + b2, s1) * (1 - phi(eta + a1, s2))) * (1 - phi(zeta, s3) + phi(eta, s4) * (phi(zeta, s3) - phi(zeta + b1, s3))) / (1 - phi(zeta + b2, s1) + phi(eta + a1, s2) * (phi(zeta + b2, s1) - phi(zeta + b1 + b2, s2))) / (phi(zeta + b1, s3) * phi(eta, s3) + phi(zeta, s3) * (1 - phi(eta, s4)))



# args2
#
# e = 0.000001
# (Phi(args2[1] + e, args2[2], args2[3], args2[4], args2[5], args2[6], args2[7], args2[8], args2[9]) - do.call(Phi, as.list(args2))) / e
#
# do.call(grad_Phi, as.list(args2))


# Derivatives of Phi (the total mediation effect)
# All formulas have been validated against a simple finite difference approximation

## Most of the formula comes directly from Maple. However, there is one piece of Maple syntax that I wasn't able to fully accommodate: expressing some phi derivatives as, e.g., diff(phi(a, s4), a). Specifically, Maple uses different syntax for the a-derivative of phi(a, g(b)) vs the a-derivative of phi(f(a), g(b)). I could fix the latter within Maple, but couldn't find a good solution for the former.
## To solve this, I did some find/replace in VSCode using regular expressions with groups. One such pattern was: diff\(phi\(a, (\w+)\), a\) -> d1_phi(a, $1), which creates a "group" consisting of whatever is matched by the \w+, and inserts this group into the replacement string as indicated by the $1.







#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Old formula for Phi was wrong (a couple typos). It's been fixed, but the partial derivatives have not been updated.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





## d Phi / d a
dPhi_da = function(a, b, a1, b1, b2, s1, s2, s3, s4){
  q1 = (phi(b + b1 + b2, s1) * d1_phi(a + a1, s2) - phi(b + b2, s1) * d1_phi(a + a1, s2))
  q2 = (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2)))
  q3 = (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3)))
  q4 = (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4)))

  remainder = (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) ^ 2 * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) * d1_phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2)) + (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) * d1_phi(a, s4) * (phi(b, s3) - phi(b + b1, s3)) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) - (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) ^ 2 * (phi(b + b1, s3) * d1_phi(a, s3) - phi(b, s3) * d1_phi(a, s4))


  q1 / q2 * q3 / q4 - remainder
}


## d Phi / d b
dPhi_db = function(a, b, a1, b1, b2, s1, s2, s3, s4){
  (d1_phi(b + b1 + b2, s1) * phi(a + a1, s2) + d1_phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) + (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (-d1_phi(b, s3) + phi(a, s4) * (d1_phi(b, s3) - d1_phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) - (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) ^ 2 / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) * (-d1_phi(b + b2, s1) + phi(a + a1, s2) * (d1_phi(b + b2, s1) - d1_phi(b + b1 + b2, s2))) - (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) ^ 2 * (d1_phi(b + b1, s3) * phi(a, s3) + d1_phi(b, s3) * (1 - phi(a, s4)))
}


## d Phi / d a1
dPhi_da1 = function(a, b, a1, b1, b2, s1, s2, s3, s4){
  (phi(b + b1 + b2, s1) * d1_phi(a + a1, s2) - phi(b + b2, s1) * d1_phi(a + a1, s2)) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) - (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) ^ 2 / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) * d1_phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))
}


## d Phi / d b1
dPhi_db1 = function(a, b, a1, b1, b2, s1, s2, s3, s4){
  d1_phi(b + b1 + b2, s1) * phi(a + a1, s2) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) - (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * phi(a, s4) * d1_phi(b + b1, s3) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) + (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) ^ 2 / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) * phi(a + a1, s2) * d1_phi(b + b1 + b2, s2) - (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) ^ 2 * d1_phi(b + b1, s3) * phi(a, s3)
}


## d Phi / d b2
dPhi_db2 = function(a, b, a1, b1, b2, s1, s2, s3, s4){
  (d1_phi(b + b1 + b2, s1) * phi(a + a1, s2) + d1_phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) - (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) ^ 2 / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) * (-d1_phi(b + b2, s1) + phi(a + a1, s2) * (d1_phi(b + b2, s1) - d1_phi(b + b1 + b2, s2)))
}


## d Phi / d s1
dPhi_ds1 = function(a, b, a1, b1, b2, s1, s2, s3, s4){
  (d2_phi(b + b1 + b2, s1) * phi(a + a1, s2) + d2_phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) - (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) ^ 2 / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) * (-d2_phi(b + b2, s1) + phi(a + a1, s2) * d2_phi(b + b2, s1))
}


## d Phi / d s2
dPhi_ds2 = function(a, b, a1, b1, b2, s1, s2, s3, s4){
  (phi(b + b1 + b2, s1) * d2_phi(a + a1, s2) - phi(b + b2, s1) * d2_phi(a + a1, s2)) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) - (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) ^ 2 / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) * (d2_phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2)) - phi(a + a1, s2) * d2_phi(b + b1 + b2, s2))
}


## d Phi / d s3
dPhi_ds3 = function(a, b, a1, b1, b2, s1, s2, s3, s4){
  (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (-d2_phi(b, s3) + phi(a, s4) * (d2_phi(b, s3) - d2_phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) - (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) ^ 2 * (d2_phi(b + b1, s3) * phi(a, s3) + phi(b + b1, s3) * d2_phi(a, s3) + d2_phi(b, s3) * (1 - phi(a, s4)))
}


## d Phi / d s4
dPhi_ds4 = function(a, b, a1, b1, b2, s1, s2, s3, s4){
  (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * d2_phi(a, s4) * (phi(b, s3) - phi(b + b1, s3)) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) + (phi(b + b1 + b2, s1) * phi(a + a1, s2) + phi(b + b2, s1) * (1 - phi(a + a1, s2))) * (1 - phi(b, s3) + phi(a, s4) * (phi(b, s3) - phi(b + b1, s3))) / (1 - phi(b + b2, s1) + phi(a + a1, s2) * (phi(b + b2, s1) - phi(b + b1 + b2, s2))) / (phi(b + b1, s3) * phi(a, s3) + phi(b, s3) * (1 - phi(a, s4))) ^ 2 * phi(b, s3) * d2_phi(a, s4)
}




# Full gradient of Phi
grad_Phi = function(a, b, a1, b1, b2, s1, s2, s3, s4){
  return(c(dPhi_da(a, b, a1, b1, b2, s1, s2, s3, s4),
            dPhi_db(a, b, a1, b1, b2, s1, s2, s3, s4),
            dPhi_da1(a, b, a1, b1, b2, s1, s2, s3, s4),
            dPhi_db1(a, b, a1, b1, b2, s1, s2, s3, s4),
            dPhi_db2(a, b, a1, b1, b2, s1, s2, s3, s4),
            dPhi_ds1(a, b, a1, b1, b2, s1, s2, s3, s4),
            dPhi_ds2(a, b, a1, b1, b2, s1, s2, s3, s4),
            dPhi_ds3(a, b, a1, b1, b2, s1, s2, s3, s4),
            dPhi_ds4(a, b, a1, b1, b2, s1, s2, s3, s4)))
}


# # Finite difference verification
# a = 1
# b = 1
# a1 = 1
# b1 = 1
# b2 = 1
# s1 = 0.5
# s2 = 0.5
# s3 = 0.5
# s4 = 0.5
#
# e = 0.000001
#
# (Phi(a, b, a1, b1, b2, s1, s2, s3, s4 + e) - Phi(a, b, a1, b1, b2, s1, s2, s3, s4)) / e
# dPhi_ds4(a, b, a1, b1, b2, s1, s2, s3, s4)
#
#
# grad_Phi(a, b, a1, b1, b2, s1, s2, s3, s4)



# Transform from GLMM parameters to arguments of Phi
## Note: This transformation also depends on x and W, although we will not worry about the gradient wrt these.

xi = function(a, theta, b, gamma, x, W){

  # Extract fixed effects
  a_0 = a[1]
  a_x = a[2]
  A_2 = a[3:5]

  b_0 = b[1]
  b_m = b[2]
  b_x = b[3]
  B_3 = b[4:6]


  # Compute linear predictors
  eta = as.numeric(a_0 + a_x * x + W %*% A_2)
  zeta = as.numeric(b_0 + b_x * x + W %*% B_3)


  # Random effects covariances
  s_M_0 = theta[1]
  s_M_x = theta[2]
  rho_M = theta[3]

  s_Y_0 = gamma[1]
  s_Y_x = gamma[2]
  rho_Y = gamma[3]


  # Sigma functions
  sigma_M1 = sigma_fun(x, s_M_0, s_M_x, rho_M)
  sigma_M2 = sigma_fun(x + 1, s_M_0, s_M_x, rho_M)

  sigma_Y1 = sigma_fun(x, s_Y_0, s_Y_x, rho_Y)
  sigma_Y2 = sigma_fun(x + 1, s_Y_0, s_Y_x, rho_Y)

  # Return vector of arguments to Phi
  return(c(eta, zeta, a_x, b_m, b_x, sigma_M2, sigma_Y2, sigma_M1, sigma_Y1))
}


# Gradient of xi
# Specifically, we need the gradient of this transformation so we can use the delta-method. Rather than giving each partial its own function, I will compute the gradient of each coordinate in the range.
# Notation: a, theta, b, gamma -> eta, zeta, a_x, b_m, b_x, s_U2, s_V2, s_U1, s_V1, where
# a = (a_0, a_x, A_2) is the fixed effects for M
# theta = (theta_1, ..., theta_q) is the random effects for M
# b = (b_0, b_m, b_x, B_3) is the fixed effects for Y
# gamma = (gamma_1, ..., gamma_p) is the random effects for Y
#
# eta is the linear predictor for M
# zeta is the linear predictor for Y
# a_x, b_m and b_x are fixed effects
# s_U2 = sigma_U(x+1), s_V2 = sigma_V(x+1), s_U1 = sigma_U(x), s_V1 = sigma_V(x). See the notes from B & B for details
## All these gradients have been validated against a finite difference approximation

grad_eta = function(a, theta, b, gamma, x, W){
  q_M = length(theta)
  q_Y = length(gamma)
  len_W = length(W)
  output = c(1,
              x,
              W,
              rep(0, times=q_M),   # theta
              rep(0, times=3 + len_W),  # b_0, b_m, b_x, B_3
              rep(0, times = q_Y))  # gamma
}

grad_zeta = function(a, theta, b, gamma, x, W){
  q_M = length(theta)
  q_Y = length(gamma)
  len_W = length(W)
  output = c(rep(0, times=2 + len_W),  # a_0, a_x, A_2
              rep(0, times = q_M),   # theta
              1,  # beta_0
              0,  # beta_M
              x,  # beta_X
              W,  # B_3
              rep(0, times = q_Y))  # gamma
}

grad_a_x = function(a, theta, b, gamma, x, W){
  q_M = length(theta)
  q_Y = length(gamma)
  len_W = length(W)
  output = c(0, # a_0
              1, #a_x
              rep(0, times = len_W),  # A_2
              rep(0, times = q_M),   # theta
              rep(0, times=3),  # b_0, b_m, b_x
              rep(0, times = len_W),  # B_3
              rep(0, times = q_Y))  # gamma
}

grad_b_m = function(a, theta, b, gamma, x, W){
  q_M = length(theta)
  q_Y = length(gamma)
  len_W = length(W)
  output = c(rep(0, times=2 + len_W),  # a_0, a_x, A_2
              rep(0, times = q_M),   # theta
              0, # b_0
              1, # b_m
              0, # b_x
              rep(0, times = len_W),  # B_3
              rep(0, times = q_Y))  # gamma
}

grad_b_x = function(a, theta, b, gamma, x, W){
  q_M = length(theta)
  q_Y = length(gamma)
  len_W = length(W)
  output = c(rep(0, times=2 + len_W),  # a_0, a_x, A_2
              rep(0, times = q_M),   # theta
              0, # b_0
              0, # b_m
              1, # b_x
              rep(0, times = len_W),  # B_3
              rep(0, times = q_Y))  # gamma
}




grad_b_m = function(a, theta, b, gamma, x, W){
  q_M = length(theta)
  q_Y = length(gamma)
  len_W = length(W)
  output = c(rep(0, times=2 + len_W),  # a_0, a_x, A_2
              rep(0, times = q_M),   # theta
              0, # b_0
              1, # b_m
              0, # b_x
              rep(0, times = len_W),  # B_3
              rep(0, times = q_Y))  # gamma
}




# Dispersion parameters

## Computes the sigma function, e.g., sigma_u(x)
sigma_fun = function(x, s_0, s_x, rho){
  A = s_0^2
  B = 2 * x * rho * s_0 * s_x
  C = x^2 * s_x^2

  return(sqrt(A + B + C))
}

## Partial derivatives of sigma_fun
### All have been validated against a simple finite difference approximation

d_sigma_d_s0 = function(x, s_0, s_x, rho){
  return((s_0 + x * rho * s_x) / sigma_fun(x, s_0, s_x, rho))
}

d_sigma_d_sx = function(x, s_0, s_x, rho){
  return((s_x * x^2 + x * rho * s_0) / sigma_fun(x, s_0, s_x, rho))
}

d_sigma_d_rho = function(x, s_0, s_x, rho){
  return((s_0 * s_x * x) / sigma_fun(x, s_0, s_x, rho))
}



# x = 1
# s_0 = 2
# s_x = 3
# rho = 4
#
# e = 0.000001
#
# (sigma_fun(x+1, s_0 + e, s_x, rho) - sigma_fun(x+1, s_0, s_x, rho)) / e
# d_sigma_d_s0(x+1, s_0, s_x, rho)
#
# d_sigma_d_rho(x+1, s_0, s_x, rho)







## Gradient of sigma_U(x+1)
grad_s_U2 = function(a, theta, b, gamma, x, W){
  # Note: All ambiguous parameters in this function are for the M model. This affects s_0, s_x and rho

  s_0 = theta[1]
  s_x = theta[2]
  rho = theta[3]

  q_M = length(theta)
  q_Y = length(gamma)
  len_W = length(W)

  # Compute partials
  d_s_U2_d_s0 = d_sigma_d_s0(x+1, s_0, s_x, rho)
  d_s_U2_d_sx = d_sigma_d_sx(x+1, s_0, s_x, rho)
  d_s_U2_d_rho = d_sigma_d_rho(x+1, s_0, s_x, rho)

  # Stack partials
  d_s_U2_d_theta = c(d_s_U2_d_s0, d_s_U2_d_sx, d_s_U2_d_rho)

  # Construct gradient
  grad = c(rep(0, times=2 + len_W),  # a_0, a_x, A_2
            d_s_U2_d_theta,  # theta
            rep(0, times=3),  # b_0, b_m, b_x
            rep(0, times = len_W),  # B_3
            rep(0, times = q_Y)  # gamma
  )
  return(grad)
}

## Gradient of sigma_V(x+1)
grad_s_V2 = function(a, theta, b, gamma, x, W){
  # Note: All ambiguous parameters in this function are for the Y model. This affects s_0, s_x and rho

  s_0 = gamma[1]
  s_x = gamma[2]
  rho = gamma[3]

  q_M = length(theta)
  q_Y = length(gamma)
  len_W = length(W)

  # Compute partials
  d_s_V2_d_s0 = d_sigma_d_s0(x+1, s_0, s_x, rho)
  d_s_V2_d_sx = d_sigma_d_sx(x+1, s_0, s_x, rho)
  d_s_V2_d_rho = d_sigma_d_rho(x+1, s_0, s_x, rho)

  # Stack partials
  d_s_V2_d_gamma = c(d_s_V2_d_s0, d_s_V2_d_sx, d_s_V2_d_rho)

  # Construct gradient
  grad = c(rep(0, times=2 + len_W),  # a_0, a_x, A_2
            rep(0, times = q_M),  # theta
            rep(0, times=3),  # b_0, b_m, b_x
            rep(0, times = len_W),  # B_3
            d_s_V2_d_gamma  # gamma
  )
  return(grad)
}



## Gradient of sigma_U(x)
grad_s_U1 = function(a, theta, b, gamma, x, W){
  # Note: All ambiguous parameters in this function are for the M model

  s_0 = theta[1]
  s_x = theta[2]
  rho = theta[3]

  q_M = length(theta)
  q_Y = length(gamma)
  len_W = length(W)

  # Compute partials
  d_s_U1_d_s0 = d_sigma_d_s0(x, s_0, s_x, rho)
  d_s_U1_d_sx = d_sigma_d_sx(x, s_0, s_x, rho)
  d_s_U1_d_rho = d_sigma_d_rho(x, s_0, s_x, rho)

  # Stack partials
  d_s_U1_d_theta = c(d_s_U1_d_s0, d_s_U1_d_sx, d_s_U1_d_rho)

  # Construct gradient
  grad = c(rep(0, times=2 + len_W),  # a_0, a_x, A_2
                d_s_U1_d_theta,  # theta
                rep(0, times=3),  # b_0, b_m, b_x
                rep(0, times = len_W),  # B_3
                rep(0, times = q_Y)  # gamma
                )
  return(grad)
}

## Gradient of sigma_V(x)
grad_s_V1 = function(a, theta, b, gamma, x, W){
  # Note: All ambiguous parameters in this function are for the Y model

  s_0 = gamma[1]
  s_x = gamma[2]
  rho = gamma[3]

  q_M = length(theta)
  q_Y = length(gamma)
  len_W = length(W)

  # Compute partials
  d_s_V1_d_s0 = d_sigma_d_s0(x, s_0, s_x, rho)
  d_s_V1_d_sx = d_sigma_d_sx(x, s_0, s_x, rho)
  d_s_V1_d_rho = d_sigma_d_rho(x, s_0, s_x, rho)

  # Stack partials
  d_s_V1_d_gamma = c(d_s_V1_d_s0, d_s_V1_d_sx, d_s_V1_d_rho)

  # Construct gradient
  grad = c(rep(0, times=2 + len_W),  # a_0, a_x, A_2
                rep(0, times = q_M),  # theta
                rep(0, times=3),  # b_0, b_m, b_x
                rep(0, times = len_W),  # B_3
                d_s_V1_d_gamma  # gamma
                )
  return(grad)
}



# Stack all gradients
grad_xi = function(a, theta, b, gamma, x, W){
  output = rbind(grad_eta(a, theta, b, gamma, x, W),
                grad_zeta(a, theta, b, gamma, x, W),
                grad_a_x(a, theta, b, gamma, x, W),
                grad_b_m(a, theta, b, gamma, x, W),
                grad_b_x(a, theta, b, gamma, x, W),
                grad_s_U2(a, theta, b, gamma, x, W),
                grad_s_V2(a, theta, b, gamma, x, W),
                grad_s_U1(a, theta, b, gamma, x, W),
                grad_s_V1(a, theta, b, gamma, x, W))
  rownames(output) = NULL
  colnames(output) = NULL

  return(output)
}























# Finally, we can combine the gradient of Phi with the gradient of the map from GLMM parameters to arguments of Phi. This gives us the gradient of Phi in terms of the parameters for which we have a covariance matrix.

Phi_of_xi = function(a, theta, b, gamma, x, W){
  Phi_args = xi(a, theta, b, gamma, x, W)
  output = do.call(Phi, as.list(Phi_args)) # Pass a vector of arguments to the function Phi
  return(output)
}




d_Phi_d_GLMM_pars = function(a, theta, b, gamma, x, W){
  Phi_args = xi(a, theta, b, gamma, x, W)
  d_xi = grad_xi(a, theta, b, gamma, x, W)

  d_Phi = do.call(grad_Phi, as.list(Phi_args))  # Pass a vector of arguments to the function Phi

  output = d_Phi %*% d_xi
  return(output)
}










#### Verify that gradients of xi, Phi and Phi(xi(.)) are correct ####


# # Finite difference verification for xi
#
# e = 0.00001
#
#
# a_0 = 1
# a_x = 2
# A_2 = c(3, 4, 5)
#
# b_0 = 6
# b_m = 7
# b_x = 8
# B_3 = c(9, 10, 11)
#
# theta_0 = 1.5
# theta_x = 2.5
# theta_rho = 0.7
#
# gamma_0 = 3.5
# gamma_x = 4.5
# gamma_rho = 0.8 + e
#
# x_ref = 1
# W_ref = c(2,3,4)
#
#
#
# a = c(a_0, a_x, A_2)
# theta = c(theta_0, theta_x, theta_rho)
# b = c(b_0, b_m, b_x, B_3)
# gamma = c(gamma_0, gamma_x, gamma_rho)
#
# # a_ref = c(a_0, a_x, A_2)
# # theta_ref = c(theta_0, theta_x, theta_rho)
# # b_ref = c(b_0, b_m, b_x, B_3)
# # gamma_ref = c(gamma_0, gamma_x, gamma_rho)
#
#
#
#
# grad_xi(a_ref, theta_ref, b_ref, gamma_ref, x_ref, W_ref)
# (xi(a, theta, b, gamma, x_ref, W_ref) - xi(a_ref, theta_ref, b_ref, gamma_ref, x_ref, W_ref)) / e




# # Finite difference verification for Phi
# # a_ref = 1
# # b_ref = 1
# # a1_ref = 1
# # b1_ref = 1
# # b2_ref = 1
# # s1_ref = 0.5
# # s2_ref = 0.5
# # s3_ref = 0.5
# # s4_ref = 0.5
#
#
#
# # e = 0.0000001
# e = sqrt(.Machine$double.eps)
#
#
#
# eta = 1 + e
# zeta = 1
# a1 = 1
# b1 = 1
# b2 = 1
# s1 = 0.5
# s2 = 0.5
# s3 = 0.5
# s4 = 0.5
#
# # eta_ref = eta
# # zeta_ref = zeta
# # a1_ref = a1
# # b1_ref = b1
# # b2_ref = b2
# # s1_ref = s1
# # s2_ref = s2
# # s3_ref = s3
# # s4_ref = s4
#
#
#
# (Phi(eta, zeta, a1, b1, b2, s1, s2, s3, s4) - Phi(eta_ref, zeta_ref, a1_ref, b1_ref, b2_ref, s1_ref, s2_ref, s3_ref, s4_ref)) / e
#
# grad_Phi(eta_ref, zeta_ref, a1_ref, b1_ref, b2_ref, s1_ref, s2_ref, s3_ref, s4_ref)



# # Finite difference verification for Phi(xi(.)) using the numDeriv package
# ## Manual version follows after. The package is more efficient, albeit took a bit of time to understand.
#
#
# library(numDeriv)
#
# a_0 = 1
# a_x = 2
# A_2 = c(3, 4, 5)
#
# b_0 = 6
# b_m = 7
# b_x = 8
# B_3 = c(9, 10, 11)
#
# theta_0 = 1.5
# theta_x = 2.5
# theta_rho = 0.7
#
# gamma_0 = 3.5
# gamma_x = 4.5
# gamma_rho = 0.8
#
# x_ref = 1
# W_ref = c(2,3,4)
#
# x_ref = 0
# W_ref = c(0, 0, 0)
#
#
#
# a = c(a_0, a_x, A_2)
# theta = c(theta_0, theta_x, theta_rho)
# b = c(b_0, b_m, b_x, B_3)
# gamma = c(gamma_0, gamma_x, gamma_rho)
#
# # a_ref = c(a_0, a_x, A_2)
# # theta_ref = c(theta_0, theta_x, theta_rho)
# # b_ref = c(b_0, b_m, b_x, B_3)
# # gamma_ref = c(gamma_0, gamma_x, gamma_rho)
#
# arg_vec = c(a_0, a_x, A_2[1], A_2[2], A_2[3], theta_0, theta_x, theta_rho, b_0, b_m, b_x, B_3[1], B_3[2], B_3[3], gamma_0, gamma_x, gamma_rho)
#
#
# my_fun = function(args){
#   a = args[1:5]
#   theta = args[6:8]
#   b = args[9:14]
#   gamma = args[15:17]
#
#   return(Phi_of_xi(a, theta, b, gamma, x, W))
#
# }
#
# jacobian(my_fun, arg_vec)
# d_Phi_d_GLMM_pars(a, theta, b, gamma, x_ref, W_ref)
#
#
# #
# #
# #
# # e = 0.00001
# # e = sqrt(.Machine$double.eps)
# #
# # (Phi_of_xi(a, theta, b, gamma, x_ref, W_ref) - Phi_of_xi(a_ref, theta_ref, b_ref, gamma_ref, x_ref, W_ref)) / e
# # d_Phi_d_GLMM_pars(a, theta, b, gamma, x_ref, W_ref)

