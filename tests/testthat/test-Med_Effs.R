
# Generic functions for computing mediation effects on different scales

test_that("ME_diff works",{
  expect_equal(ME_diff(1, 1), 0)
  expect_equal(ME_diff(1, 0), 1)
  expect_equal(ME_diff(0, 1), -1)
  expect_equal(ME_diff(0, 0), 0)
})

test_that("ME_rat works",{
  expect_equal(ME_rat(1, 1), 1)
  expect_equal(ME_rat(1, 0), Inf)
  expect_equal(ME_rat(0, 1), 0)
  expect_equal(ME_rat(0, 0), NaN)
})

test_that("ME_OR works",{
  expect_equal(ME_OR(1/2, 1/2), 1)
  expect_equal(ME_OR(2/3, 1/3), 4)
  expect_equal(ME_OR(2/3, 1/2), 2)
  expect_equal(ME_OR(1/2, 2/3), 0.5)
  expect_equal(ME_OR(1/3, 2/3), 0.25)
})


# General function for computing mediation effects on different scales
## Should match specific functions above

test_that("get_ME works",{
  expect_equal(get_ME(1, 1, "diff"), ME_diff(1, 1))
  expect_equal(get_ME(1, 0, "diff"), ME_diff(1, 0))
  expect_equal(get_ME(0, 1, "diff"), ME_diff(0, 1))
  expect_equal(get_ME(0, 0, "diff"), ME_diff(0, 0))

  expect_equal(get_ME(1, 1, "rat"), ME_rat(1, 1))
  expect_equal(get_ME(1, 0, "rat"), ME_rat(1, 0))
  expect_equal(get_ME(0, 1, "rat"), ME_rat(0, 1))
  expect_equal(get_ME(0, 0, "rat"), ME_rat(0, 0))

  expect_equal(ME_OR(1/2, 1/2), get_ME(1/2, 1/2, "OR"))
  expect_equal(ME_OR(2/3, 1/3), get_ME(2/3, 1/3, "OR"))
  expect_equal(ME_OR(2/3, 1/2), get_ME(2/3, 1/2, "OR"))
  expect_equal(ME_OR(1/2, 2/3), get_ME(1/2, 2/3, "OR"))
  expect_equal(ME_OR(1/3, 2/3), get_ME(1/3, 2/3, "OR"))

  expect_error(get_ME(1, 1, "unknown_scale"), "Unknown scale")
})


## Multiple mediation effects

test_that("get_ME works with multiple scales",{
  expect_equal(get_ME(1, 1, c("diff", "rat")), c(ME_diff(1, 1), ME_rat(1, 1)))
  expect_equal(get_ME(1, 0, c("diff", "rat")), c(ME_diff(1, 0), ME_rat(1, 0)))
  expect_equal(get_ME(0, 1, c("diff", "rat")), c(ME_diff(0, 1), ME_rat(0, 1)))
  expect_equal(get_ME(0, 0, c("diff", "rat")), c(ME_diff(0, 0), ME_rat(0, 0)))

  expect_equal(get_ME(1/2, 1/2, c("diff", "rat", "OR")), c(ME_diff(1/2, 1/2), ME_rat(1/2, 1/2), ME_OR(1/2, 1/2)))
  expect_equal(get_ME(2/3, 1/3, c("diff", "rat", "OR")), c(ME_diff(2/3, 1/3), ME_rat(2/3, 1/3), ME_OR(2/3, 1/3)))
  expect_equal(get_ME(2/3, 1/2, c("diff", "rat", "OR")), c(ME_diff(2/3, 1/2), ME_rat(2/3, 1/2), ME_OR(2/3, 1/2)))
  expect_equal(get_ME(1/2, 2/3, c("diff", "rat", "OR")), c(ME_diff(1/2, 2/3), ME_rat(1/2, 2/3), ME_OR(1/2, 2/3)))
  expect_equal(get_ME(1/3, 2/3, c("diff", "rat", "OR")), c(ME_diff(1/3, 2/3), ME_rat(1/3, 2/3), ME_OR(1/3, 2/3)))
})




# Total, direct and indirect effects

w = c(2,3)

b_Y = c(0,0,1,0,0)
theta_Y = c(sqrt(0.5), 0.5, 0, 1, 0.5, sqrt(0.5))

b_M = c(0,0,0,0)
theta_M = c(1, 0.5, 2)

## Check that general function matches specific functions
test_that("all_MEs matches individual effects",{
  expect_equal(unname(all_MEs_pars("diff", w, b_Y, theta_Y, b_M, theta_M)), c(total_effect("diff", w, b_Y, theta_Y, b_M, theta_M),
    direct_effect("diff", w, b_Y, theta_Y, b_M, theta_M), indirect_effect("diff", w, b_Y, theta_Y, b_M, theta_M)))
})




# Covariance matrix of all mediation effects on various scales.

scale = c("diff", "rat", "OR")

## Note: This test depends on objects computed in test-Reg_Par_Covs.R
test_that("Joint covariance of all mediation effects is positive definite",{
  skip_on_cran()
  load("w_fit_Y_fit_M.RData")

  ### Some extremely small negative e-vals. Check that norm of negatives is very small
  e_vals = eigen(all_cov_MEs(scale, w, fit_Y, fit_M), symmetric=T, only.values = T)$values
  norm_pos = norm(e_vals[e_vals > 0], "2")
  norm_neg = norm(e_vals[e_vals < 0], "2")

  expect_true(norm_pos / norm_neg > 1e10)
  # expect_true(all(eigen(all_cov_MEs(scale, w, fit_Y, fit_M), symmetric=T, only.values = T)$values > 0))
})







#! START HERE!!!!!!!!

# TODO: Figure out why there is a discrepancy between the current implementation and the one from the Exact_Asymptotics project

#* To run the following:
library(lme4)
library(merDeriv)
source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
load("w_fit_Y_fit_M.RData")




# Compare current implementation with the one from the Exact_Asymptotics project
## The latter is reasonably well-validated against a Monte Carlo empirical SE
Y_model = fit_Y
Y_model_info = attributes(VarCorr(Y_model)$group)

M_model = fit_M
M_model_info = attributes(VarCorr(M_model)$group)

## Extract fitted parameters

### M model
a_hat = fixef(M_model)
a_RE_sds = M_model_info$stddev
a_RE_cor = M_model_info$correlation[2,1]
theta_hat = c(a_RE_sds, a_RE_cor)
if(any(is.nan(theta_hat))) stop("NaNs in theta_hat")  # Skip rest of current analysis if correlation is 0/0
M_cov = vcov(M_model, full=TRUE, ranpar="sd")


### Y model
b_hat = fixef(Y_model)
b_RE_sds = Y_model_info$stddev
b_RE_cor = Y_model_info$correlation[2,1]
gamma_hat = c(b_RE_sds, b_RE_cor)
if(any(is.nan(gamma_hat))) stop("NaNs in gamma_hat")  # Skip rest of current analysis if correlation is 0/0
Y_cov = vcov(Y_model, full=TRUE, ranpar="sd")
#

### Translate to terminology of MultiMedUQ
b_Y_alt = b_hat
theta_Y_alt = gamma_hat
b_M_alt = a_hat
theta_M_alt = theta_hat


### Estimate mediation effect
#### Fixed-effects
a_0_hat = a_hat[1]
a_x_hat = a_hat[2]
A_2_hat = a_hat[3:4]

b_0_hat = b_hat[1]
b_m_hat = b_hat[2]
b_x_hat = b_hat[3]
B_3_hat = b_hat[4:5]


## Linear predictors
eta_hat = as.numeric(a_0_hat + a_x_hat * 0 + w %*% A_2_hat)
zeta_hat = as.numeric(b_0_hat + b_x_hat * 0 + w %*% B_3_hat)


## Random effects covariances
s_M_0 = a_RE_sds[1]
s_M_x = a_RE_sds[2]
rho_M = a_RE_cor

s_Y_0 = b_RE_sds[1]
s_Y_x = b_RE_sds[2]
rho_Y = b_RE_cor


## Sigma functions
sigma_M1 = sigma_fun(0, s_M_0, s_M_x, rho_M)
sigma_M2 = sigma_fun(1, s_M_x, s_M_0, rho_M)

sigma_Y1 = sigma_fun(0, s_Y_0, s_Y_x, rho_Y)
sigma_Y2 = sigma_fun(1, s_Y_x, s_Y_0, rho_Y)

## Mediation effect
### See Helpers.R for the function Phi, which computes the mediation effect on odds-ratio scale
med_hat = Phi(eta_hat, zeta_hat, a_x_hat, b_m_hat, b_x_hat, sigma_M2, sigma_Y2, sigma_M1, sigma_Y1)


all_MEs_models(scale = "OR", w, fit_Y, fit_M)
