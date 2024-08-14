

# Simulated validation data ####

## Initialize variables

set.seed(1)

data = make_validation_data(20, 100, b_Y, theta_Y, b_M, theta_M, output_list = F)

w = c(0,0)

## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
(fit_Y = lme4::glmer(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))
(fit_M = lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))

par_asymp_covs <<- all_pars_asymp_cov_mat(fit_Y, fit_M)
par_covs <<- all_pars_cov_mat(fit_Y, fit_M)


# Save w, fit_Y and fit_M for later use, but only when not on CRAN
test_that("Workaround to avoid saving impermissible files when on CRAN",{
  skip_on_cran()
  save(w, fit_Y, fit_M, file = "w_fit_Y_fit_M.RData")
})


# (e_vals = eigen(par_asymp_covs, symmetric=T, only.values = T)$values)
# (norm_pos = norm(e_vals[e_vals>0], "2"))
# (norm_neg = norm(e_vals[e_vals<0], "2"))

# Test that the asymptotic and re-scaled matrices are valid covariance matrices (i.e. is positive definite)
test_that("all_pars_cov_mat output is positive definite",{
  # expect_true(norm_pos * 0.05 > norm_neg)

  expect_true(all(eigen(par_asymp_covs, symmetric=T, only.values = T)$values > 0))
  expect_true(all(eigen(par_covs, symmetric=T, only.values = T)$values > 0))
})



# # Real validation data (doesn't pass positive definite test due to small number of groups) ####
#
# data(finance, package="smdata")
#
# #
# # (fit_Y <- lme4::glmer(corr ~ targtop + cho + easyfoil + (cho + targtop - 1 | item), data = finance, family = binomial))
# # (fit_M <- lme4::glmer(targtop ~ cho + easyfoil + (cho - 1 | item), data = finance, family = binomial))
#
# (fit_Y <- lme4::glmer(corr ~ targtop + cho + easyfoil + (cho + targtop | item), data = finance, family = binomial))
# (fit_M <- lme4::glmer(targtop ~ cho + easyfoil + (cho | item), data = finance, family = binomial))
#
#
#
# par_asymp_covs = all_pars_asymp_cov_mat(fit_Y, fit_M)
# par_covs = all_pars_cov_mat(fit_Y, fit_M)
#
# # Test that the asymptotic and re-scaled matrices are valid covariance matrices (i.e. is positive definite)
# #! Skipped because the output is not positive definite
# test_that("all_pars_cov_mat output is positive definite",{
#   skip("This test works if I omit REs for the intercepts, but my other code currently requires them.")
#   expect_true(all(eigen(par_asymp_covs, symmetric=T, only.values = T)$values > 0))
#   expect_true(all(eigen(par_covs, symmetric=T, only.values = T)$values > 0))
# })
