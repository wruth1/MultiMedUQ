

data(finance, package="smdata")

(fit_Y <- lme4::glmer(corr ~ targtop + cho + (cho + targtop - 1 | item), data = finance, family = binomial))
(fit_M <- lme4::glmer(targtop ~ cho + (cho - 1 | item), data = finance, family = binomial))

par_asymp_covs = all_pars_asymp_cov_mat(fit_Y, fit_M)
par_covs = all_pars_cov_mat(fit_Y, fit_M)

# Test that the asymptotic and re-scaled matrices are valid covariance matrices (i.e. is positive definite)
test_that("all_pars_cov_mat output is positive definite",{
  expect_true(all(eigen(par_asymp_covs, symmetric=T, only.values = T)$values > 0))
  expect_true(all(eigen(par_covs, symmetric=T, only.values = T)$values > 0))
})


