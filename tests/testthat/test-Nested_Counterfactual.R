
# Parameter values calibrated to get a particular value for the expected nested counterfactual
# Specifically, mu_Y = mu_M = 0, gamma_Y_1 = 1, b_Y_M = 1
# This gives an analytically tractable value for all of the psi integrals, 3/4 - exp(-1/2)/4. Note that psi(0, sigma) = 0.5 for any sigma, and psi(1,1) = 1 - exp(-1/2)/2

x = 0
x_m = 1
w = c(2,3)

b_Y = c(0,0,1,0,0)
theta_Y = c(sqrt(0.5), 0.5, 0, 1, 0.5, sqrt(0.5))

b_M = c(0,0,0,0)
theta_M = c(1, 0.5, 2)



test_that("ENC works", {
  # Easy case: No effects
  expect_equal(ENC(0, 0, c(0,0), rep(0, times=5), rep(0, times=6), rep(0, times=4), rep(0, times=3)), 0.5)

  # Harder case: Non-zero effects
  expect_equal(ENC(x, x_m, w, b_Y, theta_Y, b_M, theta_M), (3/4) - exp(-1/2)/4)
})

test_that("ENC works with a subset of REs", {
  # Easy case: No effects
  expect_equal(ENC(0, 0, c(0,0), rep(0, times=5), rep(0, times=6), rep(0, times=4), rep(0, times=3)), ENC(0, 0, c(0,0), rep(0, times=5), rep(0, times=1), rep(0, times=4), rep(0, times=1), which_REs = c("Y.Int", "M.Int")))

  # Harder case: Non-zero effects
  ## Loop over all pairs of single REs
  Y_REs = c("Y.Int", "Y.X", "Y.M")
  M_REs = c("M.Int", "M.X")
  RE_pairs = expand.grid(Y_REs, M_REs)

  Y_RE_inds = c(1, 4, 6)
  M_RE_inds = c(1,3)
  RE_ind_pairs = expand.grid(Y_RE_inds, M_RE_inds)

  for (i in seq_len(nrow(RE_pairs))){
    this_REs = as.character(unlist(RE_pairs[i,]))
    Y_RE = this_REs[1]
    M_RE = this_REs[2]

    Y_ind = RE_ind_pairs[i,1]
    M_ind = RE_ind_pairs[i,2]

    this_theta_Y = rep(0, times = 6)
    if(Y_RE == "Y.Int"){
      this_theta_Y[Y_ind] = sqrt(0.5)
    } else if(Y_RE == "Y.X"){
      this_theta_Y[Y_ind] = 1
    } else if(Y_RE == "Y.M"){
      this_theta_Y[Y_ind] = sqrt(0.5)
    }

    this_theta_M = rep(0, times = 3)
    if(M_RE == "M.Int"){
      this_theta_M[M_ind] = 1
    } else if(M_RE == "M.X"){
      this_theta_M[M_ind] = 2
    }

    ENC_zeros = ENC(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M)
    ENC_effs = ENC(x, x_m, w, b_Y, theta_Y[Y_ind], b_M, theta_M[M_ind], which_REs = this_REs)

    expect_equal(ENC_zeros, ENC_effs,
                 label = paste0(Y_RE, " with ", M_RE))
  }

  expect_equal(ENC(x, x_m, w, b_Y, c(sqrt(0.5), rep(0, times=5)), b_M, c(sqrt(0.5), 0, 0)),
               ENC(x, x_m, w, b_Y, sqrt(0.5), b_M, sqrt(0.5), which_REs = c("Y.Int", "M.Int")))
})





# Gradients

## Define versions psi for compatibility with numDeriv::grad()
test_psi_Y <- function(m_val, x_y, x_m, w, params){
  b_Y = params[1:5]
  theta_Y = params[6:11]
  b_M = params[12:15]
  theta_M = params[16:18]

  mu_Y = as.numeric(b_Y[1] + x_y * b_Y[2] + w %*% b_Y[4:length(b_Y)])
  gamma_Y = theta2gamma(c(1, x_y, m_val), theta_Y)

  return(psi(mu_Y + m_val*b_Y[3], gamma_Y))
}

test_psi_M <- function(m_val, x_y, x_m, w, params){
  b_Y = params[1:5]
  theta_Y = params[6:11]
  b_M = params[12:15]
  theta_M = params[16:18]

  mu_M = as.numeric(b_M[1] + x_m * b_M[2] + w %*% b_M[3:length(b_M)])
  gamma_M = theta2gamma(c(1, x_m), theta_M)

  return(psi((2*m_val - 1)*mu_M, gamma_M))
}



params = c(b_Y, theta_Y, b_M, theta_M)



test_that("grad_psi_Y works",{
  # Order of numeric arguments is m, x, x_m. See also arguments to numDeriv::grad()

  expect_equal(grad_psi_Y(1, 1, 1, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_Y, params, m_val=1, x_y=1, x_m=1, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_Y(1, 1, 0, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_Y, params, m_val=1, x_y=1, x_m=0, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_Y(1, 0, 1, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_Y, params, m_val=1, x_y=0, x_m=1, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_Y(1, 0, 0, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_Y, params, m_val=1, x_y=0, x_m=0, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_Y(0, 1, 1, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_Y, params, m_val=0, x_y=1, x_m=1, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_Y(0, 1, 0, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_Y, params, m_val=0, x_y=1, x_m=0, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_Y(0, 0, 1, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_Y, params, m_val=0, x_y=0, x_m=1, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_Y(0, 0, 0, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_Y, params, m_val=0, x_y=0, x_m=0, w=w), tolerance = 1e-6)

})


test_that("grad_psi_M works",{
  # Order of numeric arguments is m, x, x_m. See also arguments to numDeriv::grad()

  expect_equal(grad_psi_M(1, 1, 1, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_M, params, m_val=1, x_y=1, x_m=1, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_M(1, 1, 0, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_M, params, m_val=1, x_y=1, x_m=0, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_M(1, 0, 1, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_M, params, m_val=1, x_y=0, x_m=1, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_M(1, 0, 0, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_M, params, m_val=1, x_y=0, x_m=0, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_M(0, 1, 1, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_M, params, m_val=0, x_y=1, x_m=1, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_M(0, 1, 0, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_M, params, m_val=0, x_y=1, x_m=0, w=w), tolerance = 1e-6)
  expect_equal(grad_psi_M(0, 0, 1, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_psi_M, params, m_val=0, x_y=0, x_m=1, w=w), tolerance = 1e-6)
})



## Gradient of ENC

### Define versions of ENC for compatibility with numDeriv::grad()

test_ENC <- function(x_y, x_m, w, params){
  b_Y = params[1:5]
  theta_Y = params[6:11]
  b_M = params[12:15]
  theta_M = params[16:18]

  return(ENC(x_y, x_m, w, b_Y, theta_Y, b_M, theta_M))
}


test_that("grad_ENC works",{
  # Order of numeric arguments is x, x_m. See also arguments to numDeriv::grad()

  expect_equal(grad_ENC(1, 1, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_ENC, params, x_y=1, x_m=1, w=w), tolerance = 1e-6)
  expect_equal(grad_ENC(1, 0, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_ENC, params, x_y=1, x_m=0, w=w), tolerance = 1e-6)
  expect_equal(grad_ENC(0, 1, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_ENC, params, x_y=0, x_m=1, w=w), tolerance = 1e-6)
  expect_equal(grad_ENC(0, 0, w, b_Y, theta_Y, b_M, theta_M), numDeriv::grad(test_ENC, params, x_y=0, x_m=0, w=w), tolerance = 1e-6)
})




# Covariance matrix of ENC

## Note: This test depends on objects computed in test-Reg_Par_Covs.R
test_that("Joint covariance of ENC at all input levels is positive definite",{
  skip_on_cran()
  load("w_fit_Y_fit_M.RData")
  expect_true(all(eigen(all_covs_ENC(w, fit_Y, fit_M), symmetric=T, only.values = T)$values > 0))
})
