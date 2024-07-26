
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
