

test_that("psi works", {
  expect_equal(psi(0, 1), 0.5)
  expect_equal(psi(1,0), 1/(1+exp(-1)))
  expect_equal(psi(1,1), 1 - 0.5*exp(-1/2))
})


test_that("d1_psi works", {
  expect_equal(d1_psi(1, 1), numDeriv::grad(psi, 1, sigma=1), tolerance = 1e-6)
  expect_equal(d1_psi(1, 0), numDeriv::grad(psi, 1, sigma=0), tolerance = 1e-6)
  expect_equal(d1_psi(0, 1), numDeriv::grad(psi, 0, sigma=1), tolerance = 1e-6)
  expect_equal(d1_psi(0, 0), numDeriv::grad(psi, 0, sigma=0), tolerance = 1e-6)
})


#! These are wrong. START HERE
test_that("d2_psi works", {
  expect_equal(d2_psi(1, 1), numDeriv::grad(psi, 1, mu=1), tolerance = 1e-6)
  expect_equal(d2_psi(1, 0), numDeriv::grad(psi, 0, mu=1), tolerance = 1e-6)
  expect_equal(d2_psi(0, 1), numDeriv::grad(psi, 1, mu=0), tolerance = 1e-6)
  expect_equal(d2_psi(0, 0), numDeriv::grad(psi, 0, mu=0), tolerance = 1e-6)
})
