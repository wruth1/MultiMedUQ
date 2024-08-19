


# Test expand_REs

test_that("expand_REs works", {
  expect_equal(expand_REs(c("Y.all", "M.all")), c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X"))
  expect_equal(expand_REs(c("all")), c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X"))
})

