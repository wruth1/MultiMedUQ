
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
  expect_equal(all_MEs("diff", w, b_Y, theta_Y, b_M, theta_M), c(total_effect("diff", w, b_Y, theta_Y, b_M, theta_M),
    direct_effect("diff", w, b_Y, theta_Y, b_M, theta_M), indirect_effect("diff", w, b_Y, theta_Y, b_M, theta_M)))
})
