
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



#* Compare all_MEs vs individual effects

w = c(2,3)

b_Y = c(0,0,1,0,0)
theta_Y = c(sqrt(0.5), 0.5, 0, 1, 0.5, sqrt(0.5))

b_M = c(0,0,0,0)
theta_M = c(1, 0.5, 2)

## Check that general function matches specific functions
test_that("all_MEs matches individual effects", {
  expect_equal(unname(all_MEs_pars("diff", w, b_Y, theta_Y, b_M, theta_M)), c(total_effect("diff", w, b_Y, theta_Y, b_M, theta_M),
    direct_effect("diff", w, b_Y, theta_Y, b_M, theta_M), indirect_effect("diff", w, b_Y, theta_Y, b_M, theta_M)))
})









#* Subset of REs


## Non-trivial values for the b's and theta's. Former based on output from another MC study. Latter chosen arbitrarily.
## Crucially, no parameters are equal to zero.
b_Y = c(0.0376828219852018, 0.966486302988689, 1.99644760563721, -0.00556557712859059, 0.000826754128449799)
b_M = c(-0.0990439890654785, 1.76353928991247, 0.0128566136999183, 0.00711746366915989)

theta_Y = c(sqrt(0.5), 0.5, 0, 1, 0.5, sqrt(0.5))
theta_M = c(1, 0.5, 2)


scale = c("diff", "rat", "OR")

test_that("ME works with a subset of REs", {
  # Easy case: No effects
  expect_equal(all_MEs_pars(scale, c(0,0), rep(0, times=5), rep(0, times=6), rep(0, times=4), rep(0, times=3)), all_MEs_pars(scale, c(0,0), rep(0, times=5), rep(0, times=1), rep(0, times=4), rep(0, times=1), which_REs = c("Y.Int", "M.Int")))

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
    this_theta_Y[Y_ind] = theta_Y[Y_ind]
    # if(Y_RE == "Y.Int"){
    #   this_theta_Y[Y_ind] = sqrt(0.5)
    # } else if(Y_RE == "Y.X"){
    #   this_theta_Y[Y_ind] = 1
    # } else if(Y_RE == "Y.M"){
    #   this_theta_Y[Y_ind] = sqrt(0.5)
    # }

    this_theta_M = rep(0, times = 3)
    this_theta_M[M_ind] = theta_M[M_ind]
    # if(M_RE == "M.Int"){
    #   this_theta_M[M_ind] = 1
    # } else if(M_RE == "M.X"){
    #   this_theta_M[M_ind] = 2
    # }

    ME_zeros = all_MEs_pars(scale, w, b_Y, this_theta_Y, b_M, this_theta_M)
    ME_effs = all_MEs_pars(scale, w, b_Y, theta_Y[Y_ind], b_M, theta_M[M_ind], which_REs = this_REs)

    expect_equal(ME_zeros, ME_effs,
                 label = paste0(Y_RE, " with ", M_RE))
  }
})




#* Dimensions of gradients
all_scales = c("diff", "rat", "OR")

test_that("Dimensions of gradients are correct for subsets of REs", {

  Y_REs = c("Y.Int", "Y.X", "Y.M")
  M_REs = c("M.Int", "M.X")

  all_Y_sets = rje::powerSetCond(Y_REs)
  all_M_sets = rje::powerSetCond(M_REs)
  all_scale_sets = rje::powerSetCond(all_scales)

  for(i in seq_along(all_Y_sets)){
    for(j in seq_along(all_M_sets)){
      for(l in seq_along(all_scale_sets)){
        
        this_scale = all_scale_sets[[l]]

        this_Y_REs = all_Y_sets[[i]]
        this_M_REs = all_M_sets[[j]]
        this_REs = c(this_Y_REs, this_M_REs)

        this_theta_Y = make_theta(this_Y_REs)
        this_theta_M = make_theta(this_M_REs)


        this_grad_MEs = all_grad_MEs_pars(this_scale, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs)
        expect_equal(dim(this_grad_MEs), c(3 * length(this_scale), 4))
      
      }
    }
  }
})




#* Values of gradients
b_Y = c(0.0376828219852018, 0.966486302988689, 1.99644760563721, -0.00556557712859059, 0.000826754128449799)
b_M = c(-0.0990439890654785, 1.76353928991247, 0.0128566136999183, 0.00711746366915989)

make_theta = function(RE_names){
  num_REs = length(RE_names)

  if(num_REs == 1){
    return(1)
  } else if(num_REs == 2){
    return(c(2, 0.5, 3))
  } else if(num_REs == 3){
    return(c(4, 0.4, 0.3, 5, 0.2, 6))
  }
}

all_scales = c("diff", "rat", "OR")



test_that("Values of grads are correct for subsets of REs on difference scale", {

  scale = "diff"

  this_Y_REs = c("Y.Int", "Y.X", "Y.M")
  this_M_REs = c("M.Int", "M.X")
  this_REs = c(this_Y_REs, this_M_REs)

  this_theta_Y = make_theta(this_Y_REs)
  this_theta_M = make_theta(this_M_REs)

  this_ENCs = all_ENCs(w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs)

  this_num_pars = length(b_Y) + length(b_M) + length(this_theta_Y) + length(this_theta_M)


  # Total Effect
  test_TE_diff = function(ENCs, scale){
    get_ME(ENCs[1], ENCs[4], scale)[1]
  }

  expect_equal(grad_TE_diff(this_ENCs[1], this_ENCs[2], this_ENCs[3], this_ENCs[4]), numDeriv::grad(test_TE_diff, this_ENCs, scale = scale))


  # Direct Effect
  test_DE_diff = function(ENCs, scale){
    get_ME(ENCs[2], ENCs[4], scale)[1]
  }

  expect_equal(grad_DE_diff(this_ENCs[1], this_ENCs[2], this_ENCs[3], this_ENCs[4]), numDeriv::grad(test_DE_diff, this_ENCs, scale = scale))



  # Indirect Effect
  test_IE_diff = function(ENCs, scale){
    get_ME(ENCs[1], ENCs[3], scale)[1]
  }

  expect_equal(grad_IE_diff(this_ENCs[1], this_ENCs[2], this_ENCs[3], this_ENCs[4]), numDeriv::grad(test_IE_diff, this_ENCs, scale = scale))

})



#?  Order of output is {ENC(1,1), ENC(1,0), ENC(0,1), ENC(0,0)}.



test_that("Values of grads are correct for subsets of REs on ratio scale", {

  scale = "rat"

  Y_REs = c("Y.Int", "Y.X", "Y.M")
  M_REs = c("M.Int", "M.X")

  all_Y_sets = rje::powerSetCond(Y_REs)
  all_M_sets = rje::powerSetCond(M_REs)

  for(i in seq_along(all_Y_sets)){
    for(j in seq_along(all_M_sets)){

      this_Y_REs = all_Y_sets[[i]]
      this_M_REs = all_M_sets[[j]]
      this_REs = c(this_Y_REs, this_M_REs)

      this_theta_Y = make_theta(this_Y_REs)
      this_theta_M = make_theta(this_M_REs)

      this_ENCs = all_ENCs(w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs)



      # Total Effect
      test_TE_rat = function(ENCs, scale){
        get_ME(ENCs[1], ENCs[4], scale)[1]
      }

      expect_equal(grad_TE_rat(this_ENCs[1], this_ENCs[2], this_ENCs[3], this_ENCs[4]), numDeriv::grad(test_TE_rat, this_ENCs, scale = scale))


      # Direct Effect
      test_DE_rat = function(ENCs, scale){
        get_ME(ENCs[2], ENCs[4], scale)[1]
      }

      expect_equal(grad_DE_rat(this_ENCs[1], this_ENCs[2], this_ENCs[3], this_ENCs[4]), numDeriv::grad(test_DE_rat, this_ENCs, scale = scale))


      # Indirect Effect
      test_IE_rat = function(ENCs, scale){
        get_ME(ENCs[1], ENCs[3], scale)[1]
      }

      expect_equal(grad_IE_rat(this_ENCs[1], this_ENCs[2], this_ENCs[3], this_ENCs[4]), numDeriv::grad(test_IE_rat, this_ENCs, scale = scale))
    }
  }
})


test_that("Values of grads are correct for subsets of REs on odds-ratio scale", {

  scale = "OR"

  Y_REs = c("Y.Int", "Y.X", "Y.M")
  M_REs = c("M.Int", "M.X")

  all_Y_sets = rje::powerSetCond(Y_REs)
  all_M_sets = rje::powerSetCond(M_REs)

  for(i in seq_along(all_Y_sets)){
    for(j in seq_along(all_M_sets)){

      this_Y_REs = all_Y_sets[[i]]
      this_M_REs = all_M_sets[[j]]
      this_REs = c(this_Y_REs, this_M_REs)

      this_theta_Y = make_theta(this_Y_REs)
      this_theta_M = make_theta(this_M_REs)

      this_ENCs = all_ENCs(w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs)



      # Total Effect
      test_TE_or = function(ENCs, scale){
        get_ME(ENCs[1], ENCs[4], scale)[1]
      }

      expect_equal(grad_TE_or(this_ENCs[1], this_ENCs[2], this_ENCs[3], this_ENCs[4]), numDeriv::grad(test_TE_or, this_ENCs, scale = scale))


      # Direct Effect
      test_DE_or = function(ENCs, scale){
        get_ME(ENCs[2], ENCs[4], scale)[1]
      }

      expect_equal(grad_DE_or(this_ENCs[1], this_ENCs[2], this_ENCs[3], this_ENCs[4]), numDeriv::grad(test_DE_or, this_ENCs, scale = scale))


      # Indirect Effect
      test_IE_or = function(ENCs, scale){
        get_ME(ENCs[1], ENCs[3], scale)[1]
      }

      expect_equal(grad_IE_or(this_ENCs[1], this_ENCs[2], this_ENCs[3], this_ENCs[4]), numDeriv::grad(test_IE_or, this_ENCs, scale = scale))
    }
  }
})




#* Estimated covariance matrix should be positive definite
## Note: This test depends on objects computed in test-Reg_Par_Covs.R
test_that("Joint covariance of ENC at all input levels is positive definite",{
  skip_on_cran()
  load("w_fit_Y_fit_M.RData")
  expect_true(all(eigen(all_cov_MEs(w, fit_Y, fit_M), symmetric=T, only.values = T)$values > 0))
})
