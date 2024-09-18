
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

    ENC_zeros = ENC(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M)
    ENC_effs = ENC(x, x_m, w, b_Y, theta_Y[Y_ind], b_M, theta_M[M_ind], which_REs = this_REs)

    expect_equal(ENC_zeros, ENC_effs,
                 label = paste0(Y_RE, " with ", M_RE))
  }
})


#* Dimensions of gradients

test_that("Dimensions of gradients are correct for subsets of REs", {

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

      this_num_pars = length(b_Y) + length(b_M) + length(this_theta_Y) + length(this_theta_M)


      ## grad_mu_Y
      this_grad_mu_Y = grad_mu_Y(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs)
      expect_equal(length(this_grad_mu_Y), this_num_pars)

      ## grad_mu_M
      this_grad_mu_M = grad_mu_M(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs)
      expect_equal(length(this_grad_mu_M), this_num_pars)

      ## grad_b_Y_M
      this_grad_b_Y_M = grad_b_Y_M(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs)
      expect_equal(length(this_grad_b_Y_M), this_num_pars)

      ## grad_gamma_Y
      this_grad_gamma_Y = grad_gamma_Y(m, x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs)
      expect_equal(length(this_grad_gamma_Y), this_num_pars)

      ## grad_gamma_M
      this_grad_gamma_M = grad_gamma_M(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs)
      expect_equal(length(this_grad_gamma_M), this_num_pars)
    }
  }
})






# Values of gradients

## Non-trivial values for the b's and theta's. Former based on output from another MC study. Latter chosen arbitrarily.
## Crucially, no parameters are equal to zero.
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


test_that("Values of grad_mu_Y are correct for subsets of REs", {

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

      this_num_pars = length(b_Y) + length(b_M) + length(this_theta_Y) + length(this_theta_M)


      ## grad_mu_Y
      params = c(b_Y, this_theta_Y, b_M, this_theta_M)

      test_mu_Y = function(x_val, w, len_theta_Y, params){
        b_Y = params[1:5]
        theta_Y = params[6:(5 + len_theta_Y)]
        b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
        theta_M = params[(10 + len_theta_Y):length(params)]

        mu_Y = as.numeric(b_Y[1] + x_val * b_Y[2] + w %*% b_Y[4:length(b_Y)])
        return(mu_Y)
      }

      for(x in c(0,1)){
        for(x_m in c(0,1)){
          
          expect_equal(grad_mu_Y(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_mu_Y, params, x_val=x, w=w, len_theta_Y = length(this_theta_Y)))
        }
      
      }

    }
  }
})

test_that("Values of grad_mu_M are correct for subsets of REs", {

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

      this_num_pars = length(b_Y) + length(b_M) + length(this_theta_Y) + length(this_theta_M)


      ## grad_mu_Y
      params = c(b_Y, this_theta_Y, b_M, this_theta_M)

      test_mu_M = function(x_val, w, len_theta_Y, params){
        b_Y = params[1:5]
        theta_Y = params[6:(5 + len_theta_Y)]
        b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
        theta_M = params[(10 + len_theta_Y):length(params)]

        mu_M = as.numeric(b_M[1] + x_val * b_M[2] + w %*% b_M[3:length(b_M)])
        return(mu_M)
      }

      for(x in c(0,1)){
        for(x_m in c(0,1)){
          expect_equal(grad_mu_M(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_mu_M, params, x_val=x_m, w=w, len_theta_Y = length(this_theta_Y)))
        }
      }

    }
  }
})


test_that("Values of grad_b_Y_M are correct for subsets of REs", {

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

      this_num_pars = length(b_Y) + length(b_M) + length(this_theta_Y) + length(this_theta_M)


      ## grad_mu_Y
      params = c(b_Y, this_theta_Y, b_M, this_theta_M)

      test_b_Y_M = function(len_theta_Y, params){
        b_Y = params[1:5]
        theta_Y = params[6:(5 + len_theta_Y)]
        b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
        theta_M = params[(10 + len_theta_Y):length(params)]

        b_Y_M = b_Y[3]
        return(b_Y_M)
      }

      for(x in c(0,1)){
        for(x_m in c(0,1)){
          
          expect_equal(grad_b_Y_M(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_b_Y_M, params, len_theta_Y = length(this_theta_Y)))
        }
      
      }

    }
  }
})

test_that("Values of grad_gamma_Y are correct for subsets of REs", {

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

      this_num_pars = length(b_Y) + length(b_M) + length(this_theta_Y) + length(this_theta_M)


      ## grad_mu_Y
      params = c(b_Y, this_theta_Y, b_M, this_theta_M)

      test_gamma_Y = function(x_val, m_val, this_REs, len_theta_Y, params){
        b_Y = params[1:5]
        theta_Y = params[6:(5 + len_theta_Y)]
        b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
        theta_M = params[(10 + len_theta_Y):length(params)]

        Y_vec = Y_vec_gamma(x_val, m, this_REs)
        gamma_Y = theta2gamma(Y_vec, theta_Y)

        return(gamma_Y)
      }

      for(x in c(0,1)){
        for(x_m in c(0,1)){
          for(m in c(0,1)){
            expect_equal(grad_gamma_Y(m, x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_gamma_Y, params, x_val=x, m_val=m, this_REs = this_REs, len_theta_Y = length(this_theta_Y)))
          }
        }
      }

    }
  }
})


test_that("Values of grad_gamma_M are correct for subsets of REs", {

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

      this_num_pars = length(b_Y) + length(b_M) + length(this_theta_Y) + length(this_theta_M)


      ## grad_mu_Y
      params = c(b_Y, this_theta_Y, b_M, this_theta_M)

      test_gamma_M = function(x_m, this_REs, len_theta_Y, params){
        b_Y = params[1:5]
        theta_Y = params[6:(5 + len_theta_Y)]
        b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
        theta_M = params[(10 + len_theta_Y):length(params)]

        M_vec = M_vec_gamma(x_m, this_REs)
        gamma_M = theta2gamma(M_vec, theta_M)

        return(gamma_M)
      }

      for(x in c(0,1)){
        for(x_m in c(0,1)){
          expect_equal(grad_gamma_M(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_gamma_M, params, x_m=x_m, this_REs = this_REs, len_theta_Y = length(this_theta_Y)))
        }
      }

    }
  }
})


test_that("Values of grad_psi_Y are correct for subsets of REs", {

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

      this_num_pars = length(b_Y) + length(b_M) + length(this_theta_Y) + length(this_theta_M)


      ## grad_mu_Y
      params = c(b_Y, this_theta_Y, b_M, this_theta_M)

      test_psi_Y = function(x_val, x_m, m_val, w, this_REs, len_theta_Y, params){
        b_Y = params[1:5]
        theta_Y = params[6:(5 + len_theta_Y)]
        b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
        theta_M = params[(10 + len_theta_Y):length(params)]

        mu_Y = as.numeric(b_Y[1] + x * b_Y[2] + w %*% b_Y[4:length(b_Y)])

        Y_vec = Y_vec_gamma(x_val, m_val, this_REs)
        gamma_Y = theta2gamma(Y_vec, theta_Y)

        this_psi = psi(mu_Y + m_val * b_Y[3], gamma_Y)

        return(this_psi)
      }

      for(x in c(0,1)){
        for(x_m in c(0,1)){
          for(m in c(0,1)){
            expect_equal(grad_psi_Y(m, x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_psi_Y, params, x_val = x, x_m=x_m, m_val = m, w=w, this_REs = this_REs, len_theta_Y = length(this_theta_Y)),
            tolerance = 1e-6)
          }
        }
      }

    }
  }
})


test_that("Values of grad_psi_M are correct for subsets of REs", {

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

      this_num_pars = length(b_Y) + length(b_M) + length(this_theta_Y) + length(this_theta_M)


      ## grad_mu_Y
      params = c(b_Y, this_theta_Y, b_M, this_theta_M)


      test_psi_M = function(x_m, m_val, w, this_REs, len_theta_Y, params){
        b_Y = params[1:5]
        theta_Y = params[6:(5 + len_theta_Y)]
        b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
        theta_M = params[(10 + len_theta_Y):length(params)]

        mu_M = as.numeric(b_M[1] + x_m * b_M[2] + w %*% b_M[3:length(b_M)])
  
        M_vec = M_vec_gamma(x_m, this_REs)
        gamma_M = theta2gamma(M_vec, theta_M)

        this_psi = psi(mu_M * (2*m_val - 1), gamma_M)  # (2*m_val - 1) = 1 if m=1, -1 if m=0

        return(this_psi)
      }

      for(x in c(0,1)){
        for(x_m in c(0,1)){
          for(m in c(0,1)){
            expect_equal(grad_psi_M(m, x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_psi_M, params, x_m=x_m, m_val = m, w=w, this_REs = this_REs, len_theta_Y = length(this_theta_Y)),
            tolerance = 1e-6)
          }
        }
      }

    }
  }
})


test_that("Values of grad_ENC are correct for subsets of REs", {

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

      this_num_pars = length(b_Y) + length(b_M) + length(this_theta_Y) + length(this_theta_M)


      ## grad_mu_Y
      params = c(b_Y, this_theta_Y, b_M, this_theta_M)


      test_ENC = function(x_val, x_m, w, this_REs, len_theta_Y, params){
        b_Y = params[1:5]
        theta_Y = params[6:(5 + len_theta_Y)]
        b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
        theta_M = params[(10 + len_theta_Y):length(params)]

        this_ENC = ENC(x_val, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs = this_REs)

        return(this_ENC)
      }

      for(x in c(0,1)){
        for(x_m in c(0,1)){
          expect_equal(grad_ENC(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_ENC, params, x_val = x, x_m=x_m, w=w, this_REs = this_REs, len_theta_Y = length(this_theta_Y)),
          tolerance = 1e-6)
        }
      }

    }
  }
})


# #* Test a specific value of the gradient
# test_that("Gradient of ENC is correct at simulated value of Theta",{

#   test_ENC = function(x_val, x_m, w, this_REs, len_theta_Y, params){
#     b_Y = params[1:5]
#     theta_Y = params[6:(5 + len_theta_Y)]
#     b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
#     theta_M = params[(10 + len_theta_Y):length(params)]

#     this_ENC = ENC(x_val, x_m, w, b_Y, theta_Y, b_M, theta_M, which_REs = this_REs)

#     return(this_ENC)
#   }



#   this_b_Y = c(0.0376828219852018, 0.966486302988689, 1.99644760563721, -0.00556557712859059, 0.000826754128449799)
#   this_theta_Y = c(0.703069437670255, 0.0400470080326032, 1.02976021614724)
#   this_b_M = c(-0.0990439890654785, 1.76353928991247, 0.0128566136999183, 0.00711746366915989)
#   this_theta_M  = c(0.994339635253282, 0.480677812602227, 1.85947197356182)
#   params = c(this_b_Y, this_theta_Y, this_b_M, this_theta_M)

#   this_REs = c("Y.Int", "Y.X", "M.All")

#   for(x in c(0,1)){
#     for(x_m in c(0,1)){
#       expect_equal(grad_ENC(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_ENC, params, x_val = x, x_m=x_m, w=w, this_REs = this_REs, len_theta_Y = length(this_theta_Y)),
#       tolerance = 1e-6)
#     }
#   }
# })

# grad_hat = grad_ENC(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs)
# grad_num = numDeriv::grad(test_ENC, params, x_val = x, x_m=x_m, w=w, this_REs = this_REs, len_theta_Y = length(this_theta_Y))

# (err = grad_hat - grad_num)
# norm(err, "2") / norm(grad_num, "2")
# abs(err) / abs(grad_num)
# grad_num# Covariance matrix of ENC
# err

# par_names = c("b_Y_Int", "b_Y_X", "b_Y_M", "b_Y_W1", "b_Y_W2", "theta_Y_Int", "theta_Y_corr", "theta_Y_X", "b_M_Int", "b_M_X", "b_M_W1", "b_M_W2", "theta_M_Int", "theta_M_corr", "theta_M_X")

# data_errs = data.frame(par = par_names, est = grad_hat, num = grad_num, err = err)
# # Repeat the previous line with names in quotations



# test_psi_M = function(x_val, x_m, m_val, w, this_REs, len_theta_Y, params){
#   b_Y = params[1:5]
#   theta_Y = params[6:(5 + len_theta_Y)]
#   b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
#   theta_M = params[(10 + len_theta_Y):length(params)]

#   mu_M = as.numeric(b_M[1] + x_m * b_M[2] + w %*% b_M[3:length(b_M)])

#   M_vec = M_vec_gamma(x_m, this_REs)
#   gamma_M = theta2gamma(M_vec, theta_M)

#   this_psi = psi(mu_M * (2*m - 1), gamma_M)  # (1 - 2*m) = 1 if m=1, -1 if m=0

#   return(this_psi)
# }

# for(x in c(0,1)){
#   for(x_m in c(0,1)){
#     for(m in c(0,1)){
#       expect_equal(grad_psi_M(m, x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_psi_M, params, x_val = x, x_m=x_m, m_val = m, w=w, this_REs = this_REs, len_theta_Y = length(this_theta_Y)),
#       tolerance = 1e-6)
#     }
#   }
# }

# grad_psi_M(m, x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs)
# which(numDeriv::grad(test_psi_M, params, x_val = x, x_m=x_m, m_val = m, w=w, this_REs = this_REs, len_theta_Y = length(this_theta_Y)) != 0)





# test_psi_Y = function(x_val, x_m, m_val, w, this_REs, len_theta_Y, params){
#         b_Y = params[1:5]
#         theta_Y = params[6:(5 + len_theta_Y)]
#         b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
#         theta_M = params[(10 + len_theta_Y):length(params)]

#         mu_Y = as.numeric(b_Y[1] + x_val * b_Y[2] + w %*% b_Y[4:length(b_Y)])

#         Y_vec = Y_vec_gamma(x_val, m_val, this_REs)
#         gamma_Y = theta2gamma(Y_vec, theta_Y)

#         this_psi = psi(mu_Y + m_val * b_Y[3], gamma_Y)

#         return(this_psi)
#       }

#       for(x in c(0,1)){
#         for(x_m in c(0,1)){
#           for(m in c(0,1)){
#             expect_equal(grad_psi_Y(m, x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_psi_Y, params, x_val = x, x_m=x_m, m_val = m, w=w, this_REs = this_REs, len_theta_Y = length(this_theta_Y)),
#             tolerance = 1e-6)
#           }
#         }
#       }



# test_gamma_M = function(x_m, this_REs, len_theta_Y, params){
#         b_Y = params[1:5]
#         theta_Y = params[6:(5 + len_theta_Y)]
#         b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
#         theta_M = params[(10 + len_theta_Y):length(params)]

#         M_vec = M_vec_gamma(x_m, this_REs)
#         gamma_M = theta2gamma(M_vec, theta_M)

#         return(gamma_M)
#       }

#       for(x in c(0,1)){
#         for(x_m in c(0,1)){
#           expect_equal(grad_gamma_M(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_gamma_M, params, x_m=x_m, this_REs = this_REs, len_theta_Y = length(this_theta_Y)))
#         }
#       }



# test_gamma_Y = function(x_val, m_val, this_REs, len_theta_Y, params){
#         b_Y = params[1:5]
#         theta_Y = params[6:(5 + len_theta_Y)]
#         b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
#         theta_M = params[(10 + len_theta_Y):length(params)]

#         Y_vec = Y_vec_gamma(x_val, m, this_REs)
#         gamma_Y = theta2gamma(Y_vec, theta_Y)

#         return(gamma_Y)
#       }

#       for(x in c(0,1)){
#         for(x_m in c(0,1)){
#           for(m in c(0,1)){
#             expect_equal(grad_gamma_Y(m, x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_gamma_Y, params, x_val=x, m_val=m, this_REs = this_REs, len_theta_Y = length(this_theta_Y)))
#           }
#         }
#       }



# test_b_Y_M = function(len_theta_Y, params){
#         b_Y = params[1:5]
#         theta_Y = params[6:(5 + len_theta_Y)]
#         b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
#         theta_M = params[(10 + len_theta_Y):length(params)]

#         b_Y_M = b_Y[3]
#         return(b_Y_M)
#       }

#       for(x in c(0,1)){
#         for(x_m in c(0,1)){
          
#           expect_equal(grad_b_Y_M(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_b_Y_M, params, len_theta_Y = length(this_theta_Y)))
#         }
      
#       }





# test_mu_M = function(x_val, w, len_theta_Y, params){
#         b_Y = params[1:5]
#         theta_Y = params[6:(5 + len_theta_Y)]
#         b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
#         theta_M = params[(10 + len_theta_Y):length(params)]

#         mu_M = as.numeric(b_M[1] + x_val * b_M[2] + w %*% b_M[3:length(b_M)])
#         return(mu_M)
#       }

#       for(x in c(0,1)){
#         for(x_m in c(0,1)){
#           expect_equal(grad_mu_M(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_mu_M, params, x_val=x_m, w=w, len_theta_Y = length(this_theta_Y)))
#         }
#       }



# test_mu_Y = function(x_val, w, len_theta_Y, params){
#         b_Y = params[1:5]
#         theta_Y = params[6:(5 + len_theta_Y)]
#         b_M = params[(6 + len_theta_Y):(9 + len_theta_Y)]
#         theta_M = params[(10 + len_theta_Y):length(params)]

#         mu_Y = as.numeric(b_Y[1] + x_val * b_Y[2] + w %*% b_Y[4:length(b_Y)])
#         return(mu_Y)
#       }

#       for(x in c(0,1)){
#         for(x_m in c(0,1)){
          
#           expect_equal(grad_mu_Y(x, x_m, w, b_Y, this_theta_Y, b_M, this_theta_M, which_REs = this_REs), numDeriv::grad(test_mu_Y, params, x_val=x, w=w, len_theta_Y = length(this_theta_Y)))
#         }
      
#       }

























## Note: This test depends on objects computed in test-Reg_Par_Covs.R
test_that("Joint covariance of ENC at all input levels is positive definite",{
  skip_on_cran()
  load("w_fit_Y_fit_M.RData")
  expect_true(all(eigen(all_covs_ENC(w, fit_Y, fit_M), symmetric=T, only.values = T)$values > 0))
})
