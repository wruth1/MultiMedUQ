






# ---------------------------------------------------------------------------- #
#                       Functions to do MC-delta for MEs                       #
# ---------------------------------------------------------------------------- #

#! Try simulating a fixed number of valid Thetas instead of a fixed number of attempted Thetas
#! Not effective yet
sim_one_Theta_tilde = function(Theta_hat, cov_hat){
    success = FALSE
    num_attempts = 0
    while(!success){
        num_attempts = num_attempts + 1
        print(paste0("Number of Attempts: ", num_attempts))
        tryCatch({
            this_Theta_tilde = MASS::mvrnorm(1, mu = Theta_hat, Sigma = cov_hat)
            check_theta(this_Theta_tilde)
            success = TRUE
        }, error = function(e){
            # print(e)
        })
    }
    return(this_Theta_tilde)
}

#' Simulate parameter estimates for Monte Carlo delta-method
#'
#' @param B Number of Monte Carlo samples to generate
#' @param Theta_hat Parameter estimate from data analysis. Used as mean for Monte Carlo sample.
#' @param cov_hat Covariance estimate from data analysis. Used as covariance for Monte Carlo sample.
#'
#' @details
#' Simulate from the fitted sampling distribution of $$\hat{\theta}$$. These samples are then used to approximate the sampling distribution of some function of $$\hat{\theta}$$.
#'
#'
#' @return An iid sample of Theta values. Number of rows is B, number of cols is length of Theta_hat.
#' @export
sim_Theta_tildes = function(B, Theta_hat, cov_hat){
    MASS::mvrnorm(B, mu = Theta_hat, Sigma = cov_hat)
    
}

#' Compute a sample of mediation effects from a sample of simulated parameter estimates
#'
#' @param scale Scale(s) on which to construct mediation effects.
#' @param w A vector of values for the model confounders.
#' @param some_Theta_tildes Sample of parameter estimates. Rows index sample, columns index parameter.
#' @param which_REs Which random effects have been included in the model.
#'
#' @return A sample of mediation effects, with rows indexing sample and columns indexing effect.
#' @export
Theta_tildes_2_MEs = function(scale = c("diff", "rat", "OR"), w, some_Theta_tildes, which_REs){
    some_ME_tildes = data.frame()

    for(j in seq_len(nrow(some_Theta_tildes))){

        tryCatch({
            # print(j)

            # if(j %% 1000 == 0) cat("j=", j, " of ", nrow(some_Theta_tildes), "\n", sep="")



            this_Theta_tilde = some_Theta_tildes[j,]

            # The number of parameters in theta_Y based on the number of REs for Y
            len_theta_Y = which_REs %>%
                            num_Y_REs() %>%
                            num_REs2theta_length()


            # Extract parameters
            b_Y = this_Theta_tilde[1:5]
            theta_Y = this_Theta_tilde[6:(5 + len_theta_Y)]
            b_M = this_Theta_tilde[(6 + len_theta_Y):(9 + len_theta_Y)]
            theta_M = this_Theta_tilde[(10 + len_theta_Y):length(this_Theta_tilde)]

            check_theta(theta_Y)
            check_theta(theta_M)


            this_ME_tilde = all_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M, which_REs=which_REs)
            some_ME_tildes = rbind(some_ME_tildes, this_ME_tilde)
        }, error = function(e){
            # print(e)
        })

    }

      ME_names = get_ME_names(scale)
      colnames(some_ME_tildes) = ME_names

    return(some_ME_tildes)
}


MC_delta = function(B, Theta_hat, cov_hat, scale = c("diff", "rat", "OR"), w, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){

    some_Theta_tildes = sim_Theta_tildes(B, Theta_hat, cov_hat)

    some_MEs = Theta_tildes_2_MEs(scale, w, some_Theta_tildes, which_REs)

    this_cov = cov(some_MEs)

    return(this_cov)

}





# ---------------------------------------------------------------------------- #
#                      Simulate Using TMB Parameterization                     #
# ---------------------------------------------------------------------------- #


get_num_vars_TMB <- function(fit){
    p = length(fixef(fit)$cond)

    num_pars = length(fit$fit$par)
    num_RE_pars = num_pars - p
    num_SDs = (sqrt(1 + 8*num_RE_pars) - 1) / 2
    q = num_SDs

    return(c(p, q))
}

# Split a vector of TMB pars into its consituents: FE, log-SDs and correlation Cholesky factors. Output a length-3 list
# p: number of fixed effects covariates
# q: number of random effects covariates
TMB_par_vec_2_list <- function(par_vec, p, q){
    FE_pars = par_vec[1:p]
    SD_pars = par_vec[(p + 1):(p + q)]
    corr_pars = par_vec[(p + q + 1):length(par_vec)]
    
    return(list(TMB_FE_pars = FE_pars, TMB_SD_pars = SD_pars, TMB_corr_pars = corr_pars))
}

sim_TMB_Theta_tildes = function(B, fit_Y, fit_M){
    # Number of fixed and random effects in each model. Used later to subdivide the parameter vector
    num_vars_Y = get_num_vars_TMB(fit_Y)
    p_Y = num_vars_Y[1]
    q_Y = num_vars_Y[2]

    num_vars_M = get_num_vars_TMB(fit_M)
    p_M = num_vars_M[1]
    q_M = num_vars_M[2]


    # Extract fitted parameters and covariances
    TMB_pars_Y = TMB_pars_list(fit_Y) %>% unlist
    TMB_pars_M = TMB_pars_list(fit_M) %>% unlist

    TMB_cov_mat_Y = glmmTMB_SE_mat(fit_Y)
    TMB_cov_mat_M = glmmTMB_SE_mat(fit_M)


    # Simulate new parameter values
    some_TMB_pars_Y_tilde = MASS::mvrnorm(B, mu = TMB_pars_Y, Sigma = TMB_cov_mat_Y)
    some_TMB_pars_M_tilde = MASS::mvrnorm(B, mu = TMB_pars_M, Sigma = TMB_cov_mat_M)


    # Split up simulated parameter vectors into FEs, log-SDs, and correlation Cholesky factors
    some_par_lists_Y = apply(some_TMB_pars_Y_tilde, 1, TMB_par_vec_2_list, p=p_Y, q=q_Y)
    some_par_lists_M = apply(some_TMB_pars_M_tilde, 1, TMB_par_vec_2_list, p=p_M, q=q_M)


    # Convert TMB parameterization to our parameterization
    GLMM_pars_Y = t(sapply(some_par_lists_Y, function(this_TMB_par_list){
        this_TMB_par_list %>%
            TMB_2_GLMM_pars_list %>%    # Convert from TMB parameters to GLMM parameters
            organize_GLMM_par_list      # Re-organize REs to match our parameterization
    }))

    GLMM_pars_M = t(sapply(some_par_lists_M, function(this_TMB_par_list){
        this_TMB_par_list %>%
            TMB_2_GLMM_pars_list %>%    # Convert from TMB parameters to GLMM parameters
            organize_GLMM_par_list      # Re-organize REs to match our parameterization
    }))


    output = cbind(GLMM_pars_Y, GLMM_pars_M)

    return(output)

}



MC_delta_TMB = function(B, fit_Y, fit_M, scale = c("diff", "rat", "OR"), w, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")){

    some_TMB_Theta_tildes = sim_TMB_Theta_tildes(B, fit_Y, fit_M)
    some_Theta_tildes = t(sapply(seq_len(B), function(i) TMB_Theta_2_Theta(some_TMB_Theta_tildes[i,]) ))

    some_MEs = Theta_tildes_2_MEs(scale, w, some_Theta_tildes, which_REs)

    this_cov = cov(some_MEs)

    return(this_cov)

}