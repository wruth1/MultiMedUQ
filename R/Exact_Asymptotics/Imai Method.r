






# ---------------------------------------------------------------------------- #
#                       Functions to do MC-delta for MEs                       #
# ---------------------------------------------------------------------------- #


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
