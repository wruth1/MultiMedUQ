

#? This script contains functions for averaging ENCs over the confounder distribution

#TODO: Increase all instances of n_samples (was 5000)




# ---------------------------------------------------------------------------- #
#                                   Utilities                                  #
# ---------------------------------------------------------------------------- #

#' Extract the distribution of confounders from a dataset
#'
#' @param data A data frame containing mediation variables and confounders
#' @param outcome_name,exposure_name,mediator_name,group_name Names of the outcome, exposure, mediator, and group variables respectively. Must be strings.
#'
#' @returns A list with two components: `values`, a list of vectors of confounder values, and `freqs`, a vector of frequencies of each confounder value.
#' @export
data_2_confounder_dist <- function(data, outcome_name, exposure_name, mediator_name, group_name) {
    data_confounders = data %>% dplyr::select(-dplyr::all_of(c(outcome_name, exposure_name, mediator_name, group_name)))

    info_confounders = data_confounders %>% dplyr::group_by(dplyr::across(dplyr::everything())) %>% dplyr::summarize(.freq = dplyr::n(), .groups = "drop")

    vals_confounders = info_confounders %>% dplyr::select(-.freq) %>% split(seq_len(nrow(.))) %>% lapply(unlist)
    freqs_confounders = info_confounders$.freq

    output = list(values = vals_confounders, freqs = freqs_confounders)
    return(output)
}




# ---------------------------------------------------------------------------- #
#                             Estimator of Mean ENC                            #
# ---------------------------------------------------------------------------- #

#' Marginal ENC
#'
#' Compute the marginal ENC by averaging conditional ENCs over the provided confounder distribution
#'
#' @param b_Y,b_M Coefficient vectors for \eqn{Y}-model and \eqn{M}-model, respectively.
#' @param theta_Y,theta_M Covariance parameters of random effects in \eqn{Y}-model and \eqn{M}-model, respectively. See details.
#' @param Theta Vector of parameters from both models. Order is b_Y, theta_Y, b_M, theta_M
#' @param all_confounders A list of all possible values of the confounder variables. Each element is a vector.
#' @param confounder_probs A vector of probabilities, with element `i` corresponding to the probability of the `i`th element of `all_confounders`.
#' @param which_REs Which random effects to include in the calculation. Default is all. Shorthands are available. See details.
#' @param len_par_vecs Number of entries in each parameter vector. Order is b_Y, theta_Y, b_M, theta_M
#' @param fast Should the mean be approximated by Monte Carlo?
#' @param n_samples Number of samples to use in the Monte Carlo approximation
#'
#' @name Marginal_ENCs
#'
#' @details
#' The following shorthands for random effects are available:
#' \itemize{
#' \item "All": All REs
#' \item "Y.All": All REs for Y
#' \item "M.All": All REs for M
#' }
#' Additionally, individual REs can be specified:
#' \itemize{
#' \item "Y.Int": Intercept for Y
#' \item "Y.X": Slope for X in Y
#' \item "Y.M": Slope for M in Y
#' \item "M.Int": Intercept for M
#' \item "M.X": Slope for M
#' }
#'
#'
#' @returns A vector of marginal expected nested counterfactuals (i.e. averaged over the confounder distribution). Order of output is \code{ENC(1,1), ENC(1,0), ENC(0,1), ENC(0,0)}, where the former argument is \eqn{X} and the latter is \eqn{X_M}.
#' @export
mean_ENC_pars <- function(b_Y, theta_Y, b_M, theta_M, all_confounders, confounder_probs, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X"), fast = FALSE, n_samples = 500) {

    num_confounder_vals = length(all_confounders)

    mean_ENC_hat = rep(0, times = 4)

    if (!fast || (num_confounder_vals <= n_samples)) { # Compute expected value over confounder distribution exactly

        for (j in seq_along(all_confounders)){
            this_confounder_val = all_confounders[[j]]

            this_ENCs = all_ENCs(this_confounder_val, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

            mean_ENC_hat = mean_ENC_hat + confounder_probs[j] * this_ENCs
        }
    } else { # Compute expected value over confounder distribution by sampling

        some_confounder_inds = sample(num_confounder_vals, size = n_samples, replace = FALSE, prob = confounder_probs)

        for (j in some_confounder_inds) {
            this_confounder_val = all_confounders[[j]]

            this_ENCs = all_ENCs(this_confounder_val, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

            mean_ENC_hat = mean_ENC_hat + this_ENCs / n_samples
        }
    }

    return(mean_ENC_hat)
}

#? I chose to have this function depend on the other to avoid duplicating code. Specifically, to avoid having to remember to make changes to both functions.
#' @rdname Marginal_ENCs
#' @export
mean_ENC_Theta <- function(Theta, all_confounders, confounder_probs, len_par_vecs, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X"), fast = FALSE, n_samples = 500) {

    len_b_Y = len_par_vecs[1]
    len_theta_Y = len_par_vecs[2]
    len_b_M = len_par_vecs[3]
    len_theta_M = len_par_vecs[4]

    this_b_Y = Theta[1:len_b_Y]
    this_theta_Y = Theta[(len_b_Y + 1):(len_b_Y + len_theta_Y)]
    this_b_M = Theta[(len_b_Y + len_theta_Y + 1):(len_b_Y + len_theta_Y + len_b_M)]
    this_theta_M = Theta[(len_b_Y + len_theta_Y + len_b_M + 1):(len_b_Y + len_theta_Y + len_b_M + len_theta_M)]



    return(mean_ENC_pars(this_b_Y, this_theta_Y, this_b_M, this_theta_M, all_confounders, confounder_probs, which_REs = which_REs, fast = fast, n_samples = n_samples))
}




# ---------------------------------------------------------------------------- #
#                        Sampling Covariance of Mean ENC                       #
# ---------------------------------------------------------------------------- #


# ------------------------- Variance (diagonal) terms ------------------------ #

marginal_cov_hat_var_terms_par <- function(Sigma, b_Y, theta_Y, b_M, theta_M, all_confounders, confounder_probs, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X"), fast = FALSE, n_samples = 500) {

    num_confounder_vals = length(all_confounders)

    if  (!fast || (num_confounder_vals <= n_samples)) {     # Compute expected value over confounder distribution exactly
        some_confounder_inds = seq_along(all_confounders)
        sampled = FALSE
    } else {                                                # Compute expected value over confounder distribution by sampling
        some_confounder_inds = sample(num_confounder_vals, size = n_samples, replace = FALSE, prob = confounder_probs)
        sampled = TRUE
    }

    var_terms = matrix(0, nrow = 4, ncol = 4)

    for (j in some_confounder_inds){
        this_confounder_val = all_confounders[[j]]

        this_cov_ENCs_delta = all_covs_ENC_pars(this_confounder_val, Sigma, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

        this_prob = confounder_probs[j]

        if (sampled) {
            var_terms = var_terms + this_cov_ENCs_delta * this_prob / n_samples   # Note: Each summand has a squared probability. That's why we've got this formula that seems wrong.
        } else {
            var_terms = var_terms + confounder_probs[j]^2 * this_cov_ENCs_delta
        }
    }

    return(var_terms)
}

#? As above,I chose to have this function depend on the other to avoid duplicating code. Specifically, to avoid having to remember to make changes to both functions.
marginal_cov_hat_var_terms_Theta <- function(Sigma, Theta, all_confounders, confounder_probs, len_par_vecs, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X"), fast = FALSE, n_samples = 500) {

    len_b_Y = len_par_vecs[1]
    len_theta_Y = len_par_vecs[2]
    len_b_M = len_par_vecs[3]
    len_theta_M = len_par_vecs[4]

    this_b_Y = Theta[1:len_b_Y]
    this_theta_Y = Theta[(len_b_Y + 1):(len_b_Y + len_theta_Y)]
    this_b_M = Theta[(len_b_Y + len_theta_Y + 1):(len_b_Y + len_theta_Y + len_b_M)]
    this_theta_M = Theta[(len_b_Y + len_theta_Y + len_b_M + 1):(len_b_Y + len_theta_Y + len_b_M + len_theta_M)]

    return(marginal_cov_hat_var_terms_par(Sigma, this_b_Y, this_theta_Y, this_b_M, this_theta_M, all_confounders, confounder_probs, which_REs =  which_REs, fast = fast, n_samples = n_samples))
}



# ---------------------- Covariance (off-diagonal) terms --------------------- #

marginal_cov_hat_cov_terms_par <- function(Sigma, b_Y, theta_Y, b_M, theta_M, all_confounders, confounder_probs, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X"), fast = FALSE, n_samples = 500) {

    # all_confounders_backup = all_confounders
    # all_confounders = all_confounders[1:200]

    all_confounder_inds = seq_along(all_confounders)
    all_confounder_ind_pairs = combn(all_confounder_inds, 2, simplify = FALSE)

    num_confounder_vals = length(all_confounders)
    num_confounder_pairs = num_confounder_vals * (num_confounder_vals - 1) / 2

    # Probability of each pair of confounder values
    confounder_pair_probs = sapply(all_confounder_ind_pairs, function(x) confounder_probs[x[1]] * confounder_probs[x[2]])


    if  (!fast || (num_confounder_pairs <= n_samples)) {     # Compute expected value over confounder distribution exactly
        some_inds_for_confounder_pairs = seq_along(all_confounder_ind_pairs)
        sampled = FALSE
    } else {                                                # Compute expected value over confounder distribution by sampling
        some_inds_for_confounder_pairs = sample(num_confounder_pairs, size = n_samples, replace = FALSE, prob = confounder_pair_probs)
        sampled = TRUE
    }

    cov_terms = matrix(0, nrow = 4, ncol = 4)

    for (j in some_inds_for_confounder_pairs){
        this_confounder_ind_pair = all_confounder_ind_pairs[[j]]

        this_confounder_val_1 = all_confounders[[this_confounder_ind_pair[1]]]
        this_confounder_val_2 = all_confounders[[this_confounder_ind_pair[2]]]

        this_cross_cov = cross_cov_ENC_pars(this_confounder_val_1, this_confounder_val_2, Sigma, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)


        if (sampled) {
            cov_terms = cov_terms + 2 * this_cross_cov / n_samples                    # Monte Carlo summation -> average
        } else {
            cov_terms = cov_terms + 2 * this_cross_cov * confounder_pair_probs[j]     # Raw summation -> multiply by prob, then add (not average...technically, a weighted average)
        }
    }

    return(cov_terms)

}


marginal_cov_hat_cov_terms_Theta <- function(Sigma, Theta, all_confounders, confounder_probs, len_par_vecs, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X"), fast = FALSE, n_samples = 500) {

    len_b_Y = len_par_vecs[1]
    len_theta_Y = len_par_vecs[2]
    len_b_M = len_par_vecs[3]
    len_theta_M = len_par_vecs[4]

    this_b_Y = Theta[1:len_b_Y]
    this_theta_Y = Theta[(len_b_Y + 1):(len_b_Y + len_theta_Y)]
    this_b_M = Theta[(len_b_Y + len_theta_Y + 1):(len_b_Y + len_theta_Y + len_b_M)]
    this_theta_M = Theta[(len_b_Y + len_theta_Y + len_b_M + 1):(len_b_Y + len_theta_Y + len_b_M + len_theta_M)]

    return(marginal_cov_hat_cov_terms_par(Sigma, this_b_Y, this_theta_Y, this_b_M, this_theta_M, all_confounders, confounder_probs, which_REs =  which_REs, fast = fast, n_samples = n_samples))
}
