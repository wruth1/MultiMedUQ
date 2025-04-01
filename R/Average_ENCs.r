

#? This script contains functions for averaging ENCs over the confounder distribution



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

    #TODO: Increase n_samples (was 5000)

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


#' @rdname Marginal_ENCs
#' @export
mean_ENC_Theta <- function(Theta, all_confounders, confounder_probs, len_par_vecs, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X"), fast = FALSE, n_samples = 500) {

    #TODO: Increase n_samples (was 5000)

    num_confounder_vals = length(all_confounders)

    if  (!fast || (num_confounder_vals <= n_samples)) {     # Compute expected value over confounder distribution exactly
        some_confounder_inds = seq_along(all_confounders)
        sampled = FALSE
    } else {                                                # Compute expected value over confounder distribution by sampling
        some_confounder_inds = sample(num_confounder_vals, size = n_samples, replace = FALSE, prob = confounder_probs)
        sampled = TRUE
    }

    diag_terms = matrix(0, nrow = 4, ncol = 4)

    for (j in some_confounder_inds){
        this_confounder_val = all_confounders[[j]]

        this_cov_ENCs_delta = all_covs_ENC_Theta(this_confounder_val, Sigma, Theta, len_par_vecs, which_REs =  which_REs)

        this_prob = confounder_probs[j]

        if (sampled) {
            diag_terms = diag_terms + this_cov_ENCs_delta * this_prob / n_samples   # Note: Each summand has a squared probability. That's why we've got this formula that seems wrong.
        } else {
            diag_terms = diag_terms + confounder_probs[j]^2 * this_cov_ENCs_delta
        }
    }

    num_confounder_vals = length(all_confounders)

    mean_ENC_hat = rep(0, times = 4)

    if (!fast || (num_confounder_vals <= n_samples)) { # Compute expected value over confounder distribution exactly

        for (j in seq_along(all_confounders)){
            if (j %% 100 == 0) print(j)
            this_confounder_val = all_confounders[[j]]

            this_ENCs = all_ENCs_Theta(this_confounder_val, Theta, which_REs =  which_REs, len_par_vecs = len_par_vecs)

            mean_ENC_hat = mean_ENC_hat + confounder_probs[j] * this_ENCs
        }
    } else { # Compute expected value over confounder distribution by sampling
        some_confounder_inds = sample(num_confounder_vals, size = n_samples, replace = FALSE, prob = confounder_probs)

        for (j in some_confounder_inds){
            this_confounder_val = all_confounders[[j]]

            this_ENCs = all_ENCs_Theta(this_confounder_val, Theta, which_REs =  which_REs, len_par_vecs = len_par_vecs)

            mean_ENC_hat = mean_ENC_hat + this_ENCs / n_samples
        }
    }

    return(mean_ENC_hat)
}




# ---------------------------------------------------------------------------- #
#                        Sampling Covariance of Mean ENC                       #
# ---------------------------------------------------------------------------- #

marginal_cov_hat_diag_terms <- function(Sigma, b_Y, theta_Y, b_M, theta_M, all_confounders, confounder_probs, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X"), fast = FALSE, n_samples = 500) {

    #TODO: Increase n_samples (was 5000)

    num_confounder_vals = length(all_confounders)

    if  (!fast || (num_confounder_vals <= n_samples)) {     # Compute expected value over confounder distribution exactly
        some_confounder_inds = seq_along(all_confounders)
        sampled = FALSE
    } else {                                                # Compute expected value over confounder distribution by sampling
        some_confounder_inds = sample(num_confounder_vals, size = n_samples, replace = FALSE, prob = confounder_probs)
        sampled = TRUE
    }

    diag_terms = matrix(0, nrow = 4, ncol = 4)

    for (j in some_confounder_inds){
        this_confounder_val = all_confounders[[j]]

        this_cov_ENCs_delta = all_covs_ENC_pars(this_confounder_val, Sigma, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

        this_prob = confounder_probs[j]

        if (sampled) {
            diag_terms = diag_terms + this_cov_ENCs_delta * this_prob / n_samples   # Note: Each summand has a squared probability. That's why we've got this formula that seems wrong.
        } else {
            diag_terms = diag_terms + confounder_probs[j]^2 * this_cov_ENCs_delta
        }
    }

    return(diag_terms)
}

marginal_cov_hat_diag_terms_Theta <- function(Sigma, Theta, all_confounders, confounder_probs, len_par_vecs, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X"), fast = FALSE, n_samples = 500) {

    len_b_Y = len_par_vecs[1]
    len_theta_Y = len_par_vecs[2]
    len_b_M = len_par_vecs[3]
    len_theta_M = len_par_vecs[4]

    this_b_Y = Theta[1:len_b_Y]
    this_theta_Y = Theta[(len_b_Y + 1):(len_b_Y + len_theta_Y)]
    this_b_M = Theta[(len_b_Y + len_theta_Y + 1):(len_b_Y + len_theta_Y + len_b_M)]
    this_theta_M = Theta[(len_b_Y + len_theta_Y + len_b_M + 1):(len_b_Y + len_theta_Y + len_b_M + len_theta_M)]

    #TODO: Increase n_samples (was 5000)
    
    num_confounder_vals = length(all_confounders)

    if  (!fast || (num_confounder_vals <= n_samples)) {     # Compute expected value over confounder distribution exactly
        some_confounder_inds = seq_along(all_confounders)
        sampled = FALSE
    } else {                                                # Compute expected value over confounder distribution by sampling
        some_confounder_inds = sample(num_confounder_vals, size = n_samples, replace = FALSE, prob = confounder_probs)
        sampled = TRUE
    }

    diag_terms = matrix(0, nrow = 4, ncol = 4)

    for (j in some_confounder_inds){
        this_confounder_val = all_confounders[[j]]

        this_cov_ENCs_delta = all_covs_ENC_Theta(this_confounder_val, Sigma, Theta, len_par_vecs, which_REs =  which_REs)

        this_prob = confounder_probs[j]

        if (sampled) {
            diag_terms = diag_terms + this_cov_ENCs_delta * this_prob / n_samples   # Note: Each summand has a squared probability. That's why we've got this formula that seems wrong.
        } else {
            diag_terms = diag_terms + confounder_probs[j]^2 * this_cov_ENCs_delta
        }
    }

    return(diag_terms)
}

