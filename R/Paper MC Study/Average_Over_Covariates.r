



library(lme4)
library(merDeriv)
library(tictoc)
library(pbapply)
library(parallel)
library(magrittr)
library(dplyr)
library(kableExtra)
library(ggplot2)
library(ggmulti)
library(broom.mixed)
library(glmmTMB)
library(stringr)
library(tibble)
source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
source("R/Exact_Asymptotics/Imai Method.r")
devtools::load_all()




# ---------------------------------------------------------------------------- #
#                             Initialize Parameters                            #
# ---------------------------------------------------------------------------- #


# ----------------------------------- Note ----------------------------------- #

#? I denote regression coefficients for Y as b_Y and for M as b_M
#? Random effect parameters for Y are theta_Y and for M as theta_M
    #? These RE parameters are SDs and correlations. They are organized as, e.g. with 3 REs, SD_1, corr_12, corr_13, SD_2, corr_23, SD_3
#? We denote the single vector containing all parameters by Theta = (c(b_Y, theta_Y, b_M, theta_M))
    #? I.e. Capital Theta contains both small theta vectors




# Set parameters

B = 500     # Number of samples to generate for MC delta
scale = c("diff", "rat", "OR")      # What scales should we compute mediation effects on?
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")        # Which variables have random effects? Eventually, I will need a better way to specify this

# Number of groups
K = 100
# K=10

# Observations per group
# N = 500
N=1000
# N = 10000
n = N






# Vector of covariates
## Note that mediation effects currently require us to specify values for the confounders. I think the next move is to average these conditional effects over the observed distribution of confounders.
w = c(2,3)



## Non-trivial values for the b's and theta's. Former based on output from another MC study. Latter chosen arbitrarily.
## Crucially, no parameters are equal to zero.
##? We choose the intercepts to that the mean of the linear predictor is zero. Doing this for M makes it easier to do so for Y.
#! Scale factor for coefficients and SDs
scale_factor = 1

b_Y_X = 0.966486302988689
b_Y_M = 1.99644760563721
b_Y_C1 = -1
b_Y_C2 = 1
b_Y_int = - sum(b_Y_X, b_Y_M, b_Y_C1, b_Y_C2) / 2       # ~ -1.48
b_Y = c(b_Y_int, b_Y_X, b_Y_M, b_Y_C1, b_Y_C2) * scale_factor


b_M_X = 1.76353928991247
b_M_C1 = 1
b_M_C2 = -1
b_M_int = -sum(b_M_X, b_M_C1, b_M_C2) / 2       # ~ -0.89
b_M = c(b_M_int, b_M_X, b_M_C1, b_M_C2) * scale_factor




# Choose theta_Y and theta_M based on the values of b_Y and b_M
theta_Y = c(scale_factor*sqrt(0.5), 0.3, 0.4, scale_factor, 0.5, scale_factor*sqrt(0.8)) / 3
theta_M = c(scale_factor*sqrt(0.5), -0.5, scale_factor) / 3


all_reg_pars = c(b_Y, theta_Y, b_M, theta_M)




p_Y = length(b_Y)
p_M = length(b_M)
p = p_Y + p_M


folder_suffix = paste0("K=", K, ", N=", N)
dir.create(paste0("R/Paper MC Study/Data - ", folder_suffix), showWarnings = F)
dir.create(paste0("R/Paper MC Study/Results - ", folder_suffix), showWarnings = F)



outcome_name = "Y"
exposure_name = "X"
mediator_name = "M"
group_name = "group"


# Takes a data frame. Extracts only the confounders. 
# Returns a list with two entries. First is a list of all unique values of the confounders. Second is a list of the frequencies of each unique value of the confounders.
get_confounder_freq <- function(data, outcome_name, exposure_name, mediator_name, group_name){
    data_confounders = data %>% dplyr::select(-all_of(c(outcome_name, exposure_name, mediator_name, group_name)))

    info_confounders = data_confounders %>% group_by(across(everything())) %>% summarize(.freq = n(), .groups = "drop")

    vals_confounders = info_confounders %>% dplyr::select(-.freq) %>% split(seq_len(nrow(.))) %>% lapply(unlist)
    freqs_confounders = info_confounders$.freq

    output = list(values = vals_confounders, freqs = freqs_confounders)
    return(output)
}





# cl = makeCluster(detectCores() - 2)
# cl = makeCluster(15)
cl = makeCluster(10)
# clusterExport(cl, c("N", "b_Y", "theta_Y", "b_M", "theta_M", "which_REs"))
clusterExport(cl, c("w", "B", "scale", "which_REs", "N", "n", "K", "b_Y", "theta_Y", "b_M", "theta_M", "folder_suffix", "get_confounder_freq", "outcome_name", "exposure_name", "mediator_name", "group_name"))
clusterEvalQ(cl, {
    library(lme4)
    library(merDeriv)
    library(tictoc)
    library(pbapply)
    library(parallel)
    library(magrittr)
    library(dplyr)
    library(kableExtra)
    library(ggplot2)
    library(ggmulti)
    library(broom.mixed)
    library(glmmTMB)
    source("R/For_Bouchra/Exact_Asymptotics_Helpers.r")
    source("R/For_Bouchra/Imai Method.r")
    devtools::load_all()
})
clusterSetRNGStream(cl = cl, 123)
# clusterSetRNGStream(cl = cl, 11111111)



# Fit models, extract MEs, estimate covariance matrices and save results
#! Before running, un-comment cl=cl and both save()'s
MC_results_delta_MC_delta = pblapply(1:num_datasets, function(i) {
# MC_results_delta_MC_delta = pblapply(1:3, function(i) {
    load(paste0("R/Paper MC Study/Data - ", folder_suffix, "/", i, ".RData"), verbose = T)


    tryCatch({

        this_timings = list()


        # ---------------------------- Delta Method (ours) --------------------------- #

        # Fit models
        tic()

        fit_Y = glmmTMB(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial) #, control = glmmTMBControl(optimizer = "optim", optArgs = list(method = "BFGS", eval.max = 1e10)))
        fit_M = glmmTMB(M ~ X + C1 + C2 + (X | group), data = data, family = binomial) #, control = glmmTMBControl(optimizer = "optim", optArgs = list(method = "BFGS", eval.max = 1e8)))

        this_time = toc()
        this_timings$fit_models = this_time$toc - this_time$tic

        # diagnose(fit_Y)
        # diagnose(fit_M)



        # Extract fitted parameters

        tic()

        theta_hat_Y = get_model_pars_TMB(fit_Y)
        theta_hat_M = get_model_pars_TMB(fit_M)
        Theta_hat = c(unlist(theta_hat_Y), unlist(theta_hat_M))
        cov_hat = all_pars_cov_mat_TMB(fit_Y, fit_M)

        # cbind(Theta_hat, (diag(cov_hat)))

        b_Y = theta_hat_Y[["b"]]
        theta_Y = theta_hat_Y[["theta"]]
        b_M = theta_hat_M[["b"]]
        theta_M = theta_hat_M[["theta"]]
        len_par_vecs = sapply(list(b_Y, theta_Y, b_M, theta_M), length)


        this_time = toc()
        this_timings$get_pars = this_time$toc - this_time$tic
        # data_est = data.frame(hat = Theta_hat, SE = sqrt(diag(cov_hat)))
        # rownames(data_est) = names(Theta_hat)



        # Compute mediation effects, averaging over the EDF of the confounders
        tic()

        info_confounders = get_confounder_freq(data, outcome_name, exposure_name, mediator_name, group_name)
        all_confounders = info_confounders$values
        freqs_confounders = info_confounders$freq

        mean_ME_hat = rep(0, times = 3 * length(scale))
        mean_ME_hat_sq = rep(0, times = 3 * length(scale))
        mean_cov_hat = matrix(0, nrow = 3 * length(scale), ncol = 3 * length(scale))

        ME_hat_weights = freqs_confounders / sum(freqs_confounders)     # Weights for averaging based on observed frequencies of confounders


        for(j in seq_along(all_confounders)){
            this_confounder_val = all_confounders[[j]]

            this_MEs = all_MEs_pars(scale, this_confounder_val, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)
            this_cov_MEs_delta = all_covs_MEs_pars(scale, this_confounder_val, cov_hat, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

            mean_ME_hat = mean_ME_hat + ME_hat_weights[j] * this_MEs
            mean_ME_hat_sq = mean_ME_hat_sq + ME_hat_weights[j] * this_MEs^2
        }


        this_time = toc()
        this_timings$get_MEs = this_time$toc - this_time$tic

        # ------------------------------ MC Delta Method ----------------------------- #
        tic()
        some_Theta_tildes = sim_Theta_tildes(B, Theta_hat, cov_hat)
        some_ME_tildes = Theta_tildes_2_MEs(scale, w, some_Theta_tildes, which_REs, len_par_vecs = len_par_vecs)
        cov_MEs_MC_delta = cov(some_ME_tildes)



        this_time = toc()
        this_timings$MC_delta = this_time$toc - this_time$tic


        # ------------------------ Compile and return results ------------------------ #
        output = list(this_MEs = MEs, cov_MEs_delta = cov_MEs_delta, cov_MEs_MC_delta = cov_MEs_MC_delta, this_timings = this_timings)

        # save(output, file = paste0("R/Paper MC Study/Results - ", folder_suffix, "/", i, ".RData"))
        return(output)
    }, error = function(e){
        output = NULL

        # save(output, file = paste0("R/Paper MC Study/Results - ", folder_suffix, "/", i, ".RData"))
        return(output)
    })


# }, cl = cl)
})


stopCluster(cl)














# ---------------------------------------------------------------------------- #
#                                     ENCs                                     #
# ---------------------------------------------------------------------------- #

clusterSetRNGStream(cl = cl, 123)


num_datasets = length(list.files(paste0("R/Paper MC Study/Data - ", folder_suffix)))
num_datasets = 10

# Fit models, extract ENCs, estimate covariance matrices and save results
#! Before running, un-comment cl=cl and both save()'s
MC_results_ENC = pblapply(1:num_datasets, function(i) {
# MC_results_delta_MC_delta = pblapply(1:3, function(i) {
    load(paste0("R/Paper MC Study/Data - ", folder_suffix, "/", i, ".RData"), verbose = T)




    tryCatch({

        # ---------------------------- Delta Method (ours) --------------------------- #

        # Fit models
        fit_Y = glmmTMB(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = glmmTMBControl(optimizer = "optim", optArgs = list(method = "BFGS", eval.max = 1e10)))
        fit_M = glmmTMB(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = glmmTMBControl(optimizer = "optim", optArgs = list(method = "BFGS", eval.max = 1e10)))



        # Extract fitted parameters


        theta_hat_Y = get_model_pars_TMB(fit_Y)
        theta_hat_M = get_model_pars_TMB(fit_M)
        Theta_hat = c(unlist(theta_hat_Y), unlist(theta_hat_M))
        cov_hat = all_pars_cov_mat_TMB(fit_Y, fit_M)

        # cbind(Theta_hat, (diag(cov_hat)))

        b_Y = theta_hat_Y[["b"]]
        theta_Y = theta_hat_Y[["theta"]]
        b_M = theta_hat_M[["b"]]
        theta_M = theta_hat_M[["theta"]]
        len_par_vecs = sapply(list(b_Y, theta_Y, b_M, theta_M), length)



        # --- Compute mediation effects, averaging over the EDF of the confounders --- #

        # Also, estimate sampling covariance matrix of (wt) average mediation effect, averaging over the EDF of the confounders appropriately using the conditional variance decomposition formula.

        info_confounders = get_confounder_freq(data, outcome_name, exposure_name, mediator_name, group_name)
        all_confounders = info_confounders$values
        freqs_confounders = info_confounders$freq

        mean_ENC_hat = rep(0, times = 4)

        ENC_hat_weights = freqs_confounders / sum(freqs_confounders)     # Weights for averaging based on observed frequencies of confounders


        for(j in seq_along(all_confounders)){
            this_confounder_val = all_confounders[[j]]

            this_ENCs = all_ENCs(this_confounder_val, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

            mean_ENC_hat = mean_ENC_hat + ENC_hat_weights[j] * this_ENCs
        }



        # ------------ Estimate sampling covariance of averaged estimator ------------ #

        cov_hat_mean_ENCs = matrix(0, nrow=4, ncol=4)


        #* Diagonal elements

        diag_terms = matrix(0, nrow=4, ncol=4)

        for(j in seq_along(all_confounders)){
            this_confounder_val = all_confounders[[j]]

            this_cov_ENCs_delta = all_covs_ENC_pars(this_confounder_val, cov_hat, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

            diag_terms = diag_terms + ENC_hat_weights[j]^2 * this_cov_ENCs_delta
        }

        cov_hat_mean_ENCs = cov_hat_mean_ENCs + diag_terms

        #* Off-diagonal elements

        off_diag_terms = matrix(0, nrow=4, ncol=4)

        num_confounders = length(all_confounders)
        for(j1 in 1:(num_confounders-1)){
            for(j2 in (j1+1):num_confounders){
                this_confounder_val1 = all_confounders[[j1]]
                this_confounder_val2 = all_confounders[[j2]]

                this_cross_cov = cross_cov_ENC_pars(this_confounder_val1, this_confounder_val2, cov_hat, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

                off_diag_terms = off_diag_terms + 2 * ENC_hat_weights[j1] * ENC_hat_weights[j2] * this_cross_cov
            }
        }

        cov_hat_mean_ENCs = cov_hat_mean_ENCs + off_diag_terms








        output = list(this_ENCs = mean_ENC_hat, cov_hat_mean_ENCs = cov_hat_mean_ENCs, diag_terms = diag_terms, off_diag_terms = off_diag_terms)

        # # ------------------------------ MC Delta Method ----------------------------- #
        # tic()
        # some_Theta_tildes = sim_Theta_tildes(B, Theta_hat, cov_hat)
        # some_ME_tildes = Theta_tildes_2_MEs(scale, w, some_Theta_tildes, which_REs, len_par_vecs = len_par_vecs)
        # cov_MEs_MC_delta = cov(some_ME_tildes)



        # this_time = toc()
        # this_timings$MC_delta = this_time$toc - this_time$tic


        # # ------------------------ Compile and return results ------------------------ #
        # output = list(this_MEs = MEs, cov_MEs_delta = cov_MEs_delta, cov_MEs_MC_delta = cov_MEs_MC_delta, this_timings = this_timings)

        # save(output, file = paste0("R/Paper MC Study/Results - ", folder_suffix, "/", i, ".RData"))
        return(output)
    }, error = function(e){
        output = NULL

        # save(output, file = paste0("R/Paper MC Study/Results - ", folder_suffix, "/", i, ".RData"))
        return(output)
    })


}, cl = cl)
# })



# ------------------------------ Extract results ----------------------------- #

#* Remove iterations with negative variances
MC_results_ENC_raw = MC_results_ENC
all_covs_raw = lapply(MC_results_ENC, function(x) x$cov_hat_mean_ENCs)
inds_remove = sapply(all_covs_raw, function(x) any(diag(x) < 0))
any(inds_remove)
MC_results_ENC = MC_results_ENC[!inds_remove]

cov_hat_ew_mins = sapply(all_covs_raw, function(x) min(Re(eigen(x)$values)))
hist(cov_hat_ew_mins[inds_remove])

all_ENCs = sapply(MC_results_ENC, function(x) x$this_ENCs) %>% t()
all_covs = lapply(MC_results_ENC, function(x) x$cov_hat_mean_ENCs)


emp_cov = cov(all_ENCs)
mean_cov_hat = Reduce("+", all_covs) / length(all_covs)
mean_cov_hat = (mean_cov_hat + t(mean_cov_hat)) / 2


norm(emp_cov - mean_cov_hat, "2") / min(norm(emp_cov, "2"), norm(mean_cov_hat, "2"))


all_cov_hat_errs = sapply(all_covs, function(x) mat_rel_err(emp_cov, x))
hist(all_cov_hat_errs)


mat_rel_err(emp_cov, mean_mean_cov)


