


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
source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
source("R/Exact_Asymptotics/Imai Method.r")
devtools::load_all()



load("R/Paper MC Study/all_datasets.RData", verbose = TRUE)


num_reps = length(all_datasets)

# all_ME_hats = list()
# all_cov_hats_delta = list()
# all_cov_hats_MC_delta = list()

# total_runtime_delta = 0
# total_runtime_MC_delta = 0

# #? Note: Extracting the SE matrix for the reg pars is slow (using merDeriv::vcov.glmerMod()). I don't want to duplicate this step (or fitting the models), so I separated out model fitting/SE extraction from my method and the MC delta.


# tic()
# for(i in seq_len(num_reps)){

#     print(paste0("i = ", i, " of ", num_reps))

#     data = all_datasets[[i]]


#     tryCatch({
            
#         ## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
#         (fit_Y = suppressMessages(lme4::glmer(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))
#         (fit_M = suppressMessages(lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))

        
#         # Compute parameter estimates
#         this_MEs = all_MEs_models(w = w, fit_Y = fit_Y, fit_M = fit_M)

#         ## Extract model parameter estimates
#         info_Y = get_model_pars(fit_Y)
#         info_M = get_model_pars(fit_M)
#         this_Theta_hat = c(unlist(info_Y), unlist(info_M))
#         this_cov_hat = all_pars_cov_mat(fit_Y, fit_M)

#         # Delta Method covariance
#         time = Sys.time()

#         this_delta_cov = all_covs_MEs_models(w = w, Sigma = this_cov_hat, fit_Y = fit_Y, fit_M = fit_M)

#         this_delta_runtime = Sys.time() - time
#         total_runtime_delta = total_runtime_delta + this_delta_runtime



#         # MC delta method covariance
#         time = Sys.time()

#         some_Theta_tildes = sim_Theta_tildes(B, this_Theta_hat, this_cov_hat)
#         some_ME_tildes = Theta_tildes_2_MEs(scale, w, some_Theta_tildes, which_REs)
#         this_MC_delta_cov = cov(some_ME_tildes)

#         this_MC_delta_runtime = Sys.time() - time
#         total_runtime_MC_delta = total_runtime_MC_delta + this_MC_delta_runtime

#         }, error = function(e){
#             this_MEs = NULL
#             this_ME_cov = NULL
#     })


#     all_ME_hats[[i]] = this_MEs
#     all_cov_hats_delta[[i]] = this_delta_cov
#     all_cov_hats_MC_delta[[i]] = this_MC_delta_cov


# }

# time_info = toc()
# runtime_total = time_info$toc - time_info$tic

# save(all_ME_hats, all_cov_hats_delta, all_cov_hats_MC_delta, file = "R/Paper MC Study/Paper - Delta and MC Delta Results.RData")
# save(runtime_total, total_runtime_delta, total_runtime_MC_delta, file = "R/Paper MC Study/Paper - Delta and MC Delta Runtimes.RData")


load("R/Paper MC Study/Paper - Delta and MC Delta Results.RData", verbose = TRUE)
load("R/Paper MC Study/Paper - Delta and MC Delta Runtimes.RData", verbose = TRUE)

#! Delta method has one pathological entry. Let's remove it and proceed
delta_traces = sapply(all_cov_hats_delta, function(x) sum(diag(x)))
which_err = which(delta_traces > 10000)
all_ME_hats = all_ME_hats[-which_err]
all_cov_hats_delta = all_cov_hats_delta[-which_err]
all_cov_hats_MC_delta = all_cov_hats_MC_delta[-which_err]



# Empirical covariance matrix of parameter estimates

## Convert all_ME_hats from a list to a data frame
all_ME_hats_data = do.call(rbind, all_ME_hats)

## Compute covariance
emp_cov = cov(all_ME_hats_data)
diag(emp_cov)



# Summaries of estimated SEs (mean of whole matrix and SEs of diagonal entries)

## Elementwise means
mean_delta_SEs = Reduce("+", all_cov_hats_delta) / length(all_cov_hats_delta)
mean_MC_delta_SEs = Reduce("+", all_cov_hats_MC_delta) / length(all_cov_hats_MC_delta)

diag(mean_delta_SEs)

## SEs of diagonal entries
### Compute the SE of the diagonal entries' means from a list of matrices
diag_mean_SEs = function(some_matrices){
    some_matrices %>% lapply(diag) %>% do.call(rbind, .) %>% apply(2, sd) / length(all_cov_hats_delta)
}

delta_diag_SEs = diag_mean_SEs(all_cov_hats_delta)
MC_delta_diag_SEs = diag_mean_SEs(all_cov_hats_MC_delta)



#* MSEs for delta method
all_delta_errors = sapply(all_cov_hats_delta, function(x) diag(x) - diag(emp_cov))
all_sq_delta_errors = all_errors^2
delta_MSEs = apply(all_sq_delta_errors, 1, mean)


#* Coverage rates for delta method
get_CIs = function(ME_hats, SEs){
    lower = ME_hats - 1.96 * SEs
    upper = ME_hats + 1.96 * SEs
    return(list(lcl = lower, ucl = upper))
}

delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_delta, ~ get_CIs(.x, diag(.y)))


load("R/Paper MC Study/true_MEs.RData", verbose = TRUE)

delta_CI_checks = lapply(delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)

delta_CI_rates = Reduce("+", delta_CI_checks) / length(delta_CI_checks)




#* Repeat for MC delta method

##* MSEs for MC delta method
all_MC_delta_errors = sapply(all_cov_hats_MC_delta, function(x) diag(x) - diag(emp_cov))
all_sq_MC_delta_errors = all_MC_delta_errors^2
MC_delta_MSEs = apply(all_sq_MC_delta_errors, 1, mean)

##* Coverage rates for MC delta method
MC_delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_MC_delta, ~ get_CIs(.x, diag(.y)))
MC_delta_CI_checks = lapply(MC_delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
MC_delta_CI_rates = Reduce("+", MC_delta_CI_checks) / length(MC_delta_CI_checks)