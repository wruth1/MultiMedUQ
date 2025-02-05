


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
source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
source("R/Exact_Asymptotics/Imai Method.r")
devtools::load_all()




# Set parameters

B = 500
scale = c("diff", "rat", "OR")
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")

# N = 100
# N = 20
# N = 40
# N = 60
N=1000

#* Main value
# N = 500
n = N


# all_Ks = c(50, 100, 200, 400, 800)
# all_Ks = c(50, 100, 200)
# all_Ks = 50 * (2:6)
# K=1000

#* Main value
# K = 100
K=10







# which_REs = c("Y.Int", "Y.X", "M.All")
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")



x = 0
x_m = 1



w = c(2,3)



## Non-trivial values for the b's and theta's. Former based on output from another MC study. Latter chosen arbitrarily.
## Crucially, no parameters are equal to zero.
##? We choose the intercepts to that the mean of the linear predictor is zero. Doing this for M makes it easier to do so for Y.
#! Scale factor for coefficients and SDs
scale_factor = 1
# b_Y_int_old = 0.0376828219852018
b_Y_X = 0.966486302988689
b_Y_M = 1.99644760563721
b_Y_C1 = -1
b_Y_C2 = 1
b_Y_int = - sum(b_Y_X, b_Y_M, b_Y_C1, b_Y_C2) / 2       # ~ -1.48
b_Y = c(b_Y_int, b_Y_X, b_Y_M, b_Y_C1, b_Y_C2) * scale_factor
# b_Y = c(0.0376828219852018, 0.966486302988689, 1.99644760563721, -0.00556557712859059, 0.000826754128449799)

# b_M_int_old = -0.0990439890654785
b_M_X = 1.76353928991247
b_M_C1 = 1
b_M_C2 = -1
b_M_int = -sum(b_M_X, b_M_C1, b_M_C2) / 2       # ~ -0.89
b_M = c(b_M_int, b_M_X, b_M_C1, b_M_C2) * scale_factor
# b_M = c(-0.0990439890654785, 1.76353928991247, 0.0128566136999183, 0.00711746366915989)




# Choose theta_Y and theta_M based on the values of b_Y and b_M
theta_Y = c(scale_factor*sqrt(0.5), 0.3, 0.4, scale_factor, 0.5, scale_factor*sqrt(0.8)) / 3
# theta_Y = c(sqrt(0.5), 0.5, 1)
# theta_M = c(1, 0.5, 2)
theta_M = c(scale_factor*sqrt(0.5), -0.5, scale_factor) / 3


all_reg_pars = c(b_Y, theta_Y, b_M, theta_M)




p_Y = length(b_Y)
p_M = length(b_M)
p = p_Y + p_M


folder_suffix = paste0("K=", K, ", N=", N)
dir.create(paste0("R/Paper MC Study/Data - ", folder_suffix), showWarnings = F)
dir.create(paste0("R/Paper MC Study/Results - ", folder_suffix), showWarnings = F)




# Setup cluster
# cl = makeCluster(detectCores() - 2)
# cl = makeCluster(15)
cl = makeCluster(10)
# clusterExport(cl, c("N", "b_Y", "theta_Y", "b_M", "theta_M", "which_REs"))
clusterExport(cl, c("w", "B", "scale", "which_REs", "N", "n", "K", "b_Y", "theta_Y", "b_M", "theta_M", "folder_suffix"))
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
    source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
    source("R/Exact_Asymptotics/Imai Method.r")
    devtools::load_all()
})
clusterSetRNGStream(cl = cl, 123)
# clusterSetRNGStream(cl = cl, 11111111)


num_datasets = 200



# -------------------------- Generate and save data -------------------------- #


set.seed(1)

# First, delete any datasets currently in the target directory
unlink(paste0("R/Paper MC Study/Data - ", folder_suffix, "/*"))

# Generate and save datasets
save_data = pbsapply(1:num_datasets, function(i) {
    data = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
    save(data, file = paste0("R/Paper MC Study/Data - ", folder_suffix, "/", i, ".RData"))
})



# ------------------------ Fit models and save results ----------------------- #


# total_runtime_delta = 0
# total_runtime_MC_delta = 0




# First, delete any results currently in the target directory
unlink(paste0("R/Paper MC Study/Results - ", folder_suffix, "/*"))

# Fit models, extract MEs, estimate covariance matrices and save results
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

        this_time = toc()
        this_timings$get_pars = this_time$toc - this_time$tic
        # data_est = data.frame(hat = Theta_hat, SE = sqrt(diag(cov_hat)))
        # rownames(data_est) = names(Theta_hat)



        # Compute mediation effects
        tic()

        MEs = all_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)
        cov_MEs_delta = all_covs_MEs_pars(scale, w, cov_hat, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

        this_time = toc()
        this_timings$get_MEs = this_time$toc - this_time$tic

        # ------------------------------ MC Delta Method ----------------------------- #
        tic()
        # some_Theta_tildes = sim_Theta_tildes(B, Theta_hat, cov_hat)
        some_Theta_tildes = sim_TMB_Theta_tildes(B, fit_Y, fit_M)
        some_ME_tildes = Theta_tildes_2_MEs(scale, w, some_Theta_tildes, which_REs)
        cov_MEs_MC_delta = cov(some_ME_tildes)

        this_time = toc()
        this_timings$MC_delta = this_time$toc - this_time$tic


        # ------------------------ Compile and return results ------------------------ #
        output = list(this_MEs = MEs, cov_MEs_delta = cov_MEs_delta, cov_MEs_MC_delta = cov_MEs_MC_delta, this_timings = this_timings)

        save(output, file = paste0("R/Paper MC Study/Results - ", folder_suffix, "/", i, ".RData"))
        return(output)
    }, error = function(e){
        output = NULL

        save(output, file = paste0("R/Paper MC Study/Results - ", folder_suffix, "/", i, ".RData"))
        return(output)
    })

    # stop("Error: Loop should never reach this point.")

}, cl = cl)
# })


stopCluster(cl)





#* Build list of all output
output_names = list.files(paste0("R/Paper MC Study/Results - ", folder_suffix, "/"))
MC_results_delta_MC_delta = pblapply(seq_along(output_names), function(x) {
    load(paste0("R/Paper MC Study/Results - ", folder_suffix, "/", x, ".RData"))
    return(output)
})

## Remove NULL entries
MC_results_delta_MC_delta = MC_results_delta_MC_delta[!sapply(MC_results_delta_MC_delta, is.null)]

#* Extract results into separate lists
all_ME_hats = t(sapply(MC_results_delta_MC_delta, function(x) x$this_MEs))
all_cov_hats_delta = lapply(MC_results_delta_MC_delta, function(x) x$cov_MEs_delta)
all_cov_hats_MC_delta = lapply(MC_results_delta_MC_delta, function(x) x$cov_MEs_MC_delta)
all_timings = t(sapply(MC_results_delta_MC_delta, function(x) unlist(x$this_timings)))


#* Compute total time spent on both methods
mean_times = colMeans(all_timings)



#! Get coverage rates using functions from `Paper - GLMM Par Bias.r`
#! Note: Intervals based on the empirical covariance all use the same matrix. Those based on fitted covariances use different matrices for each estimate/dataset

true_MEs = all_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M, which_REs = which_REs)

emp_cov = cov(all_ME_hats)

cover_rate_emp = get_coverage_rates(all_ME_hats, emp_cov, true_MEs)
cover_rate_delta = get_coverage_rates_many_cov_mats(all_ME_hats, all_cov_hats_delta, true_MEs)
cover_rate_MC_delta = get_coverage_rates_many_cov_mats(all_ME_hats, all_cov_hats_MC_delta, true_MEs)

data_cover = data.frame(emp = cover_rate_emp, delta = cover_rate_delta, MC_delta = cover_rate_MC_delta)
rownames(data_cover) = names(true_MEs)
data_cover







# ---------------------------------------------------------------------------- #
#                                      Old                                     #
# ---------------------------------------------------------------------------- #

# # Construct Wald-type CIs
# get_CIs = function(ME_hats, SEs){

#     lower = ME_hats - 1.96 * SEs
#     upper = ME_hats + 1.96 * SEs
#     return(list(lcl = lower, ucl = upper))
# }

# # Empirical covariance matrix of parameter estimates

# load("R/Paper MC Study/true_MEs (new).RData", verbose = TRUE)


# ## Convert all_ME_hats from a list to a data frame
# all_ME_hats_data = do.call(rbind, all_ME_hats)

# ## Compute covariance
# emp_cov = cov(all_ME_hats_data)
# diag(emp_cov)

# ## Wald intervals from empirical sampling vars
# emp_SE_CIs = purrr::map(all_ME_hats, ~ get_CIs(.x, sqrt(diag(emp_cov))))
# emp_SE_CI_checks = lapply(emp_SE_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
# emp_SE_CI_rates = Reduce("+", emp_SE_CI_checks) / length(emp_SE_CI_checks)

# emp_SE_CI_widths = lapply(emp_SE_CIs, function(x) x$ucl - x$lcl)
# emp_SE_mean_widths = Reduce("+", emp_SE_CI_widths) / length(emp_SE_CI_widths)



# # Summaries of estimated SEs (mean of whole matrix and SEs of diagonal entries)

# ## Elementwise means
# mean_delta_SEs = Reduce("+", all_cov_hats_delta) / length(all_cov_hats_delta)
# mean_MC_delta_SEs = Reduce("+", all_cov_hats_MC_delta) / length(all_cov_hats_MC_delta)

# diag(mean_delta_SEs)

# ## SEs of diagonal entries
# ### Compute the SE of the diagonal entries' means from a list of matrices
# diag_mean_SEs = function(some_matrices){
#     some_matrices %>% lapply(diag) %>% do.call(rbind, .) %>% apply(2, sd) / length(all_cov_hats_delta)
# }

# delta_diag_SEs = diag_mean_SEs(all_cov_hats_delta)
# MC_delta_diag_SEs = diag_mean_SEs(all_cov_hats_MC_delta)



# #* MSEs for delta method
# all_delta_errors = sapply(all_cov_hats_delta, function(x) diag(x) - diag(emp_cov))
# all_sq_delta_errors = all_errors^2
# delta_MSEs = apply(all_sq_delta_errors, 1, mean)


# #* Coverage rates for delta method

# load("R/Paper MC Study/true_MEs.RData", verbose = TRUE)

# delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_delta, ~ get_CIs(.x, sqrt(diag(.y))))
# delta_CI_checks = lapply(delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
# delta_CI_rates = Reduce("+", delta_CI_checks) / length(delta_CI_checks)

# delta_CI_widths = lapply(delta_CIs, function(x) x$ucl - x$lcl)
# delta_mean_widths = Reduce("+", delta_CI_widths) / length(delta_CI_widths)



# ## Mean endpoints for CIs
# delta_CI_lcl = Reduce("+", lapply(delta_CIs, function(x) x$lcl)) / length(delta_CIs)
# delta_CI_ucl = Reduce("+", lapply(delta_CIs, function(x) x$ucl)) / length(delta_CIs)
# data_mean_CIs = data.frame(lcl = delta_CI_lcl, truth = true_MEs, ucl = delta_CI_ucl)




# #* Repeat for MC delta method

# ##* MSEs for MC delta method
# all_MC_delta_errors = sapply(all_cov_hats_MC_delta, function(x) diag(x) - diag(emp_cov))
# all_sq_MC_delta_errors = all_MC_delta_errors^2
# MC_delta_MSEs = apply(all_sq_MC_delta_errors, 1, mean)

# ##* Coverage rates for MC delta method
# MC_delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_MC_delta, ~ get_CIs(.x, sqrt(diag(.y))))
# MC_delta_CI_checks = lapply(MC_delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
# MC_delta_CI_rates = Reduce("+", MC_delta_CI_checks) / length(MC_delta_CI_checks)

# MC_delta_CI_widths = lapply(MC_delta_CIs, function(x) x$ucl - x$lcl)
# MC_delta_mean_widths = Reduce("+", MC_delta_CI_widths) / length(MC_delta_CI_widths)


# cover_data = data.frame(empirical = emp_SE_CI_rates, delta = delta_CI_rates, MC_delta = MC_delta_CI_rates)
# width_data = data.frame(empirical = emp_SE_mean_widths, delta = delta_mean_widths, MC_delta = MC_delta_mean_widths)

# CI_data = cbind(cover_data, width_data)
# apply(CI_data, 2, mean)






# #* Check coverage rates as a function of number of datasets
# data_delta_CI_checks = Reduce("rbind", delta_CI_checks)
# cumul_delta_CI_rates = apply(data_delta_CI_checks, 2, cumsum) / (1:nrow(data_delta_CI_checks))

# for(i in seq_len(ncol(cumul_delta_CI_rates))){
#     plot(cumul_delta_CI_rates[, i], type = "l", ylim = c(0.75, 1), main = paste0("Coverage Rates for i=", i, " (", colnames(cumul_delta_CI_rates)[i], ")"))
#     abline(h = cumul_delta_CI_rates[nrow(data_delta_CI_checks), i], col="red")
# }




#  # ---------------------------------------------------------------------------- #
#  #                          Repeat analysis with half K                         #
#  # ---------------------------------------------------------------------------- #
# data_names = list.files("R/Paper MC Study/Results (half K) - Delta, MC Delta/")
# MC_results_delta_MC_delta = pblapply(seq_along(data_names), function(x) {
#     load(paste0("R/Paper MC Study/Results - Delta, MC Delta/", x, ".RData"))
#     return(output)
# })

# ## Remove NULL entries
# MC_results_delta_MC_delta = MC_results_delta_MC_delta[!sapply(MC_results_delta_MC_delta, function(x) is.null(x$this_MEs))]

# #* Extract results into separate lists
# all_ME_hats = lapply(MC_results_delta_MC_delta, function(x) x$this_MEs)
# all_cov_hats_delta = lapply(MC_results_delta_MC_delta, function(x) x$this_delta_cov)
# all_cov_hats_MC_delta = lapply(MC_results_delta_MC_delta, function(x) x$this_MC_delta_cov)
# runtime_total_delta = sapply(MC_results_delta_MC_delta, function(x) x$this_delta_runtime) %>% sum()
# runtime_total_MC_delta = sapply(MC_results_delta_MC_delta, function(x) x$this_MC_delta_runtime) %>% sum()


# # save(all_ME_hats, all_cov_hats_delta, all_cov_hats_MC_delta, runtime_total_delta, runtime_total_MC_delta, file = "R/Paper MC Study/Paper - Delta and MC Delta Results (half K).RData")
# load("R/Paper MC Study/Paper - Delta and MC Delta Results (half K).RData", verbose = TRUE)


# # Estimated mediation effects
# all_ME_hats_data = do.call(rbind, all_ME_hats)

# delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_delta, ~ get_CIs(.x, sqrt(diag(.y))))
# delta_CI_checks = lapply(delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
# delta_CI_rates = Reduce("+", delta_CI_checks) / length(delta_CI_checks)

# MC_delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_MC_delta, ~ get_CIs(.x, sqrt(diag(.y))))
# MC_delta_CI_checks = lapply(MC_delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
# MC_delta_CI_rates = Reduce("+", MC_delta_CI_checks) / length(MC_delta_CI_checks)


# small_K_delta_rates = delta_CI_rates
# small_K_MC_delta_rates = MC_delta_CI_rates


# delta_CI_rates
# data.frame(large_delta = delta_CI_rates, small_delta = small_K_delta_rates, large_MC_delta = MC_delta_CI_rates, small_MC_delta = small_K_MC_delta_rates)
# delta_rate_rats = delta_CI_rates / small_K_delta_rates
# MC_delta_rate_rats = MC_delta_CI_rates / small_K_MC_delta_rates







# # ---------------------------------------------------------------------------- #
# #                          Repeat analysis with half n                         #
# # ---------------------------------------------------------------------------- #
# data_names = list.files("R/Paper MC Study/Results (half n) - Delta, MC Delta/")
# MC_results_delta_MC_delta = pblapply(seq_along(data_names), function(x) {
#     load(paste0("R/Paper MC Study/Results - Delta, MC Delta/", x, ".RData"))
#     return(output)
# })

# ## Remove NULL entries
# MC_results_delta_MC_delta = MC_results_delta_MC_delta[!sapply(MC_results_delta_MC_delta, function(x) is.null(x$this_MEs))]

# #* Extract results into separate lists
# all_ME_hats = lapply(MC_results_delta_MC_delta, function(x) x$this_MEs)
# all_cov_hats_delta = lapply(MC_results_delta_MC_delta, function(x) x$this_delta_cov)
# all_cov_hats_MC_delta = lapply(MC_results_delta_MC_delta, function(x) x$this_MC_delta_cov)
# runtime_total_delta = sapply(MC_results_delta_MC_delta, function(x) x$this_delta_runtime) %>% sum()
# runtime_total_MC_delta = sapply(MC_results_delta_MC_delta, function(x) x$this_MC_delta_runtime) %>% sum()


# # save(all_ME_hats, all_cov_hats_delta, all_cov_hats_MC_delta, runtime_total_delta, runtime_total_MC_delta, file = "R/Paper MC Study/Paper - Delta and MC Delta Results (half n).RData")
# load("R/Paper MC Study/Paper - Delta and MC Delta Results (half n).RData", verbose = TRUE)


# # Estimated mediation effects
# all_ME_hats_data = do.call(rbind, all_ME_hats)

# delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_delta, ~ get_CIs(.x, sqrt(diag(.y))))
# delta_CI_checks = lapply(delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
# delta_CI_rates = Reduce("+", delta_CI_checks) / length(delta_CI_checks)

# MC_delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_MC_delta, ~ get_CIs(.x, sqrt(diag(.y))))
# MC_delta_CI_checks = lapply(MC_delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
# MC_delta_CI_rates = Reduce("+", MC_delta_CI_checks) / length(MC_delta_CI_checks)


# small_K_delta_rates = delta_CI_rates
# small_K_MC_delta_rates = MC_delta_CI_rates


# delta_CI_rates
# data.frame(large_delta = delta_CI_rates, small_delta = small_K_delta_rates, large_MC_delta = MC_delta_CI_rates, small_MC_delta = small_K_MC_delta_rates)
# delta_rate_rats = delta_CI_rates / small_K_delta_rates
# MC_delta_rate_rats = MC_delta_CI_rates / small_K_MC_delta_rates
