

#? In this script, we adjust the estimated MEs to be unbiased, based on the mean estimates over our 1000 datasets.
#? We then re-compute the CIs and check coverage.
#? Note: This bias correction doesn't trivialize the problem, since individual datasets' estimates will still deviate from the mean.


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
n = all_datasets[[1]] %>% filter(group == "G1") %>% nrow()
K = all_datasets[[1]] %>% select(group) %>% unique() %>% nrow()

w = c(2,3)
B = 500
scale = c("diff", "rat", "OR")
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")





# save(all_ME_hats, all_cov_hats_delta, all_cov_hats_MC_delta, runtime_total_delta, runtime_total_MC_delta, file = "R/Paper MC Study/Paper - Delta and MC Delta Results.RData")
load("R/Paper MC Study/Paper - Delta and MC Delta Results.RData", verbose = TRUE)


# #! Delta method has one pathological entry. Let's remove it and proceed
# delta_traces = sapply(all_cov_hats_delta, function(x) sum(diag(x)))
# which_err = which(delta_traces > 10000)
# all_ME_hats = all_ME_hats[-which_err]
# all_cov_hats_delta = all_cov_hats_delta[-which_err]
# all_cov_hats_MC_delta = all_cov_hats_MC_delta[-which_err]

# Construct Wald-type CIs
get_CIs = function(ME_hats, SEs){

    lower = ME_hats - 1.96 * SEs
    upper = ME_hats + 1.96 * SEs
    return(list(lcl = lower, ucl = upper))
}

# Empirical covariance matrix of parameter estimates

load("R/Paper MC Study/true_MEs.RData", verbose = TRUE)


## Convert all_ME_hats from a list to a data frame
all_ME_hats_data = do.call(rbind, all_ME_hats)
ME_hat_means = colMeans(all_ME_hats_data)
ME_hat_biases = ME_hat_means - true_MEs

all_ME_hats = lapply(all_ME_hats, function(x) x - ME_hat_biases)

# all_ME_hats_data =t( t(all_ME_hats_data) -  ME_hat_biases)

## Compute covariance
emp_cov = cov(all_ME_hats_data)
diag(emp_cov)

## Wald intervals from empirical sampling vars
emp_SE_CIs = purrr::map(all_ME_hats, ~ get_CIs(.x, sqrt(diag(emp_cov))))
emp_SE_CI_checks = lapply(emp_SE_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
emp_SE_CI_rates = Reduce("+", emp_SE_CI_checks) / length(emp_SE_CI_checks)

emp_SE_CI_widths = lapply(emp_SE_CIs, function(x) x$ucl - x$lcl)
emp_SE_mean_widths = Reduce("+", emp_SE_CI_widths) / length(emp_SE_CI_widths)



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


delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_delta, ~ get_CIs(.x, sqrt(diag(.y))))
delta_CI_checks = lapply(delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
delta_CI_rates = Reduce("+", delta_CI_checks) / length(delta_CI_checks)

delta_CI_widths = lapply(delta_CIs, function(x) x$ucl - x$lcl)
delta_mean_widths = Reduce("+", delta_CI_widths) / length(delta_CI_widths)



## Mean endpoints for CIs
delta_CI_lcl = Reduce("+", lapply(delta_CIs, function(x) x$lcl)) / length(delta_CIs)
delta_CI_ucl = Reduce("+", lapply(delta_CIs, function(x) x$ucl)) / length(delta_CIs)
data_mean_CIs = data.frame(lcl = delta_CI_lcl, truth = true_MEs, ucl = delta_CI_ucl)




#* Repeat for MC delta method

##* MSEs for MC delta method
all_MC_delta_errors = sapply(all_cov_hats_MC_delta, function(x) diag(x) - diag(emp_cov))
all_sq_MC_delta_errors = all_MC_delta_errors^2
MC_delta_MSEs = apply(all_sq_MC_delta_errors, 1, mean)

##* Coverage rates for MC delta method
MC_delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_MC_delta, ~ get_CIs(.x, sqrt(diag(.y))))
MC_delta_CI_checks = lapply(MC_delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
MC_delta_CI_rates = Reduce("+", MC_delta_CI_checks) / length(MC_delta_CI_checks)

MC_delta_CI_widths = lapply(MC_delta_CIs, function(x) x$ucl - x$lcl)
MC_delta_mean_widths = Reduce("+", MC_delta_CI_widths) / length(MC_delta_CI_widths)


cover_data = data.frame(empirical = emp_SE_CI_rates, delta = delta_CI_rates, MC_delta = MC_delta_CI_rates, delta_over_MC_delta = delta_CI_rates / MC_delta_CI_rates)
width_data = data.frame(empirical = emp_SE_mean_widths, delta = delta_mean_widths, MC_delta = MC_delta_mean_widths)

CI_data = cbind(cover_data, width_data)
apply(CI_data, 2, mean)






#* Check coverage rates as a function of number of datasets
data_delta_CI_checks = Reduce("rbind", delta_CI_checks)
cumul_delta_CI_rates = apply(data_delta_CI_checks, 2, cumsum) / (1:nrow(data_delta_CI_checks))

for(i in seq_len(ncol(cumul_delta_CI_rates))){
    plot(cumul_delta_CI_rates[, i], type = "l", ylim = c(0.75, 1), main = paste0("Coverage Rates for i=", i, " (", colnames(cumul_delta_CI_rates)[i], ")"))
    abline(h = cumul_delta_CI_rates[nrow(data_delta_CI_checks), i], col="red")
    abline(h = 0.95, col = "red", lty = 2)

    # Add confidence bands at +/- 2 sqrt(p ( 1 - p) / n)
    this_UCLs = cumul_delta_CI_rates[, i] + 1.96 * sqrt(cumul_delta_CI_rates[, i] * (1 - cumul_delta_CI_rates[, i]) / (1:nrow(data_delta_CI_checks)))
    this_LCLs = cumul_delta_CI_rates[, i] - 1.96 * sqrt(cumul_delta_CI_rates[, i] * (1 - cumul_delta_CI_rates[, i]) / (1:nrow(data_delta_CI_checks)))

    lines(this_LCLs, col = "blue")
    lines(this_UCLs, col = "blue")
}








# ---------------------------------------------------------------------------- #
#                          Repeat analysis with half K                         #
# ---------------------------------------------------------------------------- #
load("R/Paper MC Study/Paper - Delta and MC Delta Results (half K).RData", verbose = TRUE)


# Estimated mediation effects
all_ME_hats_data = do.call(rbind, all_ME_hats)
ME_hat_means = colMeans(all_ME_hats_data)
ME_hat_biases = ME_hat_means - true_MEs
all_ME_hats = lapply(all_ME_hats, function(x) x - ME_hat_biases)

small_K_delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_delta, ~ get_CIs(.x, sqrt(diag(.y))))
small_K_delta_CI_checks = lapply(small_K_delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
small_K_delta_rates = Reduce("+", small_K_delta_CI_checks) / length(small_K_delta_CI_checks)

small_K_MC_delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_MC_delta, ~ get_CIs(.x, sqrt(diag(.y))))
small_K_MC_delta_CI_checks = lapply(small_K_MC_delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
small_K_MC_delta_rates = Reduce("+", small_K_MC_delta_CI_checks) / length(small_K_MC_delta_CI_checks)




delta_CI_rates
data.frame(large_delta = delta_CI_rates, small_delta = small_K_delta_rates, large_MC_delta = MC_delta_CI_rates, small_MC_delta = small_K_MC_delta_rates)
(delta_rate_rats = delta_CI_rates / small_K_delta_rates)
(MC_delta_rate_rats = MC_delta_CI_rates / small_K_MC_delta_rates)




# ---------------------------------------------------------------------------- #
#                          Repeat analysis with half n                         #
# ---------------------------------------------------------------------------- #
# save(all_ME_hats, all_cov_hats_delta, all_cov_hats_MC_delta, runtime_total_delta, runtime_total_MC_delta, file = "R/Paper MC Study/Paper - Delta and MC Delta Results (half n).RData")
load("R/Paper MC Study/Paper - Delta and MC Delta Results (half n).RData", verbose = TRUE)


# Estimated mediation effects
all_ME_hats_data = do.call(rbind, all_ME_hats)
ME_hat_means = colMeans(all_ME_hats_data)
ME_hat_biases = ME_hat_means - true_MEs
all_ME_hats = lapply(all_ME_hats, function(x) x - ME_hat_biases)


small_n_delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_delta, ~ get_CIs(.x, sqrt(diag(.y))))
small_n_delta_CI_checks = lapply(small_n_delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
small_n_delta_rates = Reduce("+", small_n_delta_CI_checks) / length(small_n_delta_CI_checks)

small_n_MC_delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_MC_delta, ~ get_CIs(.x, sqrt(diag(.y))))
small_n_MC_delta_CI_checks = lapply(small_n_MC_delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
small_n_MC_delta_rates = Reduce("+", small_n_MC_delta_CI_checks) / length(small_n_MC_delta_CI_checks)



delta_CI_rates
data.frame(large_delta = delta_CI_rates, small_K_delta = small_K_delta_rates, small_n_delta = small_n_delta_rates, large_MC_delta = MC_delta_CI_rates, small_K_MC_delta = small_K_MC_delta_rates, small_n_MC_delta = small_n_MC_delta_rates)
delta_rate_rats = delta_CI_rates / small_K_delta_rates
MC_delta_rate_rats = MC_delta_CI_rates / small_K_MC_delta_rates
