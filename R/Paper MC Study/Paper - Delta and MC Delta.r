


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

# w = c(2,3)
w = c(1,1)
B = 500
scale = c("diff", "rat", "OR")
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")



# Setup cluster
# cl = makeCluster(detectCores() - 2)
# cl = makeCluster(15)
cl = makeCluster(10)
# clusterExport(cl, c("N", "b_Y", "theta_Y", "b_M", "theta_M", "which_REs"))
clusterExport(cl, c("w", "B", "scale", "which_REs", "n", "K"))
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
    source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
    source("R/Exact_Asymptotics/Imai Method.r")
    devtools::load_all()
})
clusterSetRNGStream(cl = cl, 123)



# all_ME_hats = list()
# all_cov_hats_delta = list()
# all_cov_hats_MC_delta = list()

# total_runtime_delta = 0
# total_runtime_MC_delta = 0

# #? Note: Extracting the SE matrix for the reg pars is slow (using merDeriv::vcov.glmerMod()). I don't want to duplicate this step (or fitting the models), so I separated out model fitting/SE extraction from my method and the MC delta.


tic()
# MC_results_delta_MC_delta = pblapply(1:num_reps, function(i) {
MC_results_delta_MC_delta = pblapply(1:500, function(i) {
    load(paste0("R/Paper MC Study/Datasets/", i, ".RData"))

    # #! Remove half the groups
    # groups_keep = paste0("G", 1:(K / 2))
    # data %<>% dplyr::filter(group %in% groups_keep)

    # #! Remove half the samples in each group
    # data %<>% dplyr::group_by(group) %>%
    #     dplyr::slice_sample(n = n / 2) %>% dplyr::ungroup()


    tryCatch({

        ## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
        (fit_Y = suppressMessages(lme4::glmer(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))
        (fit_M = suppressMessages(lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))


        # Compute parameter estimates
        this_MEs = all_MEs_models(w = w, fit_Y = fit_Y, fit_M = fit_M)

        ## Extract model parameter estimates
        info_Y = get_model_pars(fit_Y)
        info_M = get_model_pars(fit_M)
        this_Theta_hat = c(unlist(info_Y), unlist(info_M))
        this_cov_hat = all_pars_cov_mat(fit_Y, fit_M)

        # Delta Method covariance
        time = Sys.time()

        this_delta_cov = all_covs_MEs_models(w = w, Sigma = this_cov_hat, fit_Y = fit_Y, fit_M = fit_M)

        this_delta_runtime = Sys.time() - time
        # total_runtime_delta = total_runtime_delta + this_delta_runtime



        # MC delta method covariance
        time = Sys.time()

        some_Theta_tildes = sim_Theta_tildes(B, this_Theta_hat, this_cov_hat)
        some_ME_tildes = Theta_tildes_2_MEs(scale, w, some_Theta_tildes, which_REs)
        this_MC_delta_cov = cov(some_ME_tildes)

        this_MC_delta_runtime = Sys.time() - time
        # total_runtime_MC_delta = total_runtime_MC_delta + this_MC_delta_runtime

        }, error = function(e){
            this_MEs = NULL
            this_delta_cov = NULL
            this_MC_delta_cov = NULL
            this_delta_runtime = NULL
            this_MC_delta_runtime = NULL
    })

    tryCatch({
    output = list(this_MEs = this_MEs, this_delta_cov = this_delta_cov, this_MC_delta_cov = this_MC_delta_cov, this_delta_runtime = this_delta_runtime, this_MC_delta_runtime = this_MC_delta_runtime)

    save(output, file = paste0("R/Paper MC Study/Results (new) - Delta, MC Delta/", i, ".RData"))
    return(output)
    }, error = function(e){
      output = list(this_MEs = NULL, this_delta_cov = NULL, this_MC_delta_cov = NULL, this_delta_runtime = NULL, this_MC_delta_runtime = NULL)

      save(output, file = paste0("R/Paper MC Study/Results (new) - Delta, MC Delta/", i, ".RData"))
      return(output)
    })

    # stop("Error: Loop should never reach this point.")

}, cl = cl)
time_info = toc()
runtime_total = time_info$toc - time_info$tic

stopCluster(cl)





#* Build list of all output
data_names = list.files("R/Paper MC Study/Results (new) - Delta, MC Delta/")
MC_results_delta_MC_delta = pblapply(seq_along(data_names), function(x) {
    load(paste0("R/Paper MC Study/Results (new) - Delta, MC Delta/", x, ".RData"))
    return(output)
})

## Remove NULL entries
MC_results_delta_MC_delta = MC_results_delta_MC_delta[!sapply(MC_results_delta_MC_delta, function(x) is.null(x$this_MEs))]

#* Extract results into separate lists
all_ME_hats = lapply(MC_results_delta_MC_delta, function(x) x$this_MEs)
all_cov_hats_delta = lapply(MC_results_delta_MC_delta, function(x) x$this_delta_cov)
all_cov_hats_MC_delta = lapply(MC_results_delta_MC_delta, function(x) x$this_MC_delta_cov)
runtime_total_delta = sapply(MC_results_delta_MC_delta, function(x) x$this_delta_runtime) %>% sum()
runtime_total_MC_delta = sapply(MC_results_delta_MC_delta, function(x) x$this_MC_delta_runtime) %>% sum()


save(all_ME_hats, all_cov_hats_delta, all_cov_hats_MC_delta, runtime_total_delta, runtime_total_MC_delta, file = "R/Paper MC Study/Paper - Delta and MC Delta Results (new).RData")
# load("R/Paper MC Study/Paper - Delta and MC Delta Results.RData", verbose = TRUE)


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

load("R/Paper MC Study/true_MEs (new).RData", verbose = TRUE)


## Convert all_ME_hats from a list to a data frame
all_ME_hats_data = do.call(rbind, all_ME_hats)

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

load("R/Paper MC Study/true_MEs.RData", verbose = TRUE)

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


cover_data = data.frame(empirical = emp_SE_CI_rates, delta = delta_CI_rates, MC_delta = MC_delta_CI_rates)
width_data = data.frame(empirical = emp_SE_mean_widths, delta = delta_mean_widths, MC_delta = MC_delta_mean_widths)

CI_data = cbind(cover_data, width_data)
apply(CI_data, 2, mean)






#* Check coverage rates as a function of number of datasets
data_delta_CI_checks = Reduce("rbind", delta_CI_checks)
cumul_delta_CI_rates = apply(data_delta_CI_checks, 2, cumsum) / (1:nrow(data_delta_CI_checks))

for(i in seq_len(ncol(cumul_delta_CI_rates))){
    plot(cumul_delta_CI_rates[, i], type = "l", ylim = c(0.75, 1), main = paste0("Coverage Rates for i=", i, " (", colnames(cumul_delta_CI_rates)[i], ")"))
    abline(h = cumul_delta_CI_rates[nrow(data_delta_CI_checks), i], col="red")
}




 # ---------------------------------------------------------------------------- #
 #                          Repeat analysis with half K                         #
 # ---------------------------------------------------------------------------- #
data_names = list.files("R/Paper MC Study/Results (half K) - Delta, MC Delta/")
MC_results_delta_MC_delta = pblapply(seq_along(data_names), function(x) {
    load(paste0("R/Paper MC Study/Results - Delta, MC Delta/", x, ".RData"))
    return(output)
})

## Remove NULL entries
MC_results_delta_MC_delta = MC_results_delta_MC_delta[!sapply(MC_results_delta_MC_delta, function(x) is.null(x$this_MEs))]

#* Extract results into separate lists
all_ME_hats = lapply(MC_results_delta_MC_delta, function(x) x$this_MEs)
all_cov_hats_delta = lapply(MC_results_delta_MC_delta, function(x) x$this_delta_cov)
all_cov_hats_MC_delta = lapply(MC_results_delta_MC_delta, function(x) x$this_MC_delta_cov)
runtime_total_delta = sapply(MC_results_delta_MC_delta, function(x) x$this_delta_runtime) %>% sum()
runtime_total_MC_delta = sapply(MC_results_delta_MC_delta, function(x) x$this_MC_delta_runtime) %>% sum()


# save(all_ME_hats, all_cov_hats_delta, all_cov_hats_MC_delta, runtime_total_delta, runtime_total_MC_delta, file = "R/Paper MC Study/Paper - Delta and MC Delta Results (half K).RData")
load("R/Paper MC Study/Paper - Delta and MC Delta Results (half K).RData", verbose = TRUE)


# Estimated mediation effects
all_ME_hats_data = do.call(rbind, all_ME_hats)

delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_delta, ~ get_CIs(.x, sqrt(diag(.y))))
delta_CI_checks = lapply(delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
delta_CI_rates = Reduce("+", delta_CI_checks) / length(delta_CI_checks)

MC_delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_MC_delta, ~ get_CIs(.x, sqrt(diag(.y))))
MC_delta_CI_checks = lapply(MC_delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
MC_delta_CI_rates = Reduce("+", MC_delta_CI_checks) / length(MC_delta_CI_checks)


small_K_delta_rates = delta_CI_rates
small_K_MC_delta_rates = MC_delta_CI_rates


delta_CI_rates
data.frame(large_delta = delta_CI_rates, small_delta = small_K_delta_rates, large_MC_delta = MC_delta_CI_rates, small_MC_delta = small_K_MC_delta_rates)
delta_rate_rats = delta_CI_rates / small_K_delta_rates
MC_delta_rate_rats = MC_delta_CI_rates / small_K_MC_delta_rates







# ---------------------------------------------------------------------------- #
#                          Repeat analysis with half n                         #
# ---------------------------------------------------------------------------- #
data_names = list.files("R/Paper MC Study/Results (half n) - Delta, MC Delta/")
MC_results_delta_MC_delta = pblapply(seq_along(data_names), function(x) {
    load(paste0("R/Paper MC Study/Results - Delta, MC Delta/", x, ".RData"))
    return(output)
})

## Remove NULL entries
MC_results_delta_MC_delta = MC_results_delta_MC_delta[!sapply(MC_results_delta_MC_delta, function(x) is.null(x$this_MEs))]

#* Extract results into separate lists
all_ME_hats = lapply(MC_results_delta_MC_delta, function(x) x$this_MEs)
all_cov_hats_delta = lapply(MC_results_delta_MC_delta, function(x) x$this_delta_cov)
all_cov_hats_MC_delta = lapply(MC_results_delta_MC_delta, function(x) x$this_MC_delta_cov)
runtime_total_delta = sapply(MC_results_delta_MC_delta, function(x) x$this_delta_runtime) %>% sum()
runtime_total_MC_delta = sapply(MC_results_delta_MC_delta, function(x) x$this_MC_delta_runtime) %>% sum()


# save(all_ME_hats, all_cov_hats_delta, all_cov_hats_MC_delta, runtime_total_delta, runtime_total_MC_delta, file = "R/Paper MC Study/Paper - Delta and MC Delta Results (half n).RData")
load("R/Paper MC Study/Paper - Delta and MC Delta Results (half n).RData", verbose = TRUE)


# Estimated mediation effects
all_ME_hats_data = do.call(rbind, all_ME_hats)

delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_delta, ~ get_CIs(.x, sqrt(diag(.y))))
delta_CI_checks = lapply(delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
delta_CI_rates = Reduce("+", delta_CI_checks) / length(delta_CI_checks)

MC_delta_CIs = purrr::map2(all_ME_hats, all_cov_hats_MC_delta, ~ get_CIs(.x, sqrt(diag(.y))))
MC_delta_CI_checks = lapply(MC_delta_CIs, function(x) x$lcl < true_MEs & true_MEs < x$ucl)
MC_delta_CI_rates = Reduce("+", MC_delta_CI_checks) / length(MC_delta_CI_checks)


small_K_delta_rates = delta_CI_rates
small_K_MC_delta_rates = MC_delta_CI_rates


delta_CI_rates
data.frame(large_delta = delta_CI_rates, small_delta = small_K_delta_rates, large_MC_delta = MC_delta_CI_rates, small_MC_delta = small_K_MC_delta_rates)
delta_rate_rats = delta_CI_rates / small_K_delta_rates
MC_delta_rate_rats = MC_delta_CI_rates / small_K_MC_delta_rates
