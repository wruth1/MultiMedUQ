



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
library(mediation)
library(outliers)       # Grubbs' test for outliers
source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
source("R/Exact_Asymptotics/Imai Method.r")
source("R/Paper MC Study/Outlier Functions.r")

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

num_bin_confounders = 1
num_cont_confounders = 1
total_num_confounders = num_bin_confounders + num_cont_confounders

B = 500     # Number of samples to generate for MC delta
scale = c("diff", "rat", "OR")      # What scales should we compute mediation effects on?
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")        # Which variables have random effects? Eventually, I will need a better way to specify this

# Number of groups
# K = 10
# K = 50
K = 100
# K = 200
# K = 500

# Observations per group
# N = 100
N = 500
# N=1000
# N = 10000
n = N


# # Values of K and n for testing
# K = 74
# N = 432
# n = N







# Vector of confounders. Used as point at which to evaluate conditional effects
w = rep(0.5, times = total_num_confounders)



## Non-trivial values for the b's and theta's. Former based on output from another MC study. Latter chosen arbitrarily.
## Crucially, no parameters are equal to zero.
##? We choose the intercepts to that the mean of the linear predictor is zero. Doing this for M makes it easier to do so for Y.
#! Scale factor for coefficients and SDs
scale_factor = 1

set.seed(123)

b_Y_X = 0.966486302988689
b_Y_M = 1.99644760563721
b_Y_Cs = sample(c(-1, 1), total_num_confounders, replace = T)
b_Y_int = - sum(b_Y_X, b_Y_M, b_Y_Cs) / 2       # ~ -1.48
b_Y = c(b_Y_int, b_Y_X, b_Y_M, b_Y_Cs) * scale_factor


set.seed(321)

b_M_X = 1.76353928991247
b_M_Cs = sample(c(-1, 1), total_num_confounders, replace = T)
b_M_int = -sum(b_M_X, b_M_Cs) / 2       # ~ -0.89
b_M = c(b_M_int, b_M_X, b_M_Cs) * scale_factor




# Choose theta_Y and theta_M based on the values of b_Y and b_M
theta_Y = c(scale_factor*sqrt(0.5), 0.3, 0.4, scale_factor, 0.5, scale_factor*sqrt(0.8)) / 3
theta_M = c(scale_factor*sqrt(0.5), -0.5, scale_factor) / 3


all_reg_pars = c(b_Y, theta_Y, b_M, theta_M)




p_Y = length(b_Y)
p_M = length(b_M)
p = p_Y + p_M


folder_suffix = paste0("K=", K, ", N=", N, ", conf=", num_bin_confounders, ",", num_cont_confounders)
dir.create(paste0("R/Paper MC Study/Data/Data - ", folder_suffix), showWarnings = F)
dir.create(paste0("R/Paper MC Study/Results/Conditional_MEs/Results - ", folder_suffix), showWarnings = F)
# dir.create(paste0("R/Paper MC Study/Results/Marginal_MEs_only_delta/Results - ", folder_suffix), showWarnings = F)



outcome_name = "Y"
exposure_name = "X"
mediator_name = "M"
group_name = "group"






# cl = makeCluster(detectCores() - 2)
# cl = makeCluster(15)
cl = makeCluster(10)
# clusterExport(cl, c("N", "b_Y", "theta_Y", "b_M", "theta_M", "which_REs"))
clusterExport(cl, c("w", "B", "scale", "which_REs", "N", "n", "K", "b_Y", "theta_Y", "b_M", "theta_M", "folder_suffix", "get_confounder_freq", "outcome_name", "exposure_name", "mediator_name", "group_name", "num_confounders"))
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
    library(mediation)
    source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
    source("R/Exact_Asymptotics/Imai Method.r")
    devtools::load_all()
})
clusterSetRNGStream(cl = cl, 123)
# clusterSetRNGStream(cl = cl, 11111111)





# -------------------------- Generate and save data -------------------------- #

num_datasets = 200

# First, delete any datasets currently in the target directory
unlink(paste0("R/Paper MC Study/Data/Data - ", folder_suffix, "/*"))

set.seed(1)

# Generate and save datasets
save_data = pbsapply(1:num_datasets, function(i) {
    data = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, num_bin_confounders = 1, num_cont_confounders = 1, output_list = F, which_REs = which_REs)
    save(data, file = paste0("R/Paper MC Study/Data/Data - ", folder_suffix, "/", i, ".RData"))
})#, cl=cl)





# ---------------------------------------------------------------------------- #
#                                     ENCs                                     #
# ---------------------------------------------------------------------------- #

clusterSetRNGStream(cl = cl, 123)


num_datasets = length(list.files(paste0("R/Paper MC Study/Data/Data - ", folder_suffix)))
# num_datasets = 10

unlink(paste0("R/Paper MC Study/Results/Conditional_MEs/Results - ", folder_suffix, "/*"))
# unlink(paste0("R/Paper MC Study/Results/Marginal_MEs_only_delta/Results - ", folder_suffix, "/*"))

# Fit models, extract ENCs, estimate covariance matrices and save results
MC_results_ME = pblapply(1:num_datasets, function(i) {
# MC_results_ME = pblapply(1:10, function(i) {
    load(paste0("R/Paper MC Study/Data/Data - ", folder_suffix, "/", i, ".RData"), verbose = T)



    tryCatch({

        this_timings = list()


        # ---------------------------- Delta Method (ours) --------------------------- #

        # Fit models
        tic()

        ## Build formulas for model fitting
        ### Y ~ X + M + C1 + ... + Cr + (X + M | group)
        pred_vars_Y <- setdiff(names(data), c("Y", "group"))
        formula_Y <- as.formula(paste("Y ~", paste(pred_vars_Y, collapse = " + "), " + (X + M | group)"))

        ### M ~ X + C1 + ... + Cr + (X | group)
        pred_vars_M <- setdiff(names(data), c("Y", "M", "group"))
        formula_M <- as.formula(paste("M ~", paste(pred_vars_M, collapse = " + "), " + (X | group)"))

        Y_fit_warning = NULL
        M_fit_warning = NULL

        #? Some model fits produce warnings. I want to capture these for later investigation
        fit_Y = withCallingHandlers(glmmTMB(formula_Y, data = data, family = binomial), 
                                     warning = function(w) {
                                         Y_fit_warning <<- w
                                     })
        fit_M = withCallingHandlers(glmmTMB(formula_M, data = data, family = binomial), 
                                     warning = function(w) {
                                         M_fit_warning <<- w
                                     }) 

        model_fit_warnings = list(Y = Y_fit_warning, M = M_fit_warning)

        # fit_Y_BFGS = glmmTMB(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = glmmTMBControl(optimizer = "optim", optArgs = list(method = "BFGS", eval.max = 1e10)))
        # fit_M_BFGS = glmmTMB(M ~ X + C1 + C2 + (X | group), data = data, family = binomial)#, control = glmmTMBControl(optimizer = "optim", optArgs = list(method = "BFGS", eval.max = 1e10)))

        this_time = toc()
        this_timings$fit_models = this_time$toc - this_time$tic


        # Extract fitted parameters
        tic()

        theta_hat_Y = get_model_pars_TMB(fit_Y)
        theta_hat_M = get_model_pars_TMB(fit_M)
        Theta_hat = c(unlist(theta_hat_Y), unlist(theta_hat_M))
        cov_hat = all_pars_cov_mat_TMB(fit_Y, fit_M)

        # eigen(cov_hat)$values

        # cbind(Theta_hat, (diag(cov_hat)))

        b_Y_hat = theta_hat_Y[["b"]]
        theta_Y_hat = theta_hat_Y[["theta"]]
        b_M_hat = theta_hat_M[["b"]]
        theta_M_hat = theta_hat_M[["theta"]]
        len_par_vecs = sapply(list(b_Y_hat, theta_Y_hat, b_M_hat, theta_M_hat), length)

        this_time = toc()
        this_timings$get_pars = this_time$toc - this_time$tic




        # Compute mediation effects
        tic()

        MEs = all_MEs_pars(scale, w, b_Y_hat, theta_Y_hat, b_M_hat, theta_M_hat, which_REs =  which_REs)
        cov_MEs_delta = all_covs_MEs_pars(scale, w, cov_hat, b_Y_hat, theta_Y_hat, b_M_hat, theta_M_hat, which_REs =  which_REs)

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

        save(output, file = paste0("R/Paper MC Study/Results/Conditional_MEs/Results - ", folder_suffix, "/", i, ".RData"))
        return(output)
    }, error = function(e){
        output = NULL

        save(output, file = paste0("R/Paper MC Study/Results/Conditional_MEs/Results - ", folder_suffix, "/", i, ".RData"))
        return(output)
    })


}, cl = cl)
# })

print(folder_suffix)

stopCluster(cl)







# ---------------------------------------------------------------------------- #
#                        Analysis of Simulation Results                        #
# ---------------------------------------------------------------------------- #



K = 74
N = 432
folder_suffix = paste0("K=", K, ", N=", N, ", conf=", num_bin_confounders, ",", num_cont_confounders)

# ------------------------------ Read-in Results ----------------------------- #
output_names = list.files(paste0("R/Paper MC Study/Results/Conditional_MEs/Results - ", folder_suffix, "/"))
MC_results_ME = pblapply(output_names, function(x) {
    load(paste0("R/Paper MC Study/Results/Conditional_MEs/Results - ", folder_suffix, "/", x))
    return(output)
})



# ----------------------------- Clean-up Results ----------------------------- #

# check_lengths = sapply(MC_results_ME, length)
# which(check_lengths < 4)
# MC_results_ME = remove_indices(MC_results_ME, which(check_lengths < 4))
# length(MC_results_ME)

#* Backup raw output from simulation
MC_results_ME_raw = MC_results_ME
# MC_results_ME = MC_results_ME_raw

#* Check which entries are null
check_null = sapply(MC_results_ME_raw, is.null)
inds_null = which(check_null)

#* Remove null entries
MC_results_ME = MC_results_ME %>% Filter(Negate(is.null), .)
length(MC_results_ME)

#* Remove iterations with negative variances
all_covs_raw = lapply(MC_results_ME, function(x) x$cov_MEs_delta)
inds_remove = sapply(all_covs_raw, function(x) any(diag(x) < 0))
any(inds_remove)
MC_results_ME = MC_results_ME[!inds_remove]
length(MC_results_ME)

cov_hat_ew_mins = sapply(all_covs_raw, function(x) min(Re(eigen(x)$values)))
# hist(cov_hat_ew_mins)
# hist(cov_hat_ew_mins[inds_remove])
min(cov_hat_ew_mins)




# ------------------------------ Extract results ----------------------------- #

all_warnings = lapply(MC_results_ME, function(x) x$model_fit_warnings)
warnings_check = sapply(all_warnings, function(x) !all(is.null(unlist(x))))
inds_warnings = which(warnings_check)
inds_warnings
length(inds_warnings)

MC_results_ME = MC_results_ME[!warnings_check]


all_MEs_raw = sapply(MC_results_ME, function(x) x$this_MEs) %>% t()
all_covs_raw = lapply(MC_results_ME, function(x) x$cov_MEs_delta)
all_cov_tildes_raw = lapply(MC_results_ME, function(x) x$cov_MEs_MC_delta)



# Timings
all_timings = t(sapply(MC_results_ME, function(x) unlist(x$this_timings)))
mean_times = colMeans(all_timings)
mean_times
sd_times = apply(all_timings, 2, sd)
sd_times

# ## Group timings by method
# both_methods_timings = all_timings[,1] + all_timings[,2] + all_timings[,3]
# delta_only_timings = all_timings[,4] + all_timings[,5]
# MC_delta_only_timings = all_timings[,6]

# mean_both_methods_timings = mean(both_methods_timings)
# mean_delta_only_timings = mean(delta_only_timings)
# mean_MC_delta_only_timings = mean(MC_delta_only_timings)
# mean_timings = c(mean_both_methods_timings, mean_delta_only_timings, mean_MC_delta_only_timings)
# names(mean_timings) = c("Both Methods", "Delta Only", "MC-Delta Only")
# mean_timings




emp_cov = cov(all_MEs_raw)

mean_cov_hat_raw = Reduce("+", all_covs_raw) / length(all_covs_raw)
mean_cov_hat_raw = (mean_cov_hat_raw + t(mean_cov_hat_raw)) / 2

mean_cov_tilde_raw = Reduce("+", all_cov_tildes_raw) / length(all_cov_tildes_raw)
mean_cov_tilde_raw = (mean_cov_tilde_raw + t(mean_cov_tilde_raw)) / 2




# --------- Remove runs with pathologically large covariance matrices -------- #

#? Significance level threshold for declaring outliers under Grubbs' Test
sig_level = 1e-10

#* Get all deviations from empirical covariance for both delta and MC-delta
all_cov_hat_errs = sapply(all_covs_raw, function(x) mat_rel_err(emp_cov, x))
all_cov_tilde_errs = sapply(all_cov_tildes_raw, function(x) mat_rel_err(emp_cov, x))



#* Identify runs corresponding to error outliers
bad_runs_delta = find_all_outliers(all_cov_hat_errs, sig_level = sig_level)
bad_runs_MC_delta = find_all_outliers(all_cov_tilde_errs, sig_level = sig_level)
all_bad_runs = unique(c(bad_runs_delta, bad_runs_MC_delta))
# all_bad_runs = bad_runs_delta

length(all_bad_runs)
length(all_bad_runs) / length(all_covs_raw)


#* Remove runs corresponding to error outliers
all_MEs = remove_rows(all_MEs_raw, all_bad_runs)
all_covs = remove_indices(all_covs_raw, all_bad_runs)
all_cov_tildes = remove_indices(all_cov_tildes_raw, all_bad_runs)


#* Diagnostic plots for covariance matrix outliers
hist(all_cov_hat_errs)
# hist(all_cov_hat_errs[-which.max(all_cov_hat_errs)])
# which(all_cov_hat_errs > 100)
# all_covs[[58]]
# all_cov_hat_errs %>% Filter(function(x) x < 15, .) %>% hist(breaks = 20)


hist(all_cov_tilde_errs)
# hist(all_cov_tilde_errs[-which.max(all_cov_tilde_errs)])
# which(all_cov_tilde_errs > 3)
# all_cov_tildes[[65]]
# all_cov_tilde_errs %>% Filter(function(x) x < 3, .) %>% hist(breaks = 20)




# ---------------------------- Get Coverage Rates ---------------------------- #

true_MEs = all_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M, which_REs = which_REs)



#* Coverage rates before removing outliers
cover_rate_emp_raw = get_coverage_rates(all_MEs_raw, emp_cov, true_MEs)
cover_rate_cov_hat_raw = get_coverage_rates_many_cov_mats(all_MEs_raw, all_covs_raw, true_MEs)
cover_rate_cov_tilde_raw = get_coverage_rates_many_cov_mats(all_MEs_raw, all_cov_tildes_raw, true_MEs)

data_cover_rate_raw = data.frame(
    emp = cover_rate_emp_raw,
    cov_hat = cover_rate_cov_hat_raw,
    cov_tilde = cover_rate_cov_tilde_raw
)
rownames(data_cover_rate_raw) = colnames(all_MEs_raw)
data_cover_rate_raw

colMeans(data_cover_rate_raw)

# Relative coverage rates compared to emp
data_cover_rate_rel_raw = data_cover_rate_raw / cover_rate_emp_raw
data_cover_rate_rel_raw

delta_vs_MC_delta_cover_rate_raw = cover_rate_cov_hat_raw / cover_rate_cov_tilde_raw
delta_vs_MC_delta_cover_rate_raw



#* Coverage rates after removing outliers
cover_rate_emp = get_coverage_rates(all_MEs, emp_cov, true_MEs)
cover_rate_cov_hat = get_coverage_rates_many_cov_mats(all_MEs, all_covs, true_MEs)
cover_rate_cov_tilde = get_coverage_rates_many_cov_mats(all_MEs, all_cov_tildes, true_MEs)

data_cover_rate = data.frame(
    emp = cover_rate_emp,
    cov_hat = cover_rate_cov_hat,
    cov_tilde = cover_rate_cov_tilde
)
rownames(data_cover_rate) = colnames(all_MEs)
data_cover_rate


# Relative coverage rates compared to emp
data_cover_rate_rel = data_cover_rate / cover_rate_emp
data_cover_rate_rel

delta_vs_MC_delta_cover_rate = cover_rate_cov_hat / cover_rate_cov_tilde
delta_vs_MC_delta_cover_rate



# ------------------------------- Get CI Widths ------------------------------ #

#* Before removing outliers
mean_CI_width_emp_raw = mean_widths_one_cov_mat(all_MEs_raw, emp_cov)
mean_CI_width_cov_hat_raw = mean_widths_many_cov_mats(all_MEs_raw, all_covs_raw)
mean_CI_width_cov_tilde_raw = mean_widths_many_cov_mats(all_MEs_raw, all_cov_tildes_raw)

data_CI_width_raw = data.frame(
    emp = mean_CI_width_emp_raw,
    cov_hat = mean_CI_width_cov_hat_raw,
    cov_tilde = mean_CI_width_cov_tilde_raw
)
rownames(data_CI_width_raw) = colnames(all_MEs_raw)
data_CI_width_raw

# Average widths relative to empirical
data_CI_rel_widths_raw = data_CI_width_raw / mean_CI_width_emp_raw
data_CI_rel_widths_raw
apply(data_CI_rel_widths_raw, 2, function(x) max(abs(x - 1)))


#* After removing outliers
mean_CI_width_emp = mean_widths_one_cov_mat(all_MEs, emp_cov)
mean_CI_width_cov_hat = mean_widths_many_cov_mats(all_MEs, all_covs)
mean_CI_width_cov_tilde = mean_widths_many_cov_mats(all_MEs, all_cov_tildes)

data_CI_width = data.frame(
    emp = mean_CI_width_emp,
    cov_hat = mean_CI_width_cov_hat,
    cov_tilde = mean_CI_width_cov_tilde
)
rownames(data_CI_width) = colnames(all_MEs)
data_CI_width

#* Average widths relative to empirical
data_CI_rel_widths = data_CI_width / mean_CI_width_emp
data_CI_rel_widths
apply(data_CI_rel_widths, 2, function(x) max(abs(x - 1)))

