


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

# N = 500
N=1000
n = N


# k = 100
K=10



# Vector of covariates
## Note that mediation effects currently require us to specify values for the confounders. I think the next move is to average these conditional effects over the observed distribution of confounders.
w = c(2,3)



## Non-trivial values for the b's and theta's. Former based on output from another MC study. Latter chosen arbitrarily.
## Crucially, no parameters are equal to zero.
##? We choose the intercepts to that the mean of the linear predictor is zero. Doing this for M makes it easier to do so for Y.

# Scale factor for coefficients and SDs
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



# Covariance parameters for random effects. See above note for details
theta_Y = c(scale_factor*sqrt(0.5), 0.3, 0.4, scale_factor, 0.5, scale_factor*sqrt(0.8)) / 3
theta_M = c(scale_factor*sqrt(0.5), -0.5, scale_factor) / 3


# This is also called Theta
all_reg_pars = c(b_Y, theta_Y, b_M, theta_M)




p_Y = length(b_Y)
p_M = length(b_M)
p = p_Y + p_M






# ----------------- Tell script where to find data and output ---------------- #
folder_suffix = paste0("K=", K, ", N=", N)
dir.create(paste0("R/For_Bouchra/Data - ", folder_suffix), showWarnings = F)
dir.create(paste0("R/For_Bouchra/Results - ", folder_suffix), showWarnings = F)







# ------------------------------- Setup cluster ------------------------------ #

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
    source("R/For_Bouchra/Exact_Asymptotics_Helpers.r")
    source("R/For_Bouchra/Imai Method.r")
    devtools::load_all()
})
clusterSetRNGStream(cl = cl, 123)
# clusterSetRNGStream(cl = cl, 11111111)


num_datasets = 200




# ---------------------------------------------------------------------------- #
#                            Generate and save data                            #
# ---------------------------------------------------------------------------- #


set.seed(1)

# First, delete any datasets currently in the target directory
unlink(paste0("R/For_Bouchra/Data - ", folder_suffix, "/*"))

# Generate and save datasets
save_data = pbsapply(1:num_datasets, function(i) {
    data = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
    save(data, file = paste0("R/For_Bouchra/Data - ", folder_suffix, "/", i, ".RData"))
})





# ---------------------------------------------------------------------------- #
#                          Fit models and save results                         #
# ---------------------------------------------------------------------------- #


# total_runtime_delta = 0
# total_runtime_MC_delta = 0




# First, delete any results currently in the target directory
unlink(paste0("R/For_Bouchra/Results - ", folder_suffix, "/*"))

# Fit models, extract MEs, estimate covariance matrices and save results
MC_results_delta_MC_delta = pblapply(1:num_datasets, function(i) {
# MC_results_delta_MC_delta = pblapply(1:3, function(i) {
    load(paste0("R/For_Bouchra/Data - ", folder_suffix, "/", i, ".RData"))


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



        # Compute mediation effects
        tic()

        MEs = all_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)
        cov_MEs_delta = all_covs_MEs_pars(scale, w, cov_hat, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

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

        save(output, file = paste0("R/For_Bouchra/Results - ", folder_suffix, "/", i, ".RData"))
        return(output)
    }, error = function(e){
        #? Sometimes, the optimization fails. I'm not too worried about this, because in practice you would just try a different optimization algorithm. I can't really do that for a Monte Carlo Study, however. Furthermore, this only really comes up for small values of k.
        output = NULL

        save(output, file = paste0("R/For_Bouchra/Results - ", folder_suffix, "/", i, ".RData"))
        return(output)
    })


}, cl = cl)
# })


stopCluster(cl)





# ---------------------------------------------------------------------------- #
#                                Clean-Up Output                               #
# ---------------------------------------------------------------------------- #


#* Build list of all output
output_names = list.files(paste0("R/For_Bouchra/Results - ", folder_suffix, "/"))
MC_results_delta_MC_delta = pblapply(seq_along(output_names), function(x) {
    load(paste0("R/For_Bouchra/Results - ", folder_suffix, "/", x, ".RData"))
    return(output)
})

## Remove NULL entries (see above note in the `error` part of the tryCatch function)
MC_results_delta_MC_delta = MC_results_delta_MC_delta[!sapply(MC_results_delta_MC_delta, is.null)]


#* Extract results into separate lists
all_ME_hats = t(sapply(MC_results_delta_MC_delta, function(x) x$this_MEs))
all_cov_hats_delta = lapply(MC_results_delta_MC_delta, function(x) x$cov_MEs_delta)
all_cov_hats_MC_delta = lapply(MC_results_delta_MC_delta, function(x) x$cov_MEs_MC_delta)
all_timings = t(sapply(MC_results_delta_MC_delta, function(x) unlist(x$this_timings)))




# ---------------------------------------------------------------------------- #
#                                Analyze Output                                #
# ---------------------------------------------------------------------------- #


#* Compute total time spent on both methods
mean_times = colMeans(all_timings)



#* Get coverage rates
# Note: Intervals based on the empirical covariance all use the same matrix. Those based on fitted covariances use different matrices for each estimate/dataset

true_MEs = all_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M, which_REs = which_REs)

emp_cov = cov(all_ME_hats)

cover_rate_emp = get_coverage_rates(all_ME_hats, emp_cov, true_MEs)
cover_rate_delta = get_coverage_rates_many_cov_mats(all_ME_hats, all_cov_hats_delta, true_MEs)
cover_rate_MC_delta = get_coverage_rates_many_cov_mats(all_ME_hats, all_cov_hats_MC_delta, true_MEs)

data_cover = data.frame(emp = cover_rate_emp, delta = cover_rate_delta, MC_delta = cover_rate_MC_delta)
rownames(data_cover) = names(true_MEs)
data_cover


