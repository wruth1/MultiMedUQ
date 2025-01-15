


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
# source("R/Paper MC Study/glmmTMB Helpers.r")
devtools::load_all("D:/William/Research/MultiMedUQ")



load("R/Paper MC Study/all_datasets.RData", verbose = TRUE)

num_reps = length(all_datasets)
n = all_datasets[[1]] %>% filter(group == "G1") %>% nrow()
K = all_datasets[[1]] %>% select(group) %>% unique() %>% nrow()

# w = c(2,3)
w = c(1,1)
B = 500
scale = c("diff", "rat", "OR")
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")





# Set parameters and fit models

# N = 100
# N = 20
# N = 40
# N = 60
# N=1000

#* Main value
N = 100
n = N

# all_Ks = c(50, 100, 200, 400, 800)
# all_Ks = c(50, 100, 200)
# all_Ks = 50 * (2:6)
# K=1000

#* Main value
K = 200

# num_reps = 30
# num_reps = 500
num_reps = 1000


# which_REs = c("Y.Int", "Y.X", "M.All")
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")



x = 0
x_m = 1



w = c(2,3)


# ## Non-trivial values for the b's and theta's. Former based on output from another MC study. Latter chosen arbitrarily.
# ## Crucially, no parameters are equal to zero.
# ##? We choose the intercepts to that the mean of the linear predictor is zero. Doing this for M makes it easier to do so for Y.
# # b_Y_int_old = 0.0376828219852018
# b_Y_X = 0.966486302988689
# b_Y_M = 1.99644760563721
# b_Y_C1 = -0.00556557712859059
# b_Y_C2 = 0.000826754128449799
# b_Y_int = - sum(b_Y_X, b_Y_M, b_Y_C1, b_Y_C2) / 2       # ~ -1.48
# b_Y = c(b_Y_int, b_Y_X, b_Y_M, b_Y_C1, b_Y_C2)
# # b_Y = c(0.0376828219852018, 0.966486302988689, 1.99644760563721, -0.00556557712859059, 0.000826754128449799)

# # b_M_int_old = -0.0990439890654785
# b_M_X = 1.76353928991247
# b_M_C1 = 0.0128566136999183
# b_M_C2 = 0.00711746366915989
# b_M_int = -sum(b_M_X, b_M_C1, b_M_C2) / 2       # ~ -0.89
# b_M = c(b_M_int, b_M_X, b_M_C1, b_M_C2)
# # b_M = c(-0.0990439890654785, 1.76353928991247, 0.0128566136999183, 0.00711746366915989)



# #! Replace extremely small coefficients on C1 and C2
# ## Non-trivial values for the b's and theta's. Former based on output from another MC study. Latter chosen arbitrarily.
# ## Crucially, no parameters are equal to zero.
# ##? We choose the intercepts to that the mean of the linear predictor is zero. Doing this for M makes it easier to do so for Y.
# # b_Y_int_old = 0.0376828219852018
# b_Y_X = 0.966486302988689
# b_Y_M = 1.99644760563721
# b_Y_C1 = -1
# b_Y_C2 = 1
# b_Y_int = - sum(b_Y_X, b_Y_M, b_Y_C1, b_Y_C2) / 2       # ~ -1.48
# b_Y = c(b_Y_int, b_Y_X, b_Y_M, b_Y_C1, b_Y_C2)
# # b_Y = c(0.0376828219852018, 0.966486302988689, 1.99644760563721, -0.00556557712859059, 0.000826754128449799)

# # b_M_int_old = -0.0990439890654785
# b_M_X = 1.76353928991247
# b_M_C1 = 1
# b_M_C2 = -1
# b_M_int = -sum(b_M_X, b_M_C1, b_M_C2) / 2       # ~ -0.89
# b_M = c(b_M_int, b_M_X, b_M_C1, b_M_C2)
# # b_M = c(-0.0990439890654785, 1.76353928991247, 0.0128566136999183, 0.00711746366915989)




# # Choose theta_Y and theta_M based on the values of b_Y and b_M
# theta_Y = c(sqrt(0.5), 0.3, 0.4, 1, 0.5, sqrt(0.8)) / 3
# # theta_Y = c(sqrt(0.5), 0.5, 1)
# # theta_M = c(1, 0.5, 2)
# theta_M = c(sqrt(0.5), -0.5, 1) / 3


# all_reg_pars = c(b_Y, theta_Y, b_M, theta_M)




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



# Setup cluster
# cl = makeCluster(detectCores() - 2)
# cl = makeCluster(15)
cl = makeCluster(10)
# clusterExport(cl, c("N", "b_Y", "theta_Y", "b_M", "theta_M", "which_REs"))
clusterExport(cl, c("w", "B", "scale", "which_REs", "N", "n", "K", "b_Y", "theta_Y", "b_M", "theta_M"))
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



# all_ME_hats = list()
# all_cov_hats_delta = list()
# all_cov_hats_MC_delta = list()

# total_runtime_delta = 0
# total_runtime_MC_delta = 0

# #? Note: Extracting the SE matrix for the reg pars is slow (using merDeriv::vcov.glmerMod()). I don't want to duplicate this step (or fitting the models), so I separated out model fitting/SE extraction from my method and the MC delta.



# MC_results_delta_MC_delta = pblapply(1:num_reps, function(i) {
all_sim_results = pblapply(1:50, function(i) {
    data = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
    # load(paste0("R/Paper MC Study/Datasets/", i, ".RData"))

    tryCatch({

        # ----------------------------------- lme4 ----------------------------------- #

        ## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
        # (fit_Y = suppressMessages(lme4::glmer(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))
        # (fit_M = suppressMessages(lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))

        # info_Y = get_model_pars(fit_Y)
        # info_M = get_model_pars(fit_M)
        # this_Theta_hat = c(unlist(info_Y), unlist(info_M))



        # ---------------------------------- glmmTMB --------------------------------- #
        fit_Y = glmmTMB(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial)
        fit_M = glmmTMB(M ~ X + C1 + C2 + (X | group), data = data, family = binomial)


        ## Extract model parameter estimates
        info_Y = get_model_pars_TMB(fit_Y)
        info_M = get_model_pars_TMB(fit_M)
        this_Theta_hat = c(unlist(info_Y), unlist(info_M))

        ## Estimate standard errors
        cov_hat_Y = TMB_2_GLMM_SE(fit_Y)
        cov_hat_M = TMB_2_GLMM_SE(fit_M)


        

        # ----------------------------- Mediation Effects ---------------------------- #
        ## Compute MEs
        b_Y_hat = info_Y$b
        theta_Y_hat = info_Y$theta
        b_M_hat = info_M$b
        theta_M_hat = info_M$theta

        this_MEs = all_MEs_pars(scale, w, b_Y_hat, theta_Y_hat, b_M_hat, theta_M_hat, which_REs = which_REs)


        }, error = function(e){
            this_Theta_hat = NULL
            this_MEs = NULL
    })

    tryCatch({
    output = list(this_Theta_hat = this_Theta_hat, this_MEs = this_MEs)

    # save(output, file = paste0("R/Paper MC Study/Results (new) - Delta, MC Delta/", i, ".RData"))
    return(output)
    }, error = function(e){
      output = list(this_Theta_hat = NULL, this_MEs = NULL)

    #   save(output, file = paste0("R/Paper MC Study/Results (new) - Delta, MC Delta/", i, ".RData"))
      return(output)
    })

    # stop("Error: Loop should never reach this point.")

# })
}, cl = cl)

stopCluster(cl)

all_theta_hats = lapply(all_sim_results, function(x) x[["this_Theta_hat"]])
all_MEs = lapply(all_sim_results, function(x) x[["this_MEs"]])
all_theta_cov_hats = lapply(all_sim_results, function(x) x[["this_Theta_cov"]])

# MC_results_delta_MC_delta = MC_results_delta_MC_delta_old

mean_reg_pars = MC_results_delta_MC_delta %>%
  lapply(function(x) x[["this_Theta_hat"]]) %>%
  Reduce(f = "rbind") %>%
  colMeans()

# load("R/Paper MC Study/True GLMM Pars.RData", verbose = TRUE)
bias_reg_pars = mean_reg_pars - all_reg_pars
rel_bias_reg_pars = abs(bias_reg_pars / all_reg_pars)


data_bias_reg_pars = data.frame(par = names(mean_reg_pars), hat = mean_reg_pars, true = all_reg_pars, rel_bias = rel_bias_reg_pars)
# new_output = data_bias_reg_pars
# new_data = MC_results_delta_MC_delta





# Bias in MEs
mean_MEs = MC_results_delta_MC_delta %>%
  lapply(function(x) x[["this_MEs"]]) %>%
  Reduce(f = "rbind") %>%
  colMeans()

true_MEs = all_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M, which_REs = which_REs)
bias_MEs = mean_MEs - true_MEs
rel_bias_MEs = abs(bias_MEs / true_MEs)

data_bias_MEs = data.frame(par = names(mean_MEs), hat = mean_MEs, true = true_MEs, rel_bias = rel_bias_MEs)




MC_results_lme4 = MC_results_delta_MC_delta




# I did double-check that this is the correct point at which to evaluate the Jacobian of the gradient
## Specifically, I evaluated fit$obj$fn at the point estimates. This matched the final objective function value reported in fit$fit.
estimates = fit$fit$par
grad_fun = fit$obj$gr
SE_mat = solve(numDeriv::jacobian(grad_fun, estimates))






TMB_pars = fit$fit$par
TMB_FE_pars = TMB_pars[names(TMB_pars)=="beta"]
TMB_RE_pars = TMB_pars[names(TMB_pars)=="theta"]

# Map from TMB internal parameterization to my parameterization of RE parameters
## Theirs is log-SDs, followed by scaled Cholesky factors of the correlation matrix. See https://github.com/glmmTMB/glmmTMB/blob/master/misc/glmmTMB_corcalcs.ipynb for details on how their correlation factors are defined and how to do the mapping.

num_RE_pars = length(TMB_RE_pars)
num_vars = (sqrt(1 + 8*num_RE_pars) - 1 ) / 2

TMB_SD_pars = TMB_RE_pars[1:num_vars]
TMB_corr_pars = TMB_RE_pars[(num_vars+1):length(TMB_RE_pars)]

TMB_SDs = exp(TMB_SD_pars)

## Compute correlation matrix
TMB_corrs = glmmTMB::get_cor(TMB_corr_pars, return_val = "vec")





# Gradient of map from TMB internal parameterization to my parameterization of RE parameters
TMB_pars = fit$fit$par
TMB_FE_pars = TMB_pars[names(TMB_pars)=="beta"]
TMB_RE_pars = TMB_pars[names(TMB_pars)=="theta"]

num_pars = length(TMB_pars)
num_FE_pars = length(TMB_FE_pars)

num_RE_pars = length(TMB_RE_pars)
num_vars = (sqrt(1 + 8*num_RE_pars) - 1 ) / 2
num_SD_pars = num_vars
num_corr_pars = num_RE_pars - num_vars

TMB_SD_pars = TMB_RE_pars[1:num_vars]
TMB_corr_pars = TMB_RE_pars[(num_vars+1):length(TMB_RE_pars)]



## Jacobian of fixed effects (FE) parameters 
d_FE_FE = diag(num_FE_pars)
d_FE_RE = matrix(0, nrow = num_FE_pars, ncol = num_RE_pars)
d_FE = cbind(d_FE_FE, d_FE_RE)

## Jacobian of re-parameterization of SD parameters
d_SD_FE = matrix(0, nrow = num_SD_pars, ncol = num_FE_pars)
d_SD_SD = diag(exp(TMB_SD_pars))
d_SD_corr = matrix(0, nrow = num_SD_pars, ncol = num_corr_pars)
d_SD = cbind(d_SD_FE, d_SD_SD, d_SD_corr)

## Jacobian of re-parameterization of correlation parameters
d_corr_FE = matrix(0, nrow = num_corr_pars, ncol = num_FE_pars)
d_corr_SD = matrix(0, nrow = num_corr_pars, ncol = num_SD_pars)
d_corr_corr = t(numDeriv::jacobian(glmmTMB::get_cor, TMB_corr_pars, return_val = "vec")) # jacobian output demension is range-by-domain
d_corr = cbind(d_corr_FE, d_corr_SD, d_corr_corr)

## Combined Jacobian
d_TMB = rbind(d_FE, d_SD, d_corr)



## Re-arrange d_SD and d_corr to get d_theta, the gradient of my parameterization of the RE parameters
### Specifically, my theta is sd_1, corr_12, corr_13,..., sd_2, corr_23, ..., sd_n
d_theta = matrix(0, nrow = num_pars, ncol = num_pars)

## Fill-in gradient for fixed effects
for(i in 1:num_FE_pars){
    d_theta[i,] = d_TMB[i,]
}

# Leave space between SDs for correlations. This quantity was derived by hand.
inds_SD = num_FE_pars + ((1:num_SD_pars) - 1) * (1+num_SD_pars - (1:num_SD_pars)/2) + 1

## Fill-in gradients for SD and corr pars
for(i in seq_len(num_SD_pars)){
    ### SDs
    theta_ind_SD = inds_SD[i]
    d_theta[theta_ind_SD,] = d_SD[i,]

    ### Correlations
    this_num_corrs = num_SD_pars - i    # Number of correlations for this variable
    num_corrs_so_far = (i-1)*(num_SD_pars - i/2)    # Number of correlations already accounted for (this helps us index d_corr)
    for(j in seq_len(this_num_corrs)){
        theta_ind_corr = theta_ind_SD + j   # Index in my theta of current correlation

        this_corr_ind = num_corrs_so_far + j    # Index in d_corr of current correlation

        d_theta[theta_ind_corr,] = d_corr[this_corr_ind,]
    }
}

d_theta - TMB_2_GLMM_grad(fit)


SE_mat_TMB = d_theta %*% SE_mat %*% t(d_theta)
TMB_2_GLMM_SE(fit) - SE_mat_TMB



 SE_mat_lme4 = merDeriv::vcov.glmerMod(fit_Y, full=TRUE, ranpar = "sd")


norm(SE_mat_lme4 - SE_mat_TMB, type = "2") / norm(SE_mat_lme4)

fixef(fit_Y)
fixef(fit)$cond

theta_lme4 = get_model_pars(fit_Y)$theta
theta_TMB = c(TMB_SDs[1],
              TMB_corrs[1:2],
              TMB_SDs[2],
              TMB_corrs[3],
              TMB_SDs[3])

cbind(theta_lme4, theta_TMB)




# ---------------------------------------------------------------------------- #
#                            MC Study of Single GLMM                           #
# ---------------------------------------------------------------------------- #

N = 100
K = 200


# Setup cluster
# cl = makeCluster(detectCores() - 2)
# cl = makeCluster(15)
cl = makeCluster(10)
# clusterExport(cl, c("N", "b_Y", "theta_Y", "b_M", "theta_M", "which_REs"))
clusterExport(cl, c("w", "B", "scale", "which_REs", "N", "n", "K", "b_Y", "theta_Y", "b_M", "theta_M"))
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
    # source("R/Paper MC Study/glmmTMB Helpers.r")
    devtools::load_all()
})
clusterSetRNGStream(cl = cl, 123)
# clusterSetRNGStream(cl = cl, 11111111)



# all_ME_hats = list()
# all_cov_hats_delta = list()
# all_cov_hats_MC_delta = list()

# total_runtime_delta = 0
# total_runtime_MC_delta = 0

# #? Note: Extracting the SE matrix for the reg pars is slow (using merDeriv::vcov.glmerMod()). I don't want to duplicate this step (or fitting the models), so I separated out model fitting/SE extraction from my method and the MC delta.


# library(optimx)
# library(dfoptim)
# Y_info_lme4 = allFit(fit_Y_lme4)
# summary(Y_info_lme4)


# MC_results_delta_MC_delta = pblapply(1:num_reps, function(i) {
all_sim_results = pblapply(1:100, function(i) {
    data = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
    # load(paste0("R/Paper MC Study/Datasets/", i, ".RData"))


    # ------------------------------- TMB Analysis ------------------------------- #
    fit_Y_TMB = glmmTMB(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial)
    fit_M_TMB = glmmTMB(M ~ X + C1 + C2 + (X | group), data = data, family = binomial)

    theta_hat_Y_TMB = get_model_pars_TMB(fit_Y_TMB)
    theta_hat_M_TMB = get_model_pars_TMB(fit_M_TMB)
    Theta_hat_TMB = c(unlist(theta_hat_Y_TMB), unlist(theta_hat_M_TMB))
    cov_hat_TMB = all_pars_cov_mat_TMB(fit_Y_TMB, fit_M_TMB)

    b_Y_TMB = theta_hat_Y_TMB[["b"]]
    theta_Y_TMB = theta_hat_Y_TMB[["theta"]]
    b_M_TMB = theta_hat_M_TMB[["b"]]
    theta_M_TMB = theta_hat_M_TMB[["theta"]]
    MEs_TMB = all_MEs_pars(scale, w, b_Y_TMB, theta_Y_TMB, b_M_TMB, theta_M_TMB, which_REs =  which_REs)
    cov_MEs_TMB = all_covs_MEs_pars(scale, w, cov_hat_TMB, b_Y_TMB, theta_Y_TMB, b_M_TMB, theta_M_TMB, which_REs =  which_REs)


    # ------------------------------- lme4 Analysis ------------------------------ #
    fit_Y_lme4 = suppressMessages(lme4::glmer(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "nlminbwrap", optCtrl = list(maxfun = 1e5))))
    fit_M_lme4 = suppressMessages(lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "nloptwrap", optCtrl = list(maxfun = 1e5, algorithm = "NLOPT_LN_NELDERMEAD"))))

    theta_hat_Y_lme4 = get_model_pars(fit_Y_lme4)
    theta_hat_M_lme4 = get_model_pars(fit_M_lme4)
    Theta_hat_lme4 = c(unlist(theta_hat_Y_lme4), unlist(theta_hat_M_lme4))
    cov_hat_lme4 = all_pars_cov_mat(fit_Y_lme4, fit_M_lme4)

    MEs_lme4 = all_MEs_models(scale, w, fit_Y_lme4, fit_M_lme4, which_REs =  which_REs)
    cov_MEs_lme4 = all_covs_MEs_models(scale, w, cov_hat_lme4, fit_Y_lme4, fit_M_lme4, which_REs =  which_REs)
    

    output = list(Theta_hat_TMB = Theta_hat_TMB, cov_hat_TMB = cov_hat_TMB, MEs_TMB = MEs_TMB, cov_MEs_TMB = cov_MEs_TMB,
                    Theta_hat_lme4 = Theta_hat_lme4, cov_hat_lme4 = cov_hat_lme4, MEs_lme4 = MEs_lme4, cov_MEs_lme4 = cov_MEs_lme4)
    save(output, file = paste0("R/Paper MC Study/Results - GLMM Par Bias/", i, ".RData"))

    return(output)
}, cl = cl)
# })

stopCluster(cl)


#* Load results and store in a single object
all_sim_results = lapply(1:10, function(i) {
    load(paste0("R/Paper MC Study/Results - GLMM Par Bias/", i, ".RData"))
    return(output)
})


#* Extract each type of estimate
## Specifically: Theta hat and ME hat, as well as their SEs, for both TMB and lme4
all_Theta_hats_TMB = t(sapply(all_sim_results, function(x) x$Theta_hat_TMB))
all_Theta_cov_hats_TMB = lapply(all_sim_results, function(x) x$cov_hat_TMB)
all_ME_hats_TMB = t(sapply(all_sim_results, function(x) x$MEs_TMB))
all_ME_cov_hats_TMB = lapply(all_sim_results, function(x) x$cov_MEs_TMB)

all_Theta_hats_lme4 = t(sapply(all_sim_results, function(x) x$Theta_hat_lme4))
all_Theta_cov_hats_lme4 = lapply(all_sim_results, function(x) x$cov_hat_lme4)
all_ME_hats_lme4 = t(sapply(all_sim_results, function(x) x$MEs_lme4))
all_ME_cov_hats_lme4 = lapply(all_sim_results, function(x) x$cov_MEs_lme4)


#* Compare findings from TMB and lme4

## Empirical covariance of Theta hat

emp_cov_Theta_TMB = cov(all_Theta_hats_TMB)
emp_cov_Theta_lme4 = cov(all_Theta_hats_lme4)
norm(emp_cov_Theta_TMB - emp_cov_Theta_lme4) / norm(emp_cov_Theta_TMB)

all_theta_hats_TMB = t(sapply(all_sim_results, function(x) x$theta_hat_Y_TMB))
all_cov_hats_TMB = lapply(all_sim_results, function(x) x$cov_hat_Y_TMB)
all_theta_hats_lme4 = t(sapply(all_sim_results, function(x) x$theta_hat_Y_lme4))

# Compare empirical and mean estimated covariance matrices
emp_cov_TMB = cov(all_theta_hats_TMB)
mean_cov_hat_TMB = Reduce("+", all_cov_hats_TMB) / length(all_cov_hats_TMB)
norm(emp_cov_TMB - mean_cov_hat_TMB) / norm(emp_cov_TMB)

# Compare empirical and mean estimated variances (i.e. diagonal entries of above matrices)
emp_vars = as.matrix(as.numeric(diag(emp_cov)))
mean_vars = as.matrix(diag(mean_cov_hat))
norm(emp_vars - mean_vars, type="2") / norm(emp_vars, type = "2")


emp_cov_lme4 = all_theta_hats_lme4 %>% na.omit() %>%  cov()
norm(emp_cov_lme4 - emp_cov_TMB, type = "2") / norm(emp_cov_lme4, type = "2")
