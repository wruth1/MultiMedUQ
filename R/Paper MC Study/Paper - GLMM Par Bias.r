


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
N=200 

#* Main value
# N = 100
n = N

# all_Ks = c(50, 100, 200, 400, 800)
# all_Ks = c(50, 100, 200)
# all_Ks = 50 * (2:6)
K=1000

#* Main value
# K = 200

# num_reps = 30
# num_reps = 500
num_reps = 1000


# which_REs = c("Y.Int", "Y.X", "M.All")
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")


x = 0
x_m = 1



w = c(2,3)


## Non-trivial values for the b's and theta's. Former based on output from another MC study. Latter chosen arbitrarily.
## Crucially, no parameters are equal to zero.
##? We choose the intercepts to that the mean of the linear predictor is zero. Doing this for M makes it easier to do so for Y.
# b_Y_int_old = 0.0376828219852018
b_Y_X = 0.966486302988689
b_Y_M = 1.99644760563721
b_Y_C1 = -0.00556557712859059
b_Y_C2 = 0.000826754128449799
b_Y_int = - sum(b_Y_X, b_Y_M, b_Y_C1, b_Y_C2) / 2       # ~ -1.48
b_Y = c(b_Y_int, b_Y_X, b_Y_M, b_Y_C1, b_Y_C2)
# b_Y = c(0.0376828219852018, 0.966486302988689, 1.99644760563721, -0.00556557712859059, 0.000826754128449799)

# b_M_int_old = -0.0990439890654785
b_M_X = 1.76353928991247
b_M_C1 = 0.0128566136999183
b_M_C2 = 0.00711746366915989
b_M_int = -sum(b_M_X, b_M_C1, b_M_C2) / 2       # ~ -0.89
b_M = c(b_M_int, b_M_X, b_M_C1, b_M_C2)
# b_M = c(-0.0990439890654785, 1.76353928991247, 0.0128566136999183, 0.00711746366915989)



## Non-trivial values for the b's and theta's. Former based on output from another MC study. Latter chosen arbitrarily.
## Crucially, no parameters are equal to zero.
##? We choose the intercepts to that the mean of the linear predictor is zero. Doing this for M makes it easier to do so for Y.
# b_Y_int_old = 0.0376828219852018
b_Y_X = 0.966486302988689
b_Y_M = 1.99644760563721
b_Y_C1 = -1
b_Y_C2 = 1
b_Y_int = - sum(b_Y_X, b_Y_M, b_Y_C1, b_Y_C2) / 2       # ~ -1.48
b_Y = c(b_Y_int, b_Y_X, b_Y_M, b_Y_C1, b_Y_C2)
# b_Y = c(0.0376828219852018, 0.966486302988689, 1.99644760563721, -0.00556557712859059, 0.000826754128449799)

# b_M_int_old = -0.0990439890654785
b_M_X = 1.76353928991247
b_M_C1 = 1
b_M_C2 = -1
b_M_int = -sum(b_M_X, b_M_C1, b_M_C2) / 2       # ~ -0.89
b_M = c(b_M_int, b_M_X, b_M_C1, b_M_C2)
# b_M = c(-0.0990439890654785, 1.76353928991247, 0.0128566136999183, 0.00711746366915989)




# Choose theta_Y and theta_M based on the values of b_Y and b_M
theta_Y = c(sqrt(0.5), 0.5, 0.5, 1, 0.5, sqrt(0.5)) / 3
# theta_Y = c(sqrt(0.5), 0.5, 1)
# theta_M = c(1, 0.5, 2)
theta_M = c(sqrt(0.5), 0.5, 1) / 3


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
MC_results_delta_MC_delta = pblapply(1:100, function(i) {
    data = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
    # load(paste0("R/Paper MC Study/Datasets/", i, ".RData"))

    tryCatch({

        ## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
        # (fit_Y = suppressMessages(lme4::glmer(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))
        # (fit_M = suppressMessages(lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))
        fit_Y = glmmTMB(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial)
        fit_M = glmmTMB(M ~ X + C1 + C2 + (X | group), data = data, family = binomial)

        # library(broom.mixed)
        # tidy(fit_Y)

        # fit_Y$sdr



        ## Extract model parameter estimates
        info_Y = get_model_pars_TMB(fit_Y)
        info_M = get_model_pars_TMB(fit_M)
        this_Theta_hat = c(unlist(info_Y), unlist(info_M))
        # info_Y = get_model_pars(fit_Y)
        # info_M = get_model_pars(fit_M)
        # this_Theta_hat = c(unlist(info_Y), unlist(info_M))


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
