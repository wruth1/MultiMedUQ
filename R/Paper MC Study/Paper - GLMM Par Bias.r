


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



N = 400
K = 500


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


# # Choose a fitting algorithm for lme4
# library(optimx)
# library(dfoptim)
# Y_info_lme4 = allFit(fit_Y_lme4)
# summary(Y_info_lme4)





# -------------------------- Generate and save data -------------------------- #

# set.seed(1)

# num_reps = 1000

# save_data = pbsapply(1:num_reps, function(i) {
#     data = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
#     save(data, file = paste0("R/Paper MC Study/Data - GLMM Par Bias/", i, ".RData"))
# })





# ---------------------------- Fit models to data ---------------------------- #

# MC_results_delta_MC_delta = pblapply(1:num_reps, function(i) {
all_sim_results = pblapply(1:num_reps, function(i) {
    # data = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
    load(paste0("R/Paper MC Study/Data - GLMM Par Bias/", i, ".RData"))
    tryCatch({


    # ------------------------------- TMB Analysis ------------------------------- #
    fit_Y = glmmTMB(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial)
    fit_M = glmmTMB(M ~ X + C1 + C2 + (X | group), data = data, family = binomial)

    theta_hat_Y = get_model_pars_TMB(fit_Y)
    theta_hat_M = get_model_pars_TMB(fit_M)
    Theta_hat = c(unlist(theta_hat_Y), unlist(theta_hat_M))
    cov_hat = all_pars_cov_mat_TMB(fit_Y, fit_M)

    b_Y = theta_hat_Y[["b"]]
    theta_Y = theta_hat_Y[["theta"]]
    b_M = theta_hat_M[["b"]]
    theta_M = theta_hat_M[["theta"]]
    MEs = all_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)
    cov_MEs = all_covs_MEs_pars(scale, w, cov_hat, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

    # TMB only
    output = list(Theta_hat = Theta_hat, cov_hat = cov_hat, MEs = MEs, cov_MEs = cov_MEs)
    save(output, file = paste0("R/Paper MC Study/Results - GLMM Par Bias/", i, ".RData"))
    return(output)
    }, error = function(e){
      output = NULL
      save(output, file = paste0("R/Paper MC Study/Results - GLMM Par Bias/", i, ".RData"))
      return(output)
    })
}, cl = cl)
# })

stopCluster(cl)


#* Load results and store in a single object
sim_file_names = list.files("R/Paper MC Study/Results - GLMM Par Bias/")
all_sim_results = pblapply(sim_file_names, function(this_name) {
    load(paste0("R/Paper MC Study/Results - GLMM Par Bias/", this_name))
    return(output)
})


#* Extract each type of estimate
## Specifically: Theta hat and ME hat, as well as their SEs, for both TMB and lme4
all_Theta_hats = t(sapply(all_sim_results, function(x) x$Theta_hat))
all_Theta_cov_hats = lapply(all_sim_results, function(x) x$cov_hat)
all_ME_hats = t(sapply(all_sim_results, function(x) x$MEs))
all_ME_cov_hats = lapply(all_sim_results, function(x) x$cov_MEs)


## Empirical covariance of Theta hat
emp_cov_Theta = cov(all_Theta_hats)

## Mean estimated covariances of Theta hat
mean_cov_hat_Theta = Reduce("+", all_Theta_cov_hats) / length(all_Theta_cov_hats)

## Empirical covariance of ME hat
emp_cov_ME = cov(all_ME_hats)

## Mean estimated covariances of ME hat
mean_cov_hat_ME = Reduce("+", all_ME_cov_hats) / length(all_ME_cov_hats)






#* Investigate performance of TMB

## Empirical vs mean estimated covariance of Theta hat
norm(emp_cov_Theta - mean_cov_hat_Theta) / norm(emp_cov_Theta)

### As above, but only variances
norm(as.matrix(diag(emp_cov_Theta) - diag(mean_cov_hat_Theta))) / norm(as.matrix(diag(emp_cov_Theta)))

data.frame(mean = diag(mean_cov_hat_Theta), emp = diag(emp_cov_Theta), diff = diag(emp_cov_Theta) - diag(mean_cov_hat_Theta))



#* Compare estimates with true values

## Theta
### Matrix - vector subtraction matches on the first index (i.e. rows). Use lots of transposes to accommodate this
Theta_hat_errs = t(t(all_Theta_hats) - all_reg_pars)
Theta_hat_err_norms = apply(Theta_hat_errs, 1, function(x) norm(x, type = "2"))
theta_hat_err_rel_norms = Theta_hat_err_norms / norm(all_reg_pars, type = "2")
mean(theta_hat_err_rel_norms)
sd(theta_hat_err_rel_norms)
hist(theta_hat_err_rel_norms)

## Componentwise errors in Theta
abs((all_Theta_hats %>% as.data.frame %>% purrr::map_dbl(mean) - all_reg_pars) / all_reg_pars)




# ---------------------------------------------------------------------------- #
#                             Confidence Intervals                             #
# ---------------------------------------------------------------------------- #


# --------------------------- Supporting Functions --------------------------- #

cov_mat_2_SEs <- function(cov_mat){
    sqrt(diag(cov_mat))
}


# One parameter at a time
build_CIs_one_par <- function(Theta_hats, SE){
    lcls = Theta_hats - qnorm(0.975) * SE
    ucls = Theta_hats + qnorm(0.975) * SE

    return(list(lcl = lcls, ucl = ucls))
}

build_many_CIs <- function(Theta_hats_data, SEs){
    lapply(seq_along(SEs), function(i){
        this_Theta_hats = Theta_hats_data[,i]
        this_SE = SEs[i]

        return(build_CIs_one_par(this_Theta_hats, this_SE))
    })
}

# One dataset at a time
build_one_CI_many_pars <- function(Theta_hats, SEs){
    lapply(seq_along(Theta_hats), function(i){
        this_Theta_hat = Theta_hats[i]
        this_SE = SEs[i]

        return(build_CIs_one_par(this_Theta_hat, this_SE))
    })
}

## SEs can either be a vector of SEs or a list of such vectors with length == nrow(Theta_hats_data)
build_many_CIs_by_dataset <- function(Theta_hats_data, SEs){
    ## Check inputs
    if(is.list(SEs)){
        if(length(SEs) != nrow(Theta_hats_data)){
            stop("Number of estimates does not match number of SE vectors.")
        } else{
            if(length(SEs[[1]]) != ncol(Theta_hats_data)){
                stop("Number of parameters does not match number of SEs.")
            }
        }
        SE_list = TRUE
    } else if(is.numeric(SEs)){
        if(length(SEs) != ncol(Theta_hats_data)){
            stop("Number of parameters does not match number of SEs.")
        }
        SE_list = FALSE
    } else{
        stop("SEs must be either a vector or a list of vectors.")
    }

    ## Build CIs
    lapply(seq_len(nrow(Theta_hats_data)), function(i){
        this_Theta_hats = Theta_hats_data[i,]
        if(SE_list){
            this_SEs = SEs[[i]]
        } else{
            this_SEs = SEs
        }

        return(build_one_CI_many_pars(this_Theta_hats, this_SEs))
    })
}


# Check coverage
coverage_checks_one_par <- function(Theta, CIs){
    some_lcls = CIs$lcl
    some_ucls = CIs$ucl

    some_checks = some_lcls < Theta & Theta < some_ucls

    return(some_checks)
}

many_coverage_checks <- function(Theta, CIs){
    lapply(seq_along(CIs), function(i){
        this_CIs = CIs[[i]]
        this_Theta = Theta[i]

        return(coverage_checks_one_par(this_Theta, this_CIs))
    })
}

coverage_checks_2_rates <- function(coverage_checks_list){
    sapply(coverage_checks_list, function(x) mean(x))
}


get_coverage_rates <- function(Theta_hats_data, cov_mat, true_Thetas){
    SEs = cov_mat_2_SEs(cov_mat)
    CIs = build_many_CIs(Theta_hats_data, SEs)
    coverage_checks = many_coverage_checks(true_Thetas, CIs)

    return(coverage_checks_2_rates(coverage_checks))
}

get_coverage_rates_many_cov_mats <- function(Theta_hats_data, cov_mat_list, true_Thetas){
    SEs_list = lapply(cov_mat_list, cov_mat_2_SEs)

    # Nesting structure is [dataset][parameter][lcl/ucl]
    # Needs to be [parameter][lcl/ucl][dataset]
    CIs_raw = build_many_CIs_by_dataset(Theta_hats_data, SEs_list)

    CIs = lapply(seq_len(ncol(Theta_hats_data)), function(i){
        this_lcls = sapply(seq_len(nrow(Theta_hats_data)), function(j) CIs_raw[[j]][[i]]$lcl) %>% unname()
        this_ucls = sapply(seq_len(nrow(Theta_hats_data)), function(j) CIs_raw[[j]][[i]]$ucl) %>% unname()

        return(list(lcl = this_lcls, ucl = this_ucls))
    })

    coverage_checks = many_coverage_checks(true_Thetas, CIs)
    coverage_rates = coverage_checks_2_rates(coverage_checks)

    return(coverage_rates)
        
}






# ------------------------------ GLMM Parameters ----------------------------- #

cover_emp = get_coverage_rates(all_Theta_hats, emp_cov_Theta, all_reg_pars)
cover_hat = get_coverage_rates(all_Theta_hats, mean_cov_hat_Theta, all_reg_pars)

## CI for each dataset based on estimated SEs, not mean estimated SEs
cov_hat_honest = get_coverage_rates_many_cov_mats(all_Theta_hats, all_Theta_cov_hats, all_reg_pars)


data_cover = data.frame(emp = cover_emp, hat = cover_hat, honest = cov_hat_honest)


cover_emp - cover_hat
cover_hat - cov_hat_honest


# ----------------------------- Mediation Effects ---------------------------- #
true_MEs = all_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M, which_REs = which_REs)


cover_emp_ME = get_coverage_rates(all_ME_hats, emp_cov_ME, true_MEs)

cover_hat_ME = get_coverage_rates(all_ME_hats, mean_cov_hat_ME, true_MEs)

cov_hat_ME_honest = get_coverage_rates_many_cov_mats(all_ME_hats, all_ME_cov_hats, true_MEs)

data_cover_ME = data.frame(emp = cover_emp_ME, hat = cover_hat_ME, honest = cov_hat_ME_honest)







# ---------------------------------------------------------------------------- #
#                                     Bias                                     #
# ---------------------------------------------------------------------------- #

# ------------------------------- Compute bias ------------------------------- #

mean_Theta_hat = colMeans(all_Theta_hats)
mean_ME_hat = colMeans(all_ME_hats)


bias_Theta = mean_Theta_hat - all_reg_pars
bias_ME = mean_ME_hat - true_MEs


rel_bias_Theta = abs(bias_Theta / all_reg_pars)
rel_bias_ME = abs(bias_ME / true_MEs)
