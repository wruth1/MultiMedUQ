




# TODO: Figure out why there is a discrepancy between the current implementation and the one from the Exact_Asymptotics project

#* To run the following:
library(lme4)
library(merDeriv)
library(tictoc)
source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
devtools::load_all()


# Set parameters and fit models

N = 100

all_Ks = c(50, 100, 200, 400)
# K = 50

num_reps = 50




x = 0
x_m = 1

b_Y = c(0,0,1,0,0)
# theta_Y = c(sqrt(0.5), 0.5, 0, 1, 0.5, sqrt(0.5))
theta_Y = c(sqrt(0.5), 0, 1)

b_M = c(0,0,0,0)
theta_M = c(1, 0.5, 2)



list_par_hats = list()
list_par_cov_hats = list()



set.seed(1)

tic()

for(j in seq_along(all_Ks)){
    K = all_Ks[j]

    all_par_hats = c()
    all_par_cov_hats = list()



    # all_med_hats = c()
    # all_cov_hats = list()



    for(i in 1:num_reps){
        # print(paste0("Rep ", i, " of ", num_reps))
        print(paste0("Rep ", i, " of ", num_reps, ", K number ", j, " of ", length(all_Ks)))

        data = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F)   #!!!!!!!!! Does this get the REs right?
                                                                                         #! No!!! #! Start here!!!
                                                                                         #! Needs to pass which_REs

        w = c(0,0)

        ## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
        (fit_Y = lme4::glmer(Y ~ X + M + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))
        (fit_M = lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))

        # Explore covariances


        ## Model parameters
        info_Y = get_model_pars(fit_Y)
        info_M = get_model_pars(fit_M)

        par_hat = c(unlist(info_Y), unlist(info_M))
        par_cov_hat = all_pars_cov_mat(fit_Y, fit_M)


        all_par_hats = rbind(all_par_hats, par_hat)
        all_par_cov_hats[[i]] = par_cov_hat




        ## Mediation effects
        # med_hats = all_MEs_models(scale = "OR", w, fit_Y, fit_M, which_REs = c("Y.Int", "Y.X", "M.All")) 

        # cov_hat = all_cov_MEs(scale = "OR", w, fit_Y, fit_M, which_REs = c("Y.Int", "Y.X", "M.All"))



        # all_med_hats = rbind(all_med_hats, med_hats)
        # all_cov_hats[[i]] = cov_hat

    }

    list_par_hats[[j]] = all_par_hats
    list_par_cov_hats[[j]] = all_par_cov_hats
}


toc()



# Extract empirical and average estimated covariances for each value of K

num_pars = ncol(list_par_hats[[1]])

list_emp_covs = list()
list_mean_covs = list()


for(j in seq_along(all_Ks)){

    list_emp_covs[[j]] = cov(list_par_hats[[j]])

    this_sum_covs = matrix(0, num_pars, num_pars)
    for(i in 1:num_reps){
        this_sum_covs = this_sum_covs + list_par_cov_hats[[j]][[i]]
    }
    list_mean_covs[[j]] = this_sum_covs / num_reps
}



# Compute summary statistics of ratios of diagonal elements across pairs of K-levels

## Empirical variances
info_rats = data.frame()
for(i in (1:(length(all_Ks)-1))){
    for(j in (i+1):length(all_Ks)){
        K1 = all_Ks[i]
        K2 = all_Ks[j]
        
        this_rat_mat = list_emp_covs[[i]] / list_emp_covs[[j]]

        this_rats = diag(this_rat_mat)

        mean_rats = mean(this_rats)
        SD_rats = sd(this_rats)

        this_info = c(K1, K2, mean_rats, SD_rats, K2/K1, SD_rats / mean_rats)

        info_rats = rbind(info_rats, this_info)
    }
}
colnames(info_rats) = c("K1", "K2", "Mean", "SD", "Expected", "COV")

info_rats

## Mean estimated variances
info_rats_mean = data.frame()
for(i in (1:(length(all_Ks)-1))){
    for(j in (i+1):length(all_Ks)){
        K1 = all_Ks[i]
        K2 = all_Ks[j]
        
        this_rat_mat = list_mean_covs[[i]] / list_mean_covs[[j]]

        this_rats = diag(this_rat_mat)

        mean_rats = mean(this_rats)
        SD_rats = sd(this_rats)

        this_info = c(K1, K2, mean_rats, SD_rats, K2/K1, SD_rats / mean_rats)

        info_rats_mean = rbind(info_rats_mean, this_info)
    }
}
colnames(info_rats_mean) = c("K1", "K2", "Mean", "SD", "Expected", "COV")

info_rats_mean



info_rel_diffs = data.frame()
for(i in seq_along(all_Ks)){

    # Compare empirical and estimated covariances
    this_emp_mat = list_emp_covs[[i]]
    this_mean_mat = list_mean_covs[[i]]

    this_diff_mat = this_emp_mat - this_mean_mat
    this_rel_diff_mat = this_diff_mat / this_emp_mat

    this_diffs = diag(this_diff_mat)
    this_rel_diffs = diag(this_rel_diff_mat)

    mean_diffs = mean(this_diffs)
    SD_diffs = sd(this_diffs)

    mean_rel_diffs = mean(this_rel_diffs)
    SD_rel_diffs = sd(this_rel_diffs)

    this_info = c(all_Ks[i], mean_diffs, SD_diffs, mean_rel_diffs, SD_rel_diffs)

    info_rel_diffs = rbind(info_rel_diffs, this_info)
}
colnames(info_rel_diffs) = c("K", "Mean-Diff", "SD-Diff", "Mean-Rel", "SD-Rel")

info_rel_diffs



data_covs = data.frame()
for(i in seq_along(all_Ks)){

    # Extract empirical and estimated variances
    this_emp_mat = list_emp_covs[[i]]
    this_mean_mat = list_mean_covs[[i]]

    this_emp_vars = diag(this_emp_mat)
    # this_mean_vars = diag(this_mean_mat)
    this_mean_vars = all_Ks[i] * diag(this_mean_mat)

    this_diff = abs(this_emp_vars - this_mean_vars)
    this_rel_diff = this_diff / this_emp_vars

    this_info = cbind(all_Ks[i], this_emp_vars, this_mean_vars, this_diff, this_rel_diff)

    data_covs = rbind(data_covs, this_info)
}
colnames(data_covs) = c("K", "Emp", "Mean", "Diff", "Rel-Diff")

data_covs
