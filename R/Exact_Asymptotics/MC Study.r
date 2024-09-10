



library(lme4)
library(merDeriv)
library(tictoc)
library(pbapply)
library(parallel)
library(magrittr)
library(dplyr)
library(kableExtra)
source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
devtools::load_all()





# Set parameters and fit models

# N = 100
# N = 20
N = 200

all_Ks = c(50, 100, 200, 400, 800)
# all_Ks = c(50, 100, 200)
# all_Ks = c(50, 100)
# K = 50

# num_reps = 50
num_reps = 500
# num_reps = 20


which_REs = c("Y.Int", "Y.X", "M.All")


x = 0
x_m = 1

b_Y = c(0,0.5,1,0,0) * 2
# theta_Y = c(sqrt(0.5), 0.5, 0, 1, 0.5, sqrt(0.5))
theta_Y = c(sqrt(0.5), 0, 1)

b_M = c(0,1,0,0) * 2
theta_M = c(1, 0.5, 2)

w = c(0,0)



list_par_hats = list()
list_par_cov_hats = list()



set.seed(11111111)


# Setup cluster
# cl = makeCluster(detectCores() - 2)
cl = makeCluster(10)
# cl = makeCluster(5)
clusterExport(cl, c("N", "b_Y", "theta_Y", "b_M", "theta_M", "which_REs"))
clusterEvalQ(cl, {
    library(lme4)
    library(merDeriv)
    library(tictoc)
    library(pbapply)
    # library(parallel)
    source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
    devtools::load_all()
})





# Run MC study
for(j in seq_along(all_Ks)){

    K = all_Ks[j]
    if(nrow(showConnections()) > 0) clusterExport(cl, "K") # Only run if there is an active cluster



    # all_med_hats = c()
    # all_cov_hats = list()
    
    
    print(paste0("K = ", K, "; number ", j, " of ", length(all_Ks)))

    some_info_par = pblapply(1:num_reps, function(i) {
        data = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)

        w = c(0,0)


        # Singularities can arise in fitted covariance matrices. Use tryCatch to skip these cases.
        tryCatch({
            

        ## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
        (fit_Y = suppressMessages(lme4::glmer(Y ~ X + M + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))
        (fit_M = suppressMessages(lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))

        # Explore covariances


        ## Model parameters
        info_Y = get_model_pars(fit_Y)
        info_M = get_model_pars(fit_M)

        par_hat = c(unlist(info_Y), unlist(info_M))
        par_cov_hat = all_pars_cov_mat(fit_Y, fit_M)


        this_output = list(par_hat = par_hat, par_cov_hat = par_cov_hat)
        }, error = function(e){
            this_output = NULL
        })

        return(this_output)
    }, cl=cl)
    # })


    some_par_hats = t(sapply(some_info_par, function(x) x$par_hat))
    some_par_cov_hats = lapply(some_info_par, function(x) x$par_cov_hat)


    list_par_hats[[j]] = some_par_hats
    list_par_cov_hats[[j]] = some_par_cov_hats
}

stopCluster(cl)


# Explore timings
# # times_15 = c(8,14,34,50,107) # N = 20
# times_15 = c(9,15,29,56,109) # N = 200
# # times_10 = c(4, 12, 23, 38, 88) # N = 20
# times_10 = c(6,12,21,33,95) # N = 200
# times_5 = c(2, 4, 12, 21, 31)

# times_15 / times_5
# times_10 / times_5
# times_15 / times_10

# per_iteration = c(sum(times_15)/15, sum(times_10)/10, sum(times_5)/5)
# per_iteration = c(sum(times_15)/15, sum(times_10)/10)

# sum(times_15) / sum(times_5)
# sum(times_10) / sum(times_5)
# sum(times_15) / sum(times_10)





# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC.RData")
# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC-Large_K.RData")
# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC-Large_K_2.RData")
# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC-Large_K_Many_Reps.RData")
# load("Par_Hat_MC.RData", verbose = TRUE)
# load("Par_Hat_MC-Large_K.RData", verbose = TRUE)
# load("Par_Hat_MC-Large_K_2.RData", verbose = TRUE)
load("Par_Hat_MC-Large_K_Many_Reps.RData", verbose = TRUE)


# Concatenate parameter estimates from multiple MC studies


all_files = c("Par_Hat_MC-Large_K.RData", "Par_Hat_MC-Large_K_2.RData", "Par_Hat_MC-Large_K_Many_Reps.RData")

load(all_files[1], verbose = TRUE)

large_list_par_hats = list_par_hats
large_list_par_cov_hats = list_par_cov_hats

for(file in all_files[2:3]){

    load(file, verbose = TRUE)

    for(i in seq_along(all_Ks)){

        large_list_par_hats[[i]] = rbind(large_list_par_hats[[i]], list_par_hats[[i]])

        large_list_par_cov_hats[[i]] = c(large_list_par_cov_hats[[i]], list_par_cov_hats[[i]])
    }
}

list_par_hats = large_list_par_hats
list_par_cov_hats = large_list_par_cov_hats



# Extract empirical and average estimated covariances for each value of K

num_pars = ncol(list_par_hats[[1]])


## Remove any runs with NA parameter estimates
for(j in seq_along(all_Ks)){
    # Covariance matrices first to preserve information in parameter estimates
    list_par_cov_hats[[j]] = list_par_cov_hats[[j]][complete.cases(list_par_hats[[j]])]
    list_par_hats[[j]] = na.omit(list_par_hats[[j]])
}

lengths(list_par_cov_hats)
sapply(list_par_hats, nrow)


list_emp_covs = list()
list_mean_covs = list()
list_all_errs = list()

# Extract empirical covariance and mean estimated covariance.
for(j in seq_along(all_Ks)){

    this_emp_cov = cov(list_par_hats[[j]])
    list_emp_covs[[j]] = this_emp_cov

    some_cov_hats = list_par_cov_hats[[j]]
    this_mean_cov = Reduce("+", some_cov_hats) / length(some_cov_hats)

    list_mean_covs[[j]] = this_mean_cov

    list_all_errs[[j]] = lapply(some_cov_hats, function(x) x - this_emp_cov)

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



# Summary statistics of absolute and relative errors in estimating variances (i.e. diagonal elements of the cov mat)

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
    norm_diffs = norm(this_diffs, type = "2")

    mean_rel_diffs = mean(this_rel_diffs)
    SD_rel_diffs = sd(this_rel_diffs)
    norm_rel_diffs = norm(this_rel_diffs, type = "2")

    this_info = c(all_Ks[i], mean_diffs, SD_diffs, norm_diffs, mean_rel_diffs, SD_rel_diffs, norm_rel_diffs)

    info_rel_diffs = rbind(info_rel_diffs, this_info)
}
colnames(info_rel_diffs) = c("K", "Mean-Diff", "SD-Diff", "Norm-Diff", "Mean-Rel", "SD-Rel", "Norm-Rel")

info_rel_diffs




# Absolute and relative matrix norms of empirical vs estimated covariances

info_rel_norms = data.frame()
for(i in seq_along(all_Ks)){

    # Compare empirical and estimated covariances
    this_emp_mat = list_emp_covs[[i]]
    this_mean_mat = list_mean_covs[[i]]

    # this_norm_diff = norm(this_emp_mat - this_mean_mat, type = "2")
    # this_norm_rel = this_norm_diff / norm(this_emp_mat, type = "2")

    this_norm_diff = norm(this_emp_mat - this_mean_mat, type = "F")
    this_norm_rel = this_norm_diff / norm(this_emp_mat, type = "F")

    this_info = c(all_Ks[i], this_norm_diff, this_norm_rel)

    info_rel_norms = rbind(info_rel_norms, this_info)
}
colnames(info_rel_norms) = c("K", "Norm_Diff", "Norm_Rel")

info_rel_norms
scaled_rel_norms = info_rel_norms %>% mutate(scale_abs_norm = Norm_Diff * K, scale_rel_norm = Norm_Rel * sqrt(K))
scaled_rel_norms %>% 
    dplyr::mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
    kbl(., format = "latex", booktabs = T, 
        caption = "Absolute and relative error between the mean estimated covariance and the (Monte Carlo) empirical covariance matrix of the estimated parameter values, as well as those same values multiplied by $K$.",
        label = "tab:err_scaled_ENC",
        col.names = c("K", "Abs", "Rel", "Scaled-Abs", "Scaled-Rel"),
        align = "c"
    )


## Mean error vs error of mean in estimating empirical covariance
info_all_err_norms = data.frame()
for(i in seq_along(all_Ks)){
    this_all_errs = list_all_errs[[i]]
    this_all_err_norms = sapply(this_all_errs, norm, type = "2")
    this_mean_err_norm = mean(this_all_err_norms)
    this_SD_err_norm = sd(this_all_err_norms)

    this_emp_mat = list_emp_covs[[i]]
    this_mean_mat = list_mean_covs[[i]]
    this_err_mean = norm(this_emp_mat - this_mean_mat, type = "2")

    this_info = c(all_Ks[i], this_mean_err_norm, this_SD_err_norm, this_err_mean)
    info_all_err_norms = rbind(info_all_err_norms, this_info)
}
colnames(info_all_err_norms) = c("K", "Mean_Error", "SD_Error", "Error_of_Mean")

info_all_err_norms

info_all_err_norms_scaled = info_all_err_norms %>% 
    mutate(Mean_Error = Mean_Error * K, Error_of_Mean = Error_of_Mean * K) %>%
    select(-SD_Error)
info_all_err_norms_scaled


# Absolute and relative difference for each parameter at each level of K
data_covs = data.frame()
for(i in seq_along(all_Ks)){

    # Extract empirical and estimated variances
    this_emp_mat = list_emp_covs[[i]]
    this_mean_mat = list_mean_covs[[i]]

    this_emp_vars = diag(this_emp_mat)
    # this_mean_vars = diag(this_mean_mat)
    this_mean_vars = diag(this_mean_mat)

    this_diff = abs(this_emp_vars - this_mean_vars)
    this_rel_diff = this_diff / this_emp_vars

    this_info = cbind(all_Ks[i], this_emp_vars, this_mean_vars, this_diff, this_rel_diff)

    data_covs = rbind(data_covs, this_info)
}
colnames(data_covs) = c("K", "Emp", "Mean", "Diff", "Rel-Diff")

data_covs





# Covariances of ENCs

list_ENC_hats = list()
list_ENC_cov_hats = list()

for(i in seq_along(all_Ks)){

    some_ENC_hats = data.frame()
    some_ENC_cov_hats = list()

    for(j in seq_len(num_reps)){

        this_par_hat = list_par_hats[[i]][j,]

        this_b_Y = this_par_hat[1:5]
        this_theta_Y = this_par_hat[6:8]
        this_b_M = this_par_hat[9:12]
        this_theta_M = this_par_hat[13:15]


        this_ENC_hat = all_ENCs(w, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
        some_ENC_hats = rbind(some_ENC_hats, this_ENC_hat)


        this_Sigma = list_par_cov_hats[[i]][[j]]
        this_ENC_cov_hat = all_covs_ENC_pars(w, this_Sigma, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
        some_ENC_cov_hats[[j]] = this_ENC_cov_hat

    }

    colnames(some_ENC_hats) = c("11", "10", "01", "00")
    list_ENC_hats[[i]] = some_ENC_hats
    list_ENC_cov_hats[[i]] = some_ENC_cov_hats
}


## Summarize covariance estimates

### Empirical

list_ENC_emp_covs = list()

for(i in seq_along(all_Ks)){
    list_ENC_emp_covs[[i]] = cov(list_ENC_hats[[i]])
}


### Estimated

list_ENC_mean_covs = list()

for(i in seq_along(all_Ks)){
    some_ENC_cov_hats = list_ENC_cov_hats[[i]]
    this_mean_ENC_cov_hat = Reduce("+", some_ENC_cov_hats) / num_reps
    list_ENC_mean_covs[[i]] = this_mean_ENC_cov_hat
}



## Compare covariance estimates

### Easier: Variances (i.e. diagonal elements)
info_var_diffs = data.frame()
for(i in seq_along(all_Ks)){
    some_var_diffs = diag(list_ENC_emp_covs[[i]] - list_ENC_mean_covs[[i]])
    some_rel_var_diffs = some_var_diffs / diag(list_ENC_emp_covs[[i]])

    this_info = c(all_Ks[i], mean(some_var_diffs), sd(some_var_diffs), mean(some_rel_var_diffs), sd(some_rel_var_diffs))

    info_var_diffs = rbind(info_var_diffs, this_info)
}
colnames(info_var_diffs) = c("K", "Mean_Diff", "SD_Diff", "Mean_Rel_Diff", "SD_Rel_Diff")

info_var_diffs

info_var_diffs %>% mutate(scale_abs_diff = Mean_Diff * K) %>% select(-"Mean_Rel_Diff", -"SD_Rel_Diff")


### Harder: Matrix norms
info_ENC_norms = data.frame()
for(i in seq_along(all_Ks)){
    diff_mat = list_ENC_emp_covs[[i]] - list_ENC_mean_covs[[i]]
    # diff_mat_norm = norm(diff_mat, type = "2")
    diff_mat_norm = norm(diff_mat, type = "F")

    # rel_norm = diff_mat_norm / norm(list_ENC_emp_covs[[i]], type = "2")
    rel_norm = diff_mat_norm / norm(list_ENC_emp_covs[[i]], type = "F")

    this_info = c(all_Ks[i], diff_mat_norm, rel_norm)
    info_ENC_norms = rbind(info_ENC_norms, this_info)
}
colnames(info_ENC_norms) = c("K", "Norm_Diff", "Rel_Norm_Diff")

info_ENC_norms

scaled_ENC_norms = info_ENC_norms %>% mutate(scale_abs_diff = Norm_Diff * K) %>% select(-"Rel_Norm_Diff")




# Covariances of mediation effects


list_ME_hats = list()
list_ME_cov_hats = list()

this_scale = c("diff", "rat", "OR")

for(i in seq_along(all_Ks)){
    print(paste0("K = ", all_Ks[i], "; number ", i, " of ", length(all_Ks)))

    some_ME_hats = data.frame()
    some_ME_cov_hats = list()

    for(j in seq_len(num_reps)){

        this_par_hat = list_par_hats[[i]][j,]

        this_b_Y = this_par_hat[1:5]
        this_theta_Y = this_par_hat[6:8]
        this_b_M = this_par_hat[9:12]
        this_theta_M = this_par_hat[13:15]


        this_ME_hat = all_MEs_pars(this_scale, w, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
        some_ME_hats = rbind(some_ME_hats, this_ME_hat)


        this_Sigma = list_par_cov_hats[[i]][[j]]
        this_ME_cov_hat = all_covs_MEs_pars(this_scale, w, this_Sigma, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
        some_ME_cov_hats[[j]] = this_ME_cov_hat

    }

    colnames(some_ME_hats) = c("11", "10", "01", "00")
    list_ME_hats[[i]] = some_ME_hats
    list_ME_cov_hats[[i]] = some_ME_cov_hats
}


## Summarize covariance estimates

### Empirical

list_ME_emp_covs = list()

for(i in seq_along(all_Ks)){
    list_ME_emp_covs[[i]] = cov(list_ME_hats[[i]])
}


### Estimated

list_ME_mean_covs = list()

for(i in seq_along(all_Ks)){
    some_ME_cov_hats = list_ME_cov_hats[[i]]
    this_mean_ME_cov_hat = Reduce("+", some_ME_cov_hats) / num_reps
    list_ME_mean_covs[[i]] = this_mean_ME_cov_hat
}



## Compare covariance estimates

### Easier: Variances (i.e. diagonal elements)
info_var_diffs = data.frame()
for(i in seq_along(all_Ks)){
    some_var_diffs = diag(list_ME_emp_covs[[i]] - list_ME_mean_covs[[i]])
    some_rel_var_diffs = some_var_diffs / diag(list_ME_emp_covs[[i]])

    this_info = c(all_Ks[i], mean(some_var_diffs), sd(some_var_diffs), mean(some_rel_var_diffs), sd(some_rel_var_diffs))

    info_var_diffs = rbind(info_var_diffs, this_info)
}
colnames(info_var_diffs) = c("K", "Mean_Diff", "SD_Diff", "Mean_Rel_Diff", "SD_Rel_Diff")

info_var_diffs

info_var_diffs %>% mutate(scale_abs_diff = Mean_Diff * K) %>% select(-"Mean_Rel_Diff", -"SD_Rel_Diff")


### Harder: Matrix norms
info_ME_norms = data.frame()
for(i in seq_along(all_Ks)){
    diff_mat = list_ME_emp_covs[[i]] - list_ME_mean_covs[[i]]
    diff_mat_norm = norm(diff_mat, type = "2")

    # rel_norm = diff_mat_norm / norm(list_ME_emp_covs[[i]], type = "2")
    rel_norm = diff_mat_norm / norm(list_ME_emp_covs[[i]], type = "F")

    this_info = c(all_Ks[i], diff_mat_norm, rel_norm)
    info_ME_norms = rbind(info_ME_norms, this_info)
}
colnames(info_ME_norms) = c("K", "Norm_Diff", "Rel_Norm_Diff")

info_ME_norms

scaled_ME_norms = info_ME_norms %>% mutate(scale_abs_diff = Norm_Diff * K) %>% select(-"Rel_Norm_Diff")



# Summarize errors for all covariance estimates

theta_mat = scaled_rel_norms %>% select(-Norm_Rel, -scale_rel_norm)
ENC_mat = scaled_ENC_norms
ME_mat = scaled_ME_norms

colnames(theta_mat) = c("K", "theta-abs", "theta-scaled")
colnames(ENC_mat) = c("K", "ENC-abs", "ENC-scaled")
colnames(ME_mat) = c("K", "ME-abs", "ME-scaled")

summary_mat = cbind(theta_mat, ENC_mat[,-1], ME_mat[,-1])

summary_mat %>% 
    dplyr::mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
    kbl(., format = "latex", booktabs = F, 
        caption = "Raw and scaled absolute errors in estimating empirical covariance matrices of estimators for the parameters, ENCs and mediation effects.",
        label = "tab:err_summary",
        col.names = c("K", "Abs $\\theta$", "Scaled $\\theta$", "Abs ENC", "Scaled ENC", "Abs ME", "Scaled ME"),
        align = "c", escape = F
    )



# Decompose cov mat into components

cov_decomp = function(Sigma){
    b_Y_comp = Sigma[1:5, 1:5]
    theta_Y_comp = Sigma[6:8, 6:8]
    cross_Y_comp = Sigma[1:5, 6:8]

    b_M_comp = Sigma[9:12, 9:12]
    theta_M_comp = Sigma[13:15, 13:15]
    cross_M_comp = Sigma[9:12, 13:15]

    output = list(
        b_Y = b_Y_comp,
        theta_Y = theta_Y_comp,
        cross_Y = cross_Y_comp,
        b_M = b_M_comp,
        theta_M = theta_M_comp,
        cross_M = cross_M_comp
    )

    return(output)
}


list_emp_decomps = lapply(list_emp_covs, cov_decomp)
list_mean_decomps = lapply(list_mean_covs, cov_decomp)

data_decomp = data.frame()

for(i in seq_along(all_Ks)){

    this_emp_decomp = list_emp_decomps[[i]]
    this_mean_decomp = list_mean_decomps[[i]]

    for(j in seq_along(this_emp_decomp)){

        abs_diff = norm(this_emp_decomp[[j]] - this_mean_decomp[[j]], type = "2")
        rel_diff = abs_diff / norm(this_emp_decomp[[j]], type = "2")


        this_info = c(all_Ks[i], names(this_emp_decomp)[j], abs_diff, rel_diff)
        data_decomp = rbind(data_decomp, this_info)
    }
}
colnames(data_decomp) = c("K", "Component", "Abs_Diff", "Rel_Diff")
data_decomp %<>% mutate(across(c(1,3,4), as.numeric))

data_decomp

data_decomp %>% 
    arrange(Component)



flavour_decomp = data_decomp %>% 
    mutate(diag = ifelse(Component %in% c("cross_Y", "cross_M"), F, T)) %>% 
    group_by(K, diag) %>% 
    summarize(mean_abs = mean(Abs_Diff), mean_rel = mean(Rel_Diff), .groups = "drop")

flavour_decomp %>%
    arrange(diag)






# Trajectory of errors across increasing MC sizes

list_emp_covs
list_all_errs
# Extract empirical covariance and mean estimated covariance.
for(j in seq_along(all_Ks)){

    this_emp_cov = cov(list_par_hats[[j]])
    list_emp_covs[[j]] = this_emp_cov

    some_cov_hats = list_par_cov_hats[[j]]
    this_mean_cov = Reduce("+", some_cov_hats) / length(some_cov_hats)

    list_mean_covs[[j]] = this_mean_cov

    list_all_errs[[j]] = lapply(some_cov_hats, function(x) x - this_emp_cov)

}



norms_by_MC_size = data.frame()

total_num_reps = nrow(list_par_hats[[1]])
# all_Bs = seq(100, 200, by=10)
all_Bs = seq(100, total_num_reps, by=20)

counter = 0

for(i in seq_along(all_Ks)){
    for(r in seq_along(all_Bs)){
        this_B = all_Bs[r]

        # this_emp_cov = cov(list_par_hats2[[i]][1:this_B,])
        # this_mean_cov = Reduce("+", list_par_cov_hats2[[i]][1:this_B]) / this_B

        # # Pool both simulations' results
        this_emp_cov = cov(large_list_par_hats[[i]][1:this_B,])
        this_mean_cov = Reduce("+", large_list_par_cov_hats[[i]][1:this_B]) / this_B

        this_err = norm(this_emp_cov - this_mean_cov, type = "2")

        this_err_scaled = this_err * all_Ks[i]
        this_err_extra_scaled = this_err_scaled * sqrt(all_Ks[i])

        this_info = c(all_Ks[i], this_B, this_err, this_err_scaled, this_err_extra_scaled)


        counter = counter + 1
        norms_by_MC_size = rbind(norms_by_MC_size, this_info)
    }
}
colnames(norms_by_MC_size) = c("K", "MC_Size", "Abs_Error", "Scaled_Error", "Extra_Scaled_Error")

norms_by_MC_size

# Sort norms_by_MC size by K
norms_by_MC_size %>% arrange(K)


library(ggplot2)


norms_by_MC_size %>% 
    ggplot(aes(x = MC_Size, y = Scaled_Error, color = as.factor(K))) +
    geom_line() +
    geom_point() +
    theme_bw()

# Grid of plots
# Absolute error scaled by K
pdf("Plots/Abs_Error_by_MC_Size.pdf", width = 10, height = 10)
norms_by_MC_size %>%
    ggplot(aes(x = MC_Size, y = Abs_Error, color = as.factor(K))) +
    geom_line() +
    geom_point() +
    facet_wrap(~K, scales="free") +
    theme_bw()
dev.off()

# ## Both scales on same axes. Redundant, since scale is constant in B
# scale_factor = 15
# norms_by_MC_size %>%
#     ggplot(aes(x = MC_Size)) +
#     geom_line(aes(y = Scaled_Error), color = "blue") +
#     geom_line(aes(y = Extra_Scaled_Error / scale_factor), color = "red") +
#     facet_wrap(~K) +
#     theme_bw() +
#     # legend for colors
#     theme(legend.position = "right") +
#     scale_color_manual(values = c("blue", "red")) +
#     scale_y_continuous(name = "K Times Absolute Error", sec.axis = sec_axis(trans =~./scale_factor, name = "K^3/2 Times Absolute Error"))


# Plot scaled error vs K for maximum MC size
norms_by_MC_size %>%
    filter(MC_Size == max(MC_Size)) %>%
    ggplot(aes(x = K, y = Abs_Error)) +
    geom_line() +
    geom_point() +
    theme_bw()


# Estimate rate at which errors scale with K
best_err_estimates = norms_by_MC_size %>%
    filter(MC_Size == max(MC_Size)) %>%
    mutate(log_K = log(K), log_err = log(Abs_Error))

fit_err_rate = lm(log_err ~ log_K, data = best_err_estimates[-1,])
summary(fit_err_rate)

plot(fit_err_rate, 1)


pdf("Plots/log_Abs_Err_vs_log_K.pdf", width = 10, height = 10)
norms_by_MC_size %>%
    filter(MC_Size == max(MC_Size)) %>%
    ggplot(aes(x = log(K), y = log(Abs_Error))) +
    geom_line() +
    geom_point() +
    theme_bw()
dev.off()

filter(norms_by_MC_size, MC_Size == max(MC_Size))
