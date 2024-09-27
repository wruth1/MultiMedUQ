



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
devtools::load_all()




# Set parameters and fit models

# N = 100
# N = 20
# N = 40
# N = 60
N = 200

all_Ks = c(50, 100, 200, 400, 800)
# all_Ks = c(50, 100, 200)
# all_Ks = 50 * (2:6)
# K = 50

# num_reps = 30
num_reps = 500
# num_reps = 1200


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



set.seed(55555555)


# Setup cluster
# cl = makeCluster(detectCores() - 2)
cl = makeCluster(15)
# cl = makeCluster(10)
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
# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC-Large_K_3.RData")
# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC-Large_K_4.RData")

# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC-Large_K_5.RData")
# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC-Large_K_6.RData")

# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC-Large_K_Many_Reps.RData")

# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC-Large_K_Pooled.RData")


# load("Par_Hat_MC.RData", verbose = TRUE)
# load("Par_Hat_MC-Large_K.RData", verbose = TRUE)
# load("Par_Hat_MC-Large_K_2.RData", verbose = TRUE)
# load("Par_Hat_MC-Large_K_3.RData", verbose = TRUE)
# load("Par_Hat_MC-Large_K_Many_Reps.RData", verbose = TRUE)
load("Par_Hat_MC-Large_K_Pooled.RData", verbose = TRUE)


# # # Concatenate parameter estimates from multiple MC studies
# # ## Alternatively, load the Pooled version of the saved data.
# all_files = c("Par_Hat_MC-Large_K_5.RData", "Par_Hat_MC-Large_K_6.RData", "Par_Hat_MC-Large_K_Pooled.RData")

# load(all_files[1], verbose = TRUE)

# large_list_par_hats = list_par_hats
# large_list_par_cov_hats = list_par_cov_hats

# for(file in all_files[2:length(all_files)]){

#     load(file, verbose = TRUE)

#     for(i in seq_along(all_Ks)){

#         large_list_par_hats[[i]] = rbind(large_list_par_hats[[i]], list_par_hats[[i]])

#         large_list_par_cov_hats[[i]] = c(large_list_par_cov_hats[[i]], list_par_cov_hats[[i]])
#     }
# }

# list_par_hats = large_list_par_hats
# list_par_cov_hats = large_list_par_cov_hats

# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC-Large_K_Pooled.RData")







#* Extract empirical and average estimated covariances for each value of K

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



#* Compute empirical and analytical covariance estimates separately on each half of the data and compare.

# Merge multiple datasets
all_files = c("Par_Hat_MC-Large_K_4.RData", "Par_Hat_MC-Large_K_6.RData")

load(all_files[1], verbose = TRUE)



large_list_par_hats = list_par_hats
large_list_par_cov_hats = list_par_cov_hats

for(file in all_files[2:length(all_files)]){

    load(file, verbose = TRUE)

    for(i in seq_along(all_Ks)){

        large_list_par_hats[[i]] = rbind(large_list_par_hats[[i]], list_par_hats[[i]])

        large_list_par_cov_hats[[i]] = c(large_list_par_cov_hats[[i]], list_par_cov_hats[[i]])
    }
}

list_par_hats = large_list_par_hats
list_par_cov_hats = large_list_par_cov_hats

length(list_par_cov_hats[[1]])

# Compute LHS and RHS empirical and analytical covariances

lhs_emp_covs = list()
rhs_emp_covs = list()
lhs_mean_covs = list()
rhs_mean_covs = list()

for(j in seq_along(all_Ks)){
    this_K = all_Ks[j]
    # Empirical covariances
    some_par_hats = list_par_hats[[j]]
    lhs_par_hats = some_par_hats[1:(nrow(some_par_hats)/2),]
    rhs_par_hats = some_par_hats[(nrow(some_par_hats)/2+1):nrow(some_par_hats),]

    lhs_emp_cov = cov(lhs_par_hats)
    rhs_emp_cov = cov(rhs_par_hats)

    lhs_emp_covs[[j]] = lhs_emp_cov
    rhs_emp_covs[[j]] = rhs_emp_cov


    # Analytical covariances
    some_cov_hats = list_par_cov_hats[[j]]
    lhs_cov_hats = some_cov_hats[1:(length(some_cov_hats)/2)]
    rhs_cov_hats = some_cov_hats[(length(some_cov_hats)/2+1):length(some_cov_hats)]

    lhs_mean_cov = Reduce("+", lhs_cov_hats) / length(lhs_cov_hats)
    rhs_mean_cov = Reduce("+", rhs_cov_hats) / length(rhs_cov_hats)

    lhs_mean_covs[[j]] = lhs_mean_cov
    rhs_mean_covs[[j]] = rhs_mean_cov
}


data_diffs = data.frame()

for(j in seq_along(all_Ks)){
    this_K = all_Ks[j]

    this_emp_diff = norm(lhs_emp_covs[[j]] - rhs_emp_covs[[j]], "2")
    this_mean_diff = norm(lhs_mean_covs[[j]] - rhs_mean_covs[[j]], "2")

    this_emp_scaled_diff = this_emp_diff * this_K
    this_mean_scaled_diff = this_mean_diff * this_K

    this_emp_rel_diff = this_emp_diff / norm(lhs_emp_covs[[j]], "2")
    this_mean_rel_diff = this_mean_diff / norm(lhs_mean_covs[[j]], "2")

    this_diff = c(this_K, this_emp_diff, this_mean_diff, this_emp_scaled_diff, this_mean_scaled_diff, this_emp_rel_diff, this_mean_rel_diff)
    data_diffs = rbind(data_diffs, this_diff)
}    
colnames(data_diffs) = c("K", "Emp", "Mean", "Emp-Scaled", "Mean-Scaled", "Emp-Rel", "Mean-Rel")

data_diffs








#* Compute summary statistics of ratios of diagonal elements across pairs of K-levels

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










#* Trajectory of errors across increasing MC sizes

# list_emp_covs
# list_all_errs

# Extract empirical covariance and mean estimated covariance.
for(j in seq_along(all_Ks)){

    this_emp_cov = cov(list_par_hats[[j]])
    list_emp_covs[[j]] = this_emp_cov

    some_cov_hats = list_par_cov_hats[[j]]
    this_mean_cov = Reduce("+", some_cov_hats) / length(some_cov_hats)

    list_mean_covs[[j]] = this_mean_cov

    list_all_errs[[j]] = lapply(some_cov_hats, function(x) x - this_emp_cov)

}

list_mean_covs[[5]]


norms_by_MC_size = data.frame()

total_num_reps = nrow(list_par_hats[[1]])
# all_Bs = seq(100, 200, by=10)
all_Bs = seq(100, total_num_reps, by=50)
# all_Bs = seq(100, total_num_reps)

list_emp_covs_by_MC_size = lapply(list_emp_covs, cov)


for(i in seq_along(all_Ks)){

    some_emp_covs = list()
    for(r in seq_along(all_Bs)){

        print(paste0("K = ", all_Ks[i], " of ", max(all_Ks), ", B = ", all_Bs[r], " of ", total_num_reps))

        this_B = all_Bs[r]

        # this_emp_cov = cov(list_par_hats2[[i]][1:this_B,])
        # this_mean_cov = Reduce("+", list_par_cov_hats2[[i]][1:this_B]) / this_B

        # # Pool both simulations' results
        this_emp_cov = cov(list_par_hats[[i]][1:this_B,])
        this_mean_cov = Reduce("+", list_par_cov_hats[[i]][1:this_B]) / this_B

        some_emp_covs[[r]] = this_emp_cov

        this_err = norm(this_emp_cov - this_mean_cov, type = "2")

        this_err_scaled = this_err * all_Ks[i]
        this_err_extra_scaled = this_err_scaled * sqrt(all_Ks[i])

        this_info = c(all_Ks[i], this_B, this_err, this_err_scaled, this_err_extra_scaled)


        norms_by_MC_size = rbind(norms_by_MC_size, this_info)
    }

    list_emp_covs_by_MC_size[[i]] = some_emp_covs
}
colnames(norms_by_MC_size) = c("K", "MC_Size", "Abs_Error", "Scaled_Error", "Extra_Scaled_Error")




emp_covs_summary = lapply(list_emp_covs_by_MC_size, function(x) sapply(x, function(y) c(norm(y, type = "2"), det(y))))
data_emp_covs_summary = data.frame()
for(i in seq_along(all_Ks)){
    this_K = all_Ks[i]

    this_summary = t(emp_covs_summary[[i]])

    some_summary_data = data.frame()

    for(j in 2:length(all_Bs)){
        this_B = all_Bs[j]

        this_var_norm = var(this_summary[1:j, 1])
        this_var_det = var(this_summary[1:j, 2])

        this_info = c(this_K, this_B, this_var_norm, this_var_det)
        some_summary_data = rbind(some_summary_data, this_info)
    }
    colnames(some_summary_data) = c("K", "MC_Size", "Norm_Var", "Det_Var")

    data_emp_covs_summary = rbind(data_emp_covs_summary, some_summary_data)
}

data_emp_covs_summary %>%
    filter(K == 800) %>%
    lm(log(Norm_Var) ~ log(MC_Size), data = .) %>%
    summary()

data_emp_covs_summary %>%
    filter(K == 800) %>%
    lm(log(Det_Var) ~ log(MC_Size), data = .) %>%
    summary()

# norms_by_MC_size




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
    # geom_point() +
    facet_wrap(~K, scales="free") +
    theme_bw()
dev.off()

pdf("Plots/Scaled_Error_by_MC_Size.pdf", width = 10, height = 10)
norms_by_MC_size %>%
    ggplot(aes(x = MC_Size, y = Scaled_Error, color = as.factor(K))) +
    geom_line() +
    # geom_point() +
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
    ggplot(aes(x = K, y = Scaled_Error)) +
    geom_line() +
    geom_point() +
    theme_bw()


# Estimate rate at which errors scale with K

pdf("Plots/log_Abs_Err_vs_log_K.pdf", width = 10, height = 10)
norms_by_MC_size %>%
    filter(MC_Size == max(MC_Size)) %>%
    ggplot(aes(x = log(K), y = log(Abs_Error))) +
    geom_line() +
    geom_point() +
    theme_bw()
dev.off()


best_err_estimates = norms_by_MC_size %>%
    filter(MC_Size == max(MC_Size)) %>%
    mutate(log_K = log(K), log_err = log(Abs_Error))

fit_err_rate = lm(log_err ~ log_K, data = best_err_estimates)
# fit_err_rate = lm(log_err ~ log_K, data = best_err_estimates[-4,])
summary(fit_err_rate)

plot(fit_err_rate, 1)




filter(norms_by_MC_size, MC_Size == max(MC_Size))





#* Simple timing study

# some_Ks = rep(all_Ks, each=3)
# some_Ns = rep(c(20, 40, 60), times = length(all_Ks))

# times = c(14, 19, 26, 19, 41, 37, 29, 37, 56, 33, 57, 78, 40, 56, 86)

# data_times = data.frame(K = some_Ks, N = some_Ns, time = times)

# fit_times = lm(time ~ K + N, data = data_times)
# summary(fit_times)

# plot(data_times$K, data_times$time)
# plot(data_times$N, data_times$time)

# data_times %>% group_by(N) %>% summarize(mean_time = mean(time), sd_time = sd(time))
# data_times %>% group_by(K) %>% summarize(mean_time = mean(time), sd_time = sd(time))

# # plot time vs K for each level of N
# ggplot(data = data_times, aes(x = K, y = time)) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~N)

# ggplot(data = data_times, aes(x = N, y = time)) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~K)

# data_pred = data.frame(N=rep(200, times=length(all_Ks)), K=all_Ks)
# sum(predict(fit_times, newdata = data_pred)) / 360






#* Explore distribution of errors

load("Par_Hat_MC-Large_K_Pooled.RData", verbose = TRUE)

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



data_all_errs = data.frame()

for(i in seq_along(all_Ks)){
    some_errs = list_all_errs[[i]]

    this_data_errs = data.frame()
    for(j in seq_along(some_errs)){

        if(j %% 100 == 0) print(paste("On K = ", all_Ks[i], ", iteration ", j, " out of ", length(some_errs),  sep=""))

        this_err = some_errs[[j]]

        # Extract upper triangle of this_err with diagonal
        this_triangle = this_err[upper.tri(this_err, diag = T)]

        this_data_errs = rbind(this_data_errs, data.frame(K = all_Ks[i], t(this_triangle)))
    }

    data_all_errs = rbind(data_all_errs, this_data_errs)
}



# pdf("Plots/Visualize_All_Errors_Theta_Hat.pdf", width = 10, height = 10)
ggplot(data_all_errs, aes(X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, X6=X6, X7=X7, X8=X8, X9=X9, X10=X10, X11=X11, X12=X12, X13=X13, X14=X14, X15=X15, X16=X16, X17=X17, X18=X18, X19=X19, X20=X20, X21=X21, X22=X22, X23=X23, X24=X24, X25=X25, X26=X26, X27=X27, X28=X28, X29=X29, X30=X30, X31=X31, X32=X32, X33=X33, X34=X34, X35=X35, X36=X36, X37=X37, X38=X38, X39=X39, X40=X40, X41=X41, X42=X42, X43=X43, X44=X44, X45=X45, X46=X46, X47=X47, X48=X48, X49=X49, X50=X50, X51=X51, X52=X52, X53=X53, X54=X54, X55=X55, X56=X56, X57=X57, X58=X58, X59=X59, X60=X60, X61=X61, X62=X62, X63=X63, 
X64=X64, X65=X65, X66=X66, X67=X67, X68=X68, X69=X69, X70=X70, X71=X71, X72=X72, X73=X73, X74=X74, X75=X75, X76=X76, X77=X77, X78=X78, X79=X79, X80=X80, X81=X81, X82=X82, X83=X83, X84=X84, X85=X85, X86=X86, X87=X87, X88=X88, X89=X89, X90=X90, X91=X91, X92=X92, X93=X93, X94=X94, X95=X95, X96=X96, X97=X97, X98=X98, X99=X99, X100=X100, X101=X101, X102=X102, X103=X103, X104=X104, X105=X105, X106=X106, X107=X107, X108=X108, X109=X109, X110=X110, X111=X111, X112=X112, X113=X113, X114=X114, X115=X115, X116=X116, X117=X117, X118=X118, X119=X119, X120=X120)) + 
    geom_path(alpha=0.1) + coord_serialaxes() + facet_wrap(~K, scales = "free_y") #+ geom_histogram()
# dev.off()

data_all_errs %>%
  filter(K == 400) %>% slice(-c(bad_rows)) %>%
    ggplot(aes(X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, X6=X6, X7=X7, X8=X8, X9=X9, X10=X10, X11=X11, X12=X12, X13=X13, X14=X14, X15=X15, X16=X16, X17=X17, X18=X18, X19=X19, X20=X20, X21=X21, X22=X22, X23=X23, X24=X24, X25=X25, X26=X26, X27=X27, X28=X28, X29=X29, X30=X30, X31=X31, X32=X32, X33=X33, X34=X34, X35=X35, X36=X36, X37=X37, X38=X38, X39=X39, X40=X40, X41=X41, X42=X42, X43=X43, X44=X44, X45=X45, X46=X46, X47=X47, X48=X48, X49=X49, X50=X50, X51=X51, X52=X52, X53=X53, X54=X54, X55=X55, X56=X56, X57=X57, X58=X58, X59=X59, X60=X60, X61=X61, X62=X62, X63=X63, 
    X64=X64, X65=X65, X66=X66, X67=X67, X68=X68, X69=X69, X70=X70, X71=X71, X72=X72, X73=X73, X74=X74, X75=X75, X76=X76, X77=X77, X78=X78, X79=X79, X80=X80, X81=X81, X82=X82, X83=X83, X84=X84, X85=X85, X86=X86, X87=X87, X88=X88, X89=X89, X90=X90, X91=X91, X92=X92, X93=X93, X94=X94, X95=X95, X96=X96, X97=X97, X98=X98, X99=X99, X100=X100, X101=X101, X102=X102, X103=X103, X104=X104, X105=X105, X106=X106, X107=X107, X108=X108, X109=X109, X110=X110, X111=X111, X112=X112, X113=X113, X114=X114, X115=X115, X116=X116, X117=X117, X118=X118, X119=X119, X120=X120)) + 
    geom_path(alpha=0.1) + coord_serialaxes()



data_errs_800_full %>%
     slice(-c(bad_rows)) %>%
    ggplot(aes(X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, X6=X6, X7=X7, X8=X8, X9=X9, X10=X10, X11=X11, X12=X12, X13=X13, X14=X14, X15=X15, X16=X16, X17=X17, X18=X18, X19=X19, X20=X20, X21=X21, X22=X22, X23=X23, X24=X24, X25=X25, X26=X26, X27=X27, X28=X28, X29=X29, X30=X30, X31=X31, X32=X32, X33=X33, X34=X34, X35=X35, X36=X36, X37=X37, X38=X38, X39=X39, X40=X40, X41=X41, X42=X42, X43=X43, X44=X44, X45=X45, X46=X46, X47=X47, X48=X48, X49=X49, X50=X50, X51=X51, X52=X52, X53=X53, X54=X54, X55=X55, X56=X56, X57=X57, X58=X58, X59=X59, X60=X60, X61=X61, X62=X62, X63=X63, 
    X64=X64, X65=X65, X66=X66, X67=X67, X68=X68, X69=X69, X70=X70, X71=X71, X72=X72, X73=X73, X74=X74, X75=X75, X76=X76, X77=X77, X78=X78, X79=X79, X80=X80, X81=X81, X82=X82, X83=X83, X84=X84, X85=X85, X86=X86, X87=X87, X88=X88, X89=X89, X90=X90, X91=X91, X92=X92, X93=X93, X94=X94, X95=X95, X96=X96, X97=X97, X98=X98, X99=X99, X100=X100, X101=X101, X102=X102, X103=X103, X104=X104, X105=X105, X106=X106, X107=X107, X108=X108, X109=X109, X110=X110, X111=X111, X112=X112, X113=X113, X114=X114, X115=X115, X116=X116, X117=X117, X118=X118, X119=X119, X120=X120)) + 
    geom_path(alpha=0.1) + coord_serialaxes()



# Detect and remove outliers
data_errs_800_full = data_all_errs %>% filter(K == 800) %>% dplyr::select(-K)

SD_errs_800_full = data_errs_800_full %>% summarise_all(sd)

hist(unlist(SD_errs_800_full))

bad_cols = which(unlist(SD_errs_800_full) > 2e-04)


hist(data_errs_800_full[,bad_cols[1]])
hist(data_errs_800_full[,bad_cols[1]], breaks = 100)
bad_rows = which(data_errs_800_full[,bad_cols[1]] < -4e-3)

hist(unlist(data_errs_800_full[bad_cols[2]]), breaks = 100)
bad_rows = c(bad_rows, which(unlist(data_errs_800_full[bad_cols[2]]) < -1.2e-3)) 
# bad_rows = c(bad_rows, which(unlist(data_errs_800_full[bad_cols[2]]) < -5e-4 | unlist(data_errs_800_full[bad_cols[2]]) > 0)) # Verify that removing lots of rows does decrease the SD

# data_errs_800_full %>% slice(-c(bad_rows)) %>% dplyr::select(X55) %>% unlist() %>% hist(breaks = 100)

hist(unlist(data_errs_800_full[bad_cols[3]]), breaks=100)
bad_rows = c(bad_rows, which(unlist(data_errs_800_full[bad_cols[3]]) < -1e-3 | unlist(data_errs_800_full[bad_cols[3]]) > 5e-4)) 
# bad_rows = c(bad_rows, which(unlist(data_errs_800_full[bad_cols[3]]) < -7e-4 | unlist(data_errs_800_full[bad_cols[3]]) > 0))# Verify that removing lots of rows does decrease the SD


bad_rows %<>% unique()


data_errs_800 = data_all_errs %>% filter(K == 800) %>% dplyr::select(-K) %>% slice(-c(bad_rows))
SD_errs_800 = data_errs_800 %>% summarise_all(sd)
hist(unlist(SD_errs_800))
which(SD_errs_800 > 2e-04)



hist(unlist(data_errs_800_full[10]))#, breaks=100)




# Detect and remove outliers - K=400
data_errs_400_full = data_all_errs %>% filter(K == 400) %>% dplyr::select(-K)

SD_errs_400_full = data_errs_400_full %>% summarise_all(sd)

hist(unlist(SD_errs_400_full))

bad_cols = which(unlist(SD_errs_400_full) > 4e-04)


hist(data_errs_400_full[,bad_cols[1]])
hist(data_errs_400_full[,bad_cols[1]], breaks = 100)
bad_rows = which(data_errs_400_full[,bad_cols[1]] < -1e-2)

hist(unlist(data_errs_400_full[bad_cols[2]]), breaks = 100)
bad_rows = c(bad_rows, which(unlist(data_errs_400_full[bad_cols[2]]) < -4e-3)) 


hist(unlist(data_errs_400_full[bad_cols[3]]), breaks=100)
bad_rows = c(bad_rows, which(unlist(data_errs_400_full[bad_cols[3]]) < -4e-3)) 

hist(unlist(data_errs_400_full[bad_cols[4]]), breaks=100, xlim = c(-3e-3, 3e-3))
bad_rows = c(bad_rows, which(unlist(data_errs_400_full[bad_cols[3]]) < -2.5e-3 | unlist(data_errs_400_full[bad_cols[3]]) > 2.3e-3)) 



bad_rows %<>% unique()


data_errs_400 = data_all_errs %>% filter(K == 400) %>% dplyr::select(-K) %>% slice(-c(bad_rows))
SD_errs_400 = data_errs_400 %>% summarise_all(sd)
hist(unlist(SD_errs_400))
which(SD_errs_400 > 2e-04)




#* Explore outlier detection for Theta hat

some_par_hats = list_par_hats[[5]]

all_grubbs_p_vals = data.frame()
for(i in seq_along(list_par_hats)){
    some_par_hats = list_par_hats[[i]]
    
    some_p_vals = sapply(seq_len(ncol(some_par_hats)), function(i) grubbs.test(some_par_hats[,i])$p.value)

    this_info = c(all_Ks[i], some_p_vals)

    all_grubbs_p_vals = rbind(all_grubbs_p_vals, this_info)
}
colnames(all_grubbs_p_vals) = c("K", "b_Y_Int", "b_Y_X", "b_Y_M", "b_Y_W1", "b_Y_W2", "theta_Y_Int", "rho_Y", "theta_Y_X", "b_M_Int", "b_M_X", "b_M_W1", "b_M_W2", "theta_M_Int", "rho_M", "theta_M_X")
all_grubbs_p_vals

conclusions_grubbs = all_grubbs_p_vals < 0.05

small_p_vals = all_grubbs_p_vals[conclusions_grubbs]

which(conclusions_grubbs, arr.ind=T)

all_grubbs_p_vals[1,10]



## Plot theta hats for each K
data_theta_hats = data.frame()
for(i in seq_along(list_par_hats)){
    some_theta_hats = list_par_hats[[i]]

    this_data = data.frame(K = rep(all_Ks[i], nrow(some_theta_hats)), some_theta_hats)
    data_theta_hats = rbind(data_theta_hats, this_data)
}
colnames(data_theta_hats) = c("K", "b_Y_Int", "b_Y_X", "b_Y_M", "b_Y_W1", "b_Y_W2", "theta_Y_Int", "rho_Y", "theta_Y_X", "b_M_Int", "b_M_X", "b_M_W1", "b_M_W2", "theta_M_Int", "rho_M", "theta_M_X")

ggplot(data_theta_hats, aes(b_Y_Int=b_Y_Int, b_Y_X=b_Y_X, b_Y_M=b_Y_M, b_Y_W1=b_Y_W1, b_Y_W2=b_Y_W2, theta_Y_Int=theta_Y_Int, rho_Y=rho_Y, theta_Y_X=theta_Y_X, b_M_Int=b_M_Int, b_M_X=b_M_X, b_M_W1=b_M_W1, b_M_W2=b_M_W2, theta_M_Int=theta_M_Int, rho_M=rho_M, theta_M_X=theta_M_X)) + geom_path(alpha=0.1) + coord_serialaxes() + facet_wrap(~K, scales = "free_y")



## Plot Sigma_0 hat for each K
data_Sigma_hat_raw = data.frame()

for(i in seq_along(all_Ks)){
    some_cov_hats = list_par_cov_hats[[i]]

    this_data_cov = data.frame()
    for(j in seq_along(some_cov_hats)){

        if(j %% 100 == 0) print(paste("On K = ", all_Ks[i], ", iteration ", j, " out of ", length(some_cov_hats),  sep=""))

        this_cov_hat = some_cov_hats[[j]]

        # Extract upper triangle of this_err with diagonal
        this_triangle = this_cov_hat[upper.tri(this_cov_hat, diag = T)]

        this_data_cov = rbind(this_data_cov, data.frame(K = all_Ks[i], t(this_triangle)))
    }

    this_data_cov$empirical = F

    this_emp_cov = cov(list_par_hats[[i]])
    this_triangle = this_emp_cov[upper.tri(this_emp_cov, diag = T)]

    this_data_cov = rbind(this_data_cov, data.frame(K = all_Ks[i], t(this_triangle), empirical = T))

    data_Sigma_hat_raw = rbind(data_Sigma_hat_raw, this_data_cov)
}

### Remove identically zero columns
find_non_zero_cols = function(mat) sapply(seq_len(ncol(mat)), function(i) !all(mat[,i] == 0))

to_keep = data_Sigma_hat_raw %>%
    filter(empirical == F) %>%
    foo()
to_keep[length(to_keep)] = T

data_Sigma_hat = data_Sigma_hat_raw[,to_keep]
# data_Sigma_hat = data_Sigma_hat_raw

# Get list of all relevant variable names (some elements of Sigma_0 are identically zero)
paste(paste(colnames(data_Sigma_hat[,-1]), colnames(data_Sigma_hat[,-1]), sep="="), collapse = ", ")

data_Sigma_hat %<>% mutate(opacity = ifelse(empirical == T, 1, 0.1))

data_baseline = data.frame()
for(i in seq_along(all_Ks)){
    for(j in -1:1){
        this_baseline = rep(j/all_Ks[i], times = ncol(data_Sigma_hat)-3)
        this_baseline = c(all_Ks[i], this_baseline, T, 1)

        data_baseline = rbind(data_baseline, this_baseline)
    }
}
colnames(data_baseline) = colnames(data_Sigma_hat)
data_Sigma_hat = rbind(data_baseline, data_Sigma_hat)

ggplot(data_Sigma_hat, aes(X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, X6=X6, X7=X7, X8=X8, X9=X9, X10=X10, X11=X11, X12=X12, X13=X13, X14=X14, X15=X15, X16=X16, X17=X17, X18=X18, X19=X19, X20=X20, X21=X21, X22=X22, X23=X23, X24=X24, X25=X25, X26=X26, X27=X27, X28=X28, X29=X29, X30=X30, X31=X31, X32=X32, X33=X33, X34=X34, X35=X35, X36=X36, X45=X45, X54=X54, X55=X55, X64=X64, X65=X65, X66=X66, X75=X75, X76=X76, X77=X77, X78=X78, X87=X87, X88=X88, X89=X89, X90=X90, X91=X91, X100=X100, X101=X101, X102=X102, X103=X103, X104=X104, X105=X105, X114=X114, X115=X115, X116=X116, X117=X117, X118=X118, X119=X119, X120=X120, color = empirical)) + geom_path(aes(alpha=opacity)) + coord_serialaxes() + facet_wrap(~K, scales = "free_y")


## Detect outliers in estimated covariance matrices
all_grubbs_p_vals = data.frame()
for(i in seq_along(list_par_hats)){

    some_cov_hats = filter(data_Sigma_hat, K == all_Ks[i]) %>% select(-K, -empirical, -opacity)
    
    some_p_vals = sapply(seq_len(ncol(some_cov_hats)), function(i) grubbs.test(some_cov_hats[,i])$p.value)

    this_info = c(all_Ks[i], some_p_vals)

    all_grubbs_p_vals = rbind(all_grubbs_p_vals, this_info)
}

all_grubbs_p_vals * ncol(all_grubbs_p_vals)

all_grubbs_p_vals %>% select(-1) %>% apply(.,1, mean)


## Remove rows with any value less than -1/K

rows_keep = data_Sigma_hat %>% select(-empirical, -opacity) %>% apply(., 1, function(x) !any(x < -1/x[1]) & !any(x[-1] > 10/x[1]))

data_Sigma_hat_clean = data_Sigma_hat[rows_keep,]

ggplot(data_Sigma_hat_clean, aes(X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, X6=X6, X7=X7, X8=X8, X9=X9, X10=X10, X11=X11, X12=X12, X13=X13, X14=X14, X15=X15, X16=X16, X17=X17, X18=X18, X19=X19, X20=X20, X21=X21, X22=X22, X23=X23, X24=X24, X25=X25, X26=X26, X27=X27, X28=X28, X29=X29, X30=X30, X31=X31, X32=X32, X33=X33, X34=X34, X35=X35, X36=X36, X45=X45, X54=X54, X55=X55, X64=X64, X65=X65, X66=X66, X75=X75, X76=X76, X77=X77, X78=X78, X87=X87, X88=X88, X89=X89, X90=X90, X91=X91, X100=X100, X101=X101, X102=X102, X103=X103, X104=X104, X105=X105, X114=X114, X115=X115, X116=X116, X117=X117, X118=X118, X119=X119, X120=X120, color = empirical)) + geom_path(aes(alpha=opacity)) + coord_serialaxes() + facet_wrap(~K, scales = "free_y")



# ## Remove outliers from original datasets
# inds_remove = which(!rows_keep)
# K_inds_remove = (inds_remove %/% total_num_reps) + 1
# row_inds_remove = inds_remove %% total_num_reps

# list_par_hats_clean = list_par_hats
# list_par_cov_hats_clean = list_par_cov_hats

# for(ind in seq_along(inds_remove)){
#     i = K_inds_remove[ind]
#     j = row_inds_remove[ind]
#     list_par_hats_clean[[i]] = list_par_hats_clean[[i]][-j,]
#     list_par_cov_hats_clean[[i]] = list_par_cov_hats_clean[[i]][-j]
# }


# save(list_par_hats_clean, list_par_cov_hats_clean, num_reps, all_Ks, file = "Par_Hat_MC_Clean.RData")
load(file = "Par_Hat_MC_Clean.RData", verbose = TRUE)


list_emp_covs_clean = list()
list_mean_covs_clean = list()
list_all_errs_clean = list()

# Extract empirical covariance and mean estimated covariance.
for(j in seq_along(all_Ks)){

    this_emp_cov = cov(list_par_hats_clean[[j]])
    list_emp_covs_clean[[j]] = this_emp_cov

    some_cov_hats = list_par_cov_hats_clean[[j]]
    this_mean_cov = Reduce("+", some_cov_hats) / length(some_cov_hats)

    list_mean_covs_clean[[j]] = this_mean_cov

    list_all_errs_clean[[j]] = lapply(some_cov_hats, function(x) x - this_emp_cov)

}


all_mean_errs = sapply(list_all_errs_clean, function(some_errs) mean(sapply(some_errs, norm, type="2")))
# all_mean_errs = sapply(list_all_errs, function(some_errs) mean(sapply(some_errs, norm, type="2")))

plot(log(all_Ks), log(all_mean_errs))

data_mean_errs = data.frame(K = all_Ks, err = all_mean_errs)

data_mean_errs %>% lm(log(err) ~ log(K), data=.) %>% summary()
data_mean_errs %>% slice(-4) %>% lm(log(err) ~ log(K), data=.) %>% summary()




#* Explore UQ for ENCs with clean data


# list_ENC_hats_clean = list()
# list_ENC_cov_hats_clean = list()

# for(i in seq_along(all_Ks)){

#     some_ENC_hats = data.frame()
#     some_ENC_cov_hats = list()

#     this_num_reps = nrow(list_par_hats_clean[[i]])

#     for(j in seq_len(this_num_reps)){
#         if(j %% 50 == 0) {
#             print(paste0("j = ", j, " of ", this_num_reps, ", K = ", all_Ks[i], " (number ", i, " of ", length(all_Ks), ")"))
#         }

#         this_par_hat = list_par_hats_clean[[i]][j,]

#         this_b_Y = this_par_hat[1:5]
#         this_theta_Y = this_par_hat[6:8]
#         this_b_M = this_par_hat[9:12]
#         this_theta_M = this_par_hat[13:15]


#         this_ENC_hat = all_ENCs(w, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
#         some_ENC_hats = rbind(some_ENC_hats, this_ENC_hat)


#         this_Sigma = list_par_cov_hats_clean[[i]][[j]]
#         this_ENC_cov_hat = all_covs_ENC_pars(w, this_Sigma, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
#         some_ENC_cov_hats[[j]] = this_ENC_cov_hat

#     }

#     colnames(some_ENC_hats) = c("11", "10", "01", "00")
#     list_ENC_hats_clean[[i]] = some_ENC_hats
#     list_ENC_cov_hats_clean[[i]] = some_ENC_cov_hats
# }

# save(list_ENC_hats_clean, list_ENC_cov_hats_clean, num_reps, all_Ks, file = "ENC_Hat_MC_Clean.RData")
load(file = "ENC_Hat_MC_Clean.RData", verbose = TRUE)


list_ENC_emp_covs = lapply(list_ENC_hats_clean, cov)
list_ENC_mean_covs = lapply(list_ENC_cov_hats_clean, function(x) Reduce("+", x) / length(x))
list_ENC_errs = lapply(seq_along(list_ENC_cov_hats_clean), function(i) lapply(list_ENC_cov_hats_clean[[i]], function(y) y - list_ENC_emp_covs[[i]]))

ENC_mean_errs = sapply(list_ENC_errs, function(some_errs) mean(sapply(some_errs, norm, type="2")))

plot(log(all_Ks), log(ENC_mean_errs))

data_ENC_mean_errs = data.frame(K = all_Ks, err = ENC_mean_errs)    

data_ENC_mean_errs %>% lm(log(err) ~ log(K), data=.) %>% summary()
# data_ENC_mean_errs %>% slice(-5) %>% lm(log(err) ~ log(K), data=.) %>% summary()



#* Mediation Effects



list_ME_hats_clean = list()
list_ME_cov_hats_clean = list()

this_scale = c("diff", "rat", "OR")
ME_names = expand.grid(flavour = c("total", "direct", "indirect"), scale = this_scale) %>% 
        arrange(flavour) %>% apply(1, paste, collapse = "_")

for(i in seq_along(all_Ks)){

    some_ME_hats = data.frame()
    some_ME_cov_hats = list()

    some_par_hats = list_par_hats_clean[[i]]
    this_num_reps = nrow(some_par_hats)

    for(j in seq_len(this_num_reps)){
        if(j %% 100 == 0) {
            print(paste0("j = ", j, " of ", this_num_reps, ", K = ", all_Ks[i], " (number ", i, " of ", length(all_Ks), ")"))
        }

        this_par_hat = some_par_hats[j,]

        this_b_Y = this_par_hat[1:5]
        this_theta_Y = this_par_hat[6:8]
        this_b_M = this_par_hat[9:12]
        this_theta_M = this_par_hat[13:15]


        this_ME_hat = all_MEs_pars(this_scale, w, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
        some_ME_hats = rbind(some_ME_hats, this_ME_hat)

        q = unlist(all_MEs_ENCs(this_scale, list_ENC_hats_clean[[i]][j,], which_REs=which_REs))


        this_Sigma = list_par_cov_hats_clean[[i]][[j]]
        this_ME_cov_hat = all_covs_MEs_pars(this_scale, w, this_Sigma, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
        some_ME_cov_hats[[j]] = this_ME_cov_hat

    }

    
    colnames(some_ME_hats) = ME_names
    list_ME_hats_clean[[i]] = some_ME_hats
    list_ME_cov_hats_clean[[i]] = some_ME_cov_hats
}



# save(list_ME_hats_clean, list_ME_cov_hats_clean, num_reps, all_Ks, file = "ME_Hat_MC_Clean.RData")
load(file = "ME_Hat_MC_Clean.RData", verbose = TRUE)


list_ME_emp_covs = lapply(list_ME_hats_clean, cov)
list_ME_mean_covs = lapply(list_ME_cov_hats_clean, function(x) Reduce("+", x) / length(x))
list_ME_errs = lapply(seq_along(list_ME_cov_hats_clean), function(i) lapply(list_ME_cov_hats_clean[[i]], function(y) y - list_ME_emp_covs[[i]]))

ME_mean_errs = sapply(list_ME_errs, function(some_errs) mean(sapply(some_errs, norm, type="2")))

plot(log(all_Ks), log(ME_mean_errs))

data_ME_mean_errs = data.frame(K = all_Ks, err = ME_mean_errs)    

data_ME_mean_errs %>% lm(log(err) ~ log(K), data=.) %>% summary()
data_ME_mean_errs %>% lm(log(err) ~ log(K), data=.) %>% abline()

