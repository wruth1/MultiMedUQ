



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



# load("Par_Hat_MC-Large_K_Pooled.RData", verbose = TRUE)




# # Extract empirical and average estimated covariances for each value of K

# num_pars = ncol(list_par_hats[[1]])


# ## Remove any runs with NA parameter estimates
# for(j in seq_along(all_Ks)){
#     # Covariance matrices first to preserve information in parameter estimates
#     list_par_cov_hats[[j]] = list_par_cov_hats[[j]][complete.cases(list_par_hats[[j]])]
#     list_par_hats[[j]] = na.omit(list_par_hats[[j]])
# }

# lengths(list_par_cov_hats)
# sapply(list_par_hats, nrow)


# list_emp_covs = list()
# list_mean_covs = list()
# list_all_errs = list()

# # Extract empirical covariance and mean estimated covariance.
# for(j in seq_along(all_Ks)){

#     this_emp_cov = cov(list_par_hats[[j]])
#     list_emp_covs[[j]] = this_emp_cov

#     some_cov_hats = list_par_cov_hats[[j]]
#     this_mean_cov = Reduce("+", some_cov_hats) / length(some_cov_hats)

#     list_mean_covs[[j]] = this_mean_cov

#     list_all_errs[[j]] = lapply(some_cov_hats, function(x) x - this_emp_cov)

# }



# #* Covariances of ENCs

# list_ENC_hats = list()
# list_ENC_cov_hats = list()

# num_reps = nrow(list_par_hats[[1]])

# for(i in seq_along(all_Ks)){

#     some_ENC_hats = data.frame()
#     some_ENC_cov_hats = list()

#     for(j in seq_len(num_reps)){
#         if(j %% 50 == 0) {
#             print(paste0("j = ", j, " of ", num_reps, ", K = ", all_Ks[i], " (number ", i, " of ", length(all_Ks), ")"))
#         }

#         this_par_hat = list_par_hats[[i]][j,]

#         this_b_Y = this_par_hat[1:5]
#         this_theta_Y = this_par_hat[6:8]
#         this_b_M = this_par_hat[9:12]
#         this_theta_M = this_par_hat[13:15]


#         this_ENC_hat = all_ENCs(w, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
#         some_ENC_hats = rbind(some_ENC_hats, this_ENC_hat)


#         this_Sigma = list_par_cov_hats[[i]][[j]]
#         this_ENC_cov_hat = all_covs_ENC_pars(w, this_Sigma, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
#         some_ENC_cov_hats[[j]] = this_ENC_cov_hat

#     }

#     colnames(some_ENC_hats) = c("11", "10", "01", "00")
#     list_ENC_hats[[i]] = some_ENC_hats
#     list_ENC_cov_hats[[i]] = some_ENC_cov_hats
# }

# save(all_Ks, num_reps, list_ENC_hats, list_ENC_cov_hats, file = "ENC_Cov_Hats.RData")
load("ENC_Cov_Hats.RData", verbose = TRUE)



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

info_ENC_norms = data.frame()
for(i in seq_along(all_Ks)){
    diff_mat = list_ENC_emp_covs[[i]] - list_ENC_mean_covs[[i]]
    # diff_mat_norm = norm(diff_mat, type = "2")
    diff_mat_norm = norm(diff_mat, type = "F")

    scaled_norm = diff_mat_norm * all_Ks[i]

    this_info = c(all_Ks[i], diff_mat_norm, scaled_norm)
    info_ENC_norms = rbind(info_ENC_norms, this_info)
}
colnames(info_ENC_norms) = c("K", "Abs_Err", "Scaled_Err")

info_ENC_norms

scaled_ENC_norms = info_ENC_norms






# Trajectory of errors across increasing MC sizes


# i = 5
# (Sigma_emp = list_emp_covs[[i]] * all_Ks[i])
# (Sigma_mean = list_mean_covs[[i]] * all_Ks[i])
# (Delta = Sigma_emp - Sigma_mean)

# eigen(Sigma_emp, symmetric=T, only.values = T)$values
# eigen(Sigma_mean, symmetric=T, only.values = T)$values



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
        this_emp_cov = cov(list_ENC_hats[[i]][1:this_B,])
        this_mean_cov = Reduce("+", list_ENC_cov_hats[[i]][1:this_B]) / this_B

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
# pdf("Plots/Abs_Error_by_MC_Size.pdf", width = 10, height = 10)
norms_by_MC_size %>%
    ggplot(aes(x = MC_Size, y = Abs_Error, color = as.factor(K))) +
    geom_line() +
    geom_point() +
    facet_wrap(~K, scales="free") +
    theme_bw()
# dev.off()

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

# pdf("Plots/log_Abs_Err_vs_log_K.pdf", width = 10, height = 10)
norms_by_MC_size %>%
    filter(MC_Size == max(MC_Size)) %>%
    ggplot(aes(x = log(K), y = log(Abs_Error))) +
    geom_line() +
    geom_point() +
    theme_bw()
# dev.off()


best_err_estimates = norms_by_MC_size %>%
    filter(MC_Size == max(MC_Size)) %>%
    mutate(log_K = log(K), log_err = log(Abs_Error))

fit_err_rate = lm(log_err ~ log_K, data = best_err_estimates[-1,])
summary(fit_err_rate)

plot(fit_err_rate, 1)




filter(norms_by_MC_size, MC_Size == max(MC_Size))




# Investigate scaling of Jacobian matrix


list_ENC_Jacobians = list()

num_reps = nrow(list_par_hats[[1]])

for(i in seq_along(all_Ks)){

    some_ENC_Jacobians = list()

    for(j in seq_len(num_reps)){
        if(j %% 50 == 0) {
            print(paste0("j = ", j, " of ", num_reps, ", K = ", all_Ks[i], " (number ", i, " of ", length(all_Ks), ")"))
        }

        this_par_hat = list_par_hats[[i]][j,]

        this_b_Y = this_par_hat[1:5]
        this_theta_Y = this_par_hat[6:8]
        this_b_M = this_par_hat[9:12]
        this_theta_M = this_par_hat[13:15]


        this_ENC_Jacobian = Jacob_ENC_pars(w, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
        some_ENC_Jacobians[[j]] = this_ENC_Jacobian
    }

    list_ENC_Jacobians[[i]] = some_ENC_Jacobians
}


list_ENC_Jacob_norms = list()

for(i in seq_along(all_Ks)){
    some_ENC_Jacob_norms = sapply(list_ENC_Jacobians[[i]], function(x) norm(x, "2"))

    list_ENC_Jacob_norms[[i]] = some_ENC_Jacob_norms
}

mean_ENC_Jacob_norms = sapply(list_ENC_Jacob_norms, mean)
SD_ENC_Jacob_norms = sapply(list_ENC_Jacob_norms, sd)

(data_ENC_Jacob_norms = data.frame(all_Ks, mean_ENC_Jacob_norms, SD_ENC_Jacob_norms))

fit_Jacob_norm_SDs = lm(log(SD_ENC_Jacob_norms) ~ log(all_Ks), data = data_ENC_Jacob_norms)
summary(fit_Jacob_norm_SDs)
