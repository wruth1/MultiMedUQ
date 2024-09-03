



library(lme4)
library(merDeriv)
library(tictoc)
source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
devtools::load_all()


cond_num_thresh = 300   # Threshold for removal of sigma hat due to large condition number



# Set parameters and fit models

# N = 100
# N = 20
N = 200

# all_Ks = c(50, 100, 200, 400, 800)
# all_Ks = c(50, 100, 200)
all_Ks = c(50, 100)
# K = 50

# num_reps = 50
num_reps = 10
# num_reps = 5


which_REs = c("Y.Int", "Y.X", "M.All")


x = 0
x_m = 1

b_Y = c(0,0.5,1,0,0) * 2
# theta_Y = c(sqrt(0.5), 0.5, 0, 1, 0.5, sqrt(0.5))
theta_Y = c(sqrt(0.5), 0, 1)

b_M = c(0,1,0,0) * 2
theta_M = c(1, 0.5, 2)



list_par_hats = list()
list_par_cov_hats = list()

list_bad_data = list()

all_runtimes = c()

set.seed(1)

# tic()

for(j in seq_along(all_Ks)){
    tic()

    K = all_Ks[j]

    all_par_hats = c()
    all_par_cov_hats = list()

    some_bad_data = list()


    # all_med_hats = c()
    # all_cov_hats = list()



    for(i in 1:num_reps){

        # print(paste0("Rep ", i, " of ", num_reps))
        print(paste0("Rep ", i, " of ", num_reps, ", K number ", j, " of ", length(all_Ks)))

        data = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)

        w = c(0,0)


        # Singularities can arise in fitted covariance matrices. Use tryCatch to skip these cases.
        tryCatch({
            

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

        if(kappa(par_cov_hat) > cond_num_thresh) some_bad_data[[i]] = data


        }, error = function(e){
            all_par_cov_hats[[i]] = NULL
        })



        ## Mediation effects
        # med_hats = all_MEs_models(scale = "OR", w, fit_Y, fit_M, which_REs = c("Y.Int", "Y.X", "M.All")) 

        # cov_hat = all_cov_MEs(scale = "OR", w, fit_Y, fit_M, which_REs = c("Y.Int", "Y.X", "M.All"))



        # all_med_hats = rbind(all_med_hats, med_hats)
        # all_cov_hats[[i]] = cov_hat

    }

    list_par_hats[[j]] = all_par_hats
    list_par_cov_hats[[j]] = all_par_cov_hats
    list_bad_data[[j]] = some_bad_data

    this_runtime = toc()
    all_runtimes = c(all_runtimes, this_runtime$toc - this_runtime$tic)
}

all_runtimes / all_Ks

# runtimes_5 = all_runtimes
# runtimes_10 = all_runtimes
# runtimes_10 / runtimes_5

runtimes_5[2] / runtimes_5[1]
runtimes_10[2] / runtimes_10[1]
# toc()

# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC.RData")
# save(num_reps,all_Ks, list_par_hats, list_par_cov_hats, file = "Par_Hat_MC-Large_K.RData")
load("Par_Hat_MC.RData", verbose = TRUE)
load("Par_Hat_MC-Large_K.RData", verbose = TRUE)
# Extract empirical and average estimated covariances for each value of K

num_pars = ncol(list_par_hats[[1]])

list_emp_covs = list()
list_mean_covs = list()

list_par_cov_hats_clean = list()

#* Occasionally, the estimated covariance has an astronomical condition number. I suspect that this is due to poor convergence of the glmer function, but I haven't checked this.
#* I'm going to drop any covariance estimates that are too ill-conditioned (i.e. kappa > 300).


for(j in seq_along(all_Ks)){

    num_pars = ncol(list_par_hats[[1]])

    this_par_cov_hats_clean = list()

    list_emp_covs[[j]] = cov(list_par_hats[[j]])

    this_sum_covs = matrix(0, num_pars, num_pars)
    this_num_covs = 0

    for(i in 1:num_reps){
        print(i)
        this_cov_hat = list_par_cov_hats[[j]][[i]]

        if(kappa(this_cov_hat, exact = T) <= cond_num_thresh){
            this_sum_covs = this_sum_covs + list_par_cov_hats[[j]][[i]]
            this_num_covs = this_num_covs + 1

            this_par_cov_hats_clean[[this_num_covs]] = this_cov_hat
        }
    }
    list_mean_covs[[j]] = this_sum_covs / num_reps


    list_par_cov_hats_clean[[j]] = this_par_cov_hats_clean
}

sapply(list_par_cov_hats_clean, length)


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


#* Explore problematic covariances
##* Problem is with estimated covs for K = 100

q50 = list_mean_covs[[1]]
q100 = list_mean_covs[[2]]
q200 = list_mean_covs[[3]]
q400 = list_mean_covs[[4]]

eigen(q50, only.values = T)$values
eigen(q100, only.values = T)$values
eigen(q200, only.values = T)$values
eigen(q400, only.values = T)$values

kappa(q50)
kappa(q100)
kappa(q200)
kappa(q400)



bad_covs = list_par_cov_hats[[2]]
bad_cond_nums = sapply(bad_covs, kappa, exact=T)

sapply(list_par_cov_hats, function(some_covs) sum(sapply(some_covs, is.null)))

list_cond_nums = lapply(list_par_cov_hats, function(some_covs){
    sapply(some_covs, function(this_cov){
        if(!is.null(this_cov)){
            return(kappa(this_cov, exact = T))
        } else{
            return(-1)
        }
    })
} )
sapply(seq_along(list_cond_nums), function(i){
    some_cond_nums = list_cond_nums[[i]]
    return(c(all_Ks[i],mean(some_cond_nums), sd(some_cond_nums), max(some_cond_nums)))
})



### Investigate dataset leading to largest condition number
j_bad = which.max(list_cond_nums[[2]])
data_bad = list_bad_data[[2]][[j_bad]]

(fit_Y = lme4::glmer(Y ~ X + M + C1 + C2 + (X | group), data = data_bad, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6), tolPwrss = 1e-8)))
(fit_M = lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data_bad, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6), tolPwrss = 1e-8)))

bad_par_cov_hat = all_pars_cov_mat(fit_Y, fit_M)
kappa(bad_par_cov_hat)



par(mfrow = c(2,2))
for(i in seq_along(list_cond_nums)){
    hist(list_cond_nums[[i]], main = paste0("K = ", all_Ks[i]), xlab = "Condition number")
}
par(mfrow = c(1,1))


## Plot condition number vs difference between estimated and empirical covariance
### Note: There is only one empirical covariance for each K. We're exploring deviations of estimates around that fixed reference
### Note: Matrix difference can be quantified in several ways. One is mean absolute difference along the diagonal. Another is norm of the difference.

par(mfrow = c(2,2))
for(i in seq_along(all_Ks)){
        
    this_emp_cov = list_emp_covs[[i]]
    some_cov_hats = list_par_cov_hats[[i]]
    some_diff_mats = lapply(some_cov_hats, function(X) X - this_emp_cov)

    some_diag_diffs = sapply(some_diff_mats, function(X) mean(abs(diag(X))))
    some_norm_diffs = sapply(some_diff_mats, function(X) norm(X, type = "2"))

    cor(some_diag_diffs, some_norm_diffs)
    cor(list_cond_nums[[i]], some_norm_diffs)

    plot(list_cond_nums[[i]], some_norm_diffs, xlab = "Condition number", ylab = "2-Norm Difference", main = paste0("K = ", all_Ks[i], ", Corr = ", round(cor(list_cond_nums[[i]], some_norm_diffs), 2)))
}
par(mfrow = c(1,1))


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




info_rel_norms = data.frame()
for(i in seq_along(all_Ks)){

    # Compare empirical and estimated covariances
    this_emp_mat = list_emp_covs[[i]]
    this_mean_mat = list_mean_covs[[i]]

    this_norm_diff = norm(this_emp_mat - this_mean_mat, type = "2")
    this_norm_rel = this_norm_diff / norm(this_emp_mat, type = "2")

    this_info = c(all_Ks[i], this_norm_diff, this_norm_rel)

    info_rel_norms = rbind(info_rel_norms, this_info)
}
colnames(info_rel_norms) = c("K", "Norm-Diff", "Norm-Rel")

info_rel_norms


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
colnames(info_var_diffs) = c("K", "Mean-Diff", "SD-Diff", "Mean-Rel-Diff", "SD-Rel-Diff")

info_var_diffs


### Harder: Matrix norms
info_ENC_norms = data.frame()
for(i in seq_along(all_Ks)){
    diff_mat = list_ENC_emp_covs[[i]] - list_ENC_mean_covs[[i]]
    diff_mat_norm = norm(diff_mat, type = "2")

    rel_norm = diff_mat_norm / norm(list_ENC_emp_covs[[i]], type = "2")

    this_info = c(all_Ks[i], diff_mat_norm, rel_norm)
    info_ENC_norms = rbind(info_ENC_norms, this_info)
}
colnames(info_ENC_norms) = c("K", "Norm-Diff", "Rel-Norm-Diff")

info_ENC_norms





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



