

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


load(file = "Par_Hat_MC_Clean.RData", verbose = TRUE)
load(file = "ENC_Hat_MC_Clean.RData", verbose = TRUE)
load(file = "ME_Hat_MC_Clean.RData", verbose = TRUE)

which_REs = c("Y.Int", "Y.X", "M.All")



#* Theta hats
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

plot(log(all_Ks), log(all_mean_errs))

data_mean_errs = data.frame(K = all_Ks, err = all_mean_errs) %>% mutate(rel_err = err * K)

data_mean_errs %>% lm(log(err) ~ log(K), data=.) %>% summary()
data_mean_errs %>% slice(-4) %>% lm(log(err) ~ log(K), data=.) %>% summary()



#* ENC hats

list_ENC_emp_covs = lapply(list_ENC_hats_clean, cov)
list_ENC_mean_covs = lapply(list_ENC_cov_hats_clean, function(x) Reduce("+", x) / length(x))
list_ENC_errs = lapply(seq_along(list_ENC_cov_hats_clean), function(i) lapply(list_ENC_cov_hats_clean[[i]], function(y) y - list_ENC_emp_covs[[i]]))

ENC_mean_errs = sapply(list_ENC_errs, function(some_errs) mean(sapply(some_errs, norm, type="2")))

plot(log(all_Ks), log(ENC_mean_errs))

data_ENC_mean_errs = data.frame(K = all_Ks, err = ENC_mean_errs)

# data_ENC_mean_errs %>% lm(log(err) ~ log(K), data=.) %>% summary()
data_ENC_mean_errs %>% slice(-5) %>% lm(log(err) ~ log(K), data=.) %>% summary()




#* ME hats

list_ME_emp_covs = lapply(list_ME_hats_clean, cov)
list_ME_mean_covs = lapply(list_ME_cov_hats_clean, function(x) Reduce("+", x) / length(x))
list_ME_errs = lapply(seq_along(list_ME_cov_hats_clean), function(i) lapply(list_ME_cov_hats_clean[[i]], function(y) y - list_ME_emp_covs[[i]]))

ME_mean_errs = sapply(list_ME_errs, function(some_errs) mean(sapply(some_errs, norm, type="2")))

plot(log(all_Ks), log(ME_mean_errs))

data_ME_mean_errs = data.frame(K = all_Ks, err = ME_mean_errs)

data_ME_mean_errs %>% lm(log(err) ~ log(K), data=.) %>% summary()
# data_ME_mean_errs %>% slice(-5) %>% lm(log(err) ~ log(K), data=.) %>% summary()






#* Simplified MC study for MEs
set.seed(1)

w = c(1,2)

list_err_norms = list()
scale = c("diff", "rat", "OR")
# scale = "rat"

B = 1000

for(i in seq_along(all_Ks)){
        
    Gamma_1 = list_ENC_emp_covs[[i]]
    ENC_0 = apply(list_ENC_hats_clean[[i]], 2, mean)

    # Point estimates
    some_ENCs = MASS::mvrnorm(n = B, mu = ENC_0, Sigma = Gamma_1)
    some_MEs = t(apply(some_ENCs, 1, all_MEs_ENCs, scale = scale, which_REs = which_REs))

    # Estimated covariances
    some_Jacobs = lapply(seq_len(B), function(i) all_grad_MEs_ENCs(some_ENCs[i,], scale = scale, which_REs = which_REs))
    some_gamma_2_hats = lapply(some_Jacobs, function(x) x %*% Gamma_1 %*% t(x))

    # Empirical covariance
    emp_gamma_2 = cov(some_MEs)

    # Investigate empirical vs estimated Gamma 2
    some_err_norms = sapply(some_gamma_2_hats, function(x) norm(x - emp_gamma_2, type="2"))

    list_err_norms[[i]] = some_err_norms
}

mean_err_norms = sapply(list_err_norms, mean)
sd_err_norms = sapply(list_err_norms, sd)

(data_err_norms = data.frame(K = all_Ks, err = mean_err_norms, se = sd_err_norms / sqrt(B)))
data_err_norms %>% mutate(scaled_err = err * K)

plot(log(all_Ks), log(mean_err_norms))

data_err_norms %>% lm(log(err) ~ log(K), data=.) %>% summary()
data_err_norms %>% lm(log(err) ~ log(K), data=.) %>% abline()






#* Simplified MC study for ENCs (for reference)
set.seed(1)

w = c(1,2)

list_err_norms = list()
list_ENCs = list()
list_ENC_covs = list()
list_Thetas = list()

B = 1000

for(i in seq_along(all_Ks)){
    print(paste0("K = ", all_Ks[i], "; number ", i, " of ", length(all_Ks)))
        
    Gamma_0 = list_emp_covs_clean[[i]]
    Theta_0 = apply(list_par_hats_clean[[i]], 2, mean)

    # Point estimates
    some_Thetas = MASS::mvrnorm(n = B, mu = Theta_0, Sigma = Gamma_0)
    some_ENCs = t(apply(some_Thetas, 1, all_ENCs_theta, w = w, which_REs = which_REs))

    list_Thetas[[i]] = some_Thetas
    list_ENCs[[i]] = some_ENCs

    # Estimated covariances
    some_Jacobs = lapply(seq_len(B), function(i) Jacob_ENC_Theta(w, some_Thetas[i,], which_REs = which_REs))
    some_gamma_1_hats = lapply(some_Jacobs, function(x) x %*% Gamma_0 %*% t(x))
    list_ENC_covs[[i]] = some_gamma_1_hats

    # Empirical covariance
    emp_gamma_1 = cov(some_ENCs)

    # Investigate empirical vs estimated Gamma 2
    some_err_norms = sapply(some_gamma_1_hats, function(x) norm(x - emp_gamma_1, type="2"))

    list_err_norms[[i]] = some_err_norms
}

mean_err_norms = sapply(list_err_norms, mean)
sd_err_norms = sapply(list_err_norms, sd)

(data_err_norms = data.frame(K = all_Ks, err = mean_err_norms, sd = sd_err_norms))

plot(log(all_Ks), log(mean_err_norms))

data_err_norms %>% lm(log(err) ~ log(K), data=.) %>% summary()
data_err_norms %>% lm(log(err) ~ log(K), data=.) %>% abline()

data_err_norms %>% slice(-3) %>% lm(log(err) ~ log(K), data=.) %>% summary()
data_err_norms %>% slice(-3) %>% lm(log(err) ~ log(K), data=.) %>% abline()

expand_REs(which_REs)

q = list_Thetas[[1]]
str(q)



##** Plot values of simulated objects

### Theta
data_Thetas = data.frame()
for(i in seq_along(all_Ks)){
    some_Thetas = list_Thetas[[i]]
    this_data = cbind(all_Ks[i], some_Thetas)

    data_Thetas = rbind(data_Thetas, this_data)
}

colnames(data_Thetas) = c("K", "b_Y_Int", "b_Y_X", "b_Y_M", "b_Y_W1", "b_Y_W2", "theta_Y_Int", "theta_Y_corr", "theta_Y_X", "b_M_Int", "b_M_X", "b_M_W1", "b_M_W2", "theta_M_Int", "theta_M_corr", "theta_M_X")

ggplot(data_Thetas, aes(b_Y_Int=b_Y_Int, b_Y_X=b_Y_X, b_Y_M=b_Y_M, b_Y_W1=b_Y_W1, b_Y_W2=b_Y_W2, theta_Y_Int=theta_Y_Int, theta_Y_corr=theta_Y_corr, theta_Y_X=theta_Y_X, b_M_Int=b_M_Int, b_M_X=b_M_X, b_M_W1=b_M_W1, b_M_W2=b_M_W2, theta_M_Int=theta_M_Int, theta_M_corr=theta_M_corr, theta_M_X=theta_M_X)) + geom_path(alpha=0.1) + coord_serialaxes() + facet_wrap(~K, scales = "free_y")



### ENCs
data_ENCs = data.frame()
for(i in seq_along(all_Ks)){
    some_ENCs = list_ENCs[[i]]
    this_data = cbind(all_Ks[i], some_ENCs)

    data_ENCs = rbind(data_ENCs, this_data)
}

colnames(data_ENCs) = c("K", "V11", "V10", "V01", "V00")

ggplot(data_ENCs, aes(V11=V11, V10=V10, V01=V01, V00=V00)) + geom_path(alpha=0.1) + coord_serialaxes() + facet_wrap(~K, scales = "free_y")




### ENC covariances
data_ENC_covs = data.frame()
data_ENC_cov_diffs = data.frame()
for(i in seq_along(all_Ks)){
    some_ENC_cov_mats = list_ENC_covs[[i]]
    some_ENC_cov_vecs = t(sapply(some_ENC_cov_mats, gdata::upperTriangle, diag=T))
    this_data = cbind(all_Ks[i], some_ENC_cov_vecs)

    data_ENC_covs = rbind(data_ENC_covs, this_data)


    this_emp_cov = cov(list_ENCs[[i]])
    this_data_diff = cbind(all_Ks[i], some_ENC_cov_vecs - gdata::upperTriangle(this_emp_cov, diag=T))

    data_ENC_cov_diffs = rbind(data_ENC_cov_diffs, this_data_diff)
}

colnames(data_ENC_covs) = c("K", paste0("V", 1:10))
colnames(data_ENC_cov_diffs) = c("K", paste0("V", 1:10))

ggplot(data_ENC_covs, aes(V1=V1, V2=V2, V3=V3, V4=V4, V5=V5, V6=V6, V7=V7, V8=V8, V9=V9, V10=V10)) + geom_path(alpha=0.1) + coord_serialaxes() + facet_wrap(~K, scales = "free_y")


#### Add a line at zero for each level of K
data_ENC_cov_diffs %<>% mutate(color_label = "1", opacity = 0.1)

data_zero = data.frame()
for(i in seq_along(all_Ks)){
    this_data = data.frame(K = all_Ks[i], V1 = 0, V2 = 0, V3 = 0, V4 = 0, V5 = 0, V6 = 0, V7 = 0, V8 = 0, V9 = 0, V10 = 0, color_label = "2", opacity = 1)
    data_zero = rbind(data_zero, this_data)
}
data_ENC_cov_diffs = rbind(data_ENC_cov_diffs, data_zero)




ggplot(data_ENC_cov_diffs, aes(V1=V1, V2=V2, V3=V3, V4=V4, V5=V5, V6=V6, V7=V7, V8=V8, V9=V9, V10=V10, color = color_label)) + geom_path(aes(alpha=opacity)) + coord_serialaxes() + facet_wrap(~K, scales = "free_y")
