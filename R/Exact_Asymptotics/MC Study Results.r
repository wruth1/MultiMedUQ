

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
# source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
devtools::load_all()



load(file = "Par_Hat_MC_Clean.RData", verbose = TRUE)
load(file = "ENC_Hat_MC_Clean.RData", verbose = TRUE)
load(file = "ME_Hat_MC_Clean.RData", verbose = TRUE)




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




#* Plot trajectories wrt MC size

