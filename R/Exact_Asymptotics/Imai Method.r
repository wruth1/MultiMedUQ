


library(lme4)
library(merDeriv)
library(tictoc)
library(pbapply)
library(parallel)
library(magrittr)
library(dplyr)
library(kableExtra)
library(ggplot2)
library(MASS)
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




load("Par_Hat_MC-Large_K_Pooled.RData", verbose = TRUE)
load("ENC_Cov_Hats.RData", verbose = TRUE)



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





B = 50000


list_all_theta_errs = list()
for(i in seq_along(all_Ks)){
    some_theta_hats = list_par_hats[[i]]

    theta_bar = apply(some_theta_hats, 2, mean)

    some_errs = c()
    for(j in 1:nrow(some_theta_hats)){
        this_err = some_theta_hats[j,] - theta_bar

        some_errs = rbind(some_errs, norm(this_err, type="2"))
    }

    list_all_theta_errs[[i]] = some_errs
}

data_all_theta_errs = data.frame(K = rep(all_Ks, each = num_reps), err = Reduce("c", list_all_theta_errs))

# histogram of errors, faceted by K
ggplot(data_all_theta_errs, aes(x = err)) + geom_histogram() + facet_wrap(~K)


library(ggmulti)

i=1
q = data.frame(list_par_hats[[i]])
colnames(q) = c(paste0("b_Y", 1:5), paste0("theta_Y", 1:3), paste0("b_M", 1:4), paste0("theta_M", 1:3))
q_sub = q[sample(1:nrow(q), 200),]

ggplot(q_sub, aes(b_Y1 = b_Y1,
                                b_Y2 = b_Y2,
                                b_Y3 = b_Y3,
                                b_Y4 = b_Y4,
                                b_Y5 = b_Y5,
                                theta_Y1 = theta_Y1,
                                theta_Y2 = theta_Y2,
                                theta_Y3 = theta_Y3,
                                b_M1 = b_M1,
                                b_M2 = b_M2,
                                b_M3 = b_M3,
                                b_M4 = b_M4,
                                theta_M1 = theta_M1,
                                theta_M2 = theta_M2,
                                theta_M3 = theta_M3
                                )) + 
    geom_path(alpha=0.1) + coord_serialaxes() + geom_histogram()
 

# MC delta method for ENC

# list_cov_ENC_tildes = list()
# list_all_ENC_tildes = list()

# for(i in seq_along(all_Ks)){

#     this_theta_hat = apply(list_par_hats[[i]], 2, mean)

#     this_cov_hat = list_mean_covs[[i]]


#     some_theta_tildes = mvrnorm(B, mu = this_theta_hat, Sigma = this_cov_hat)

#     some_ENC_tildes = data.frame()

#     for(j in 1:B){

#         if(j%%1000 == 0) cat("i = ", i, " of ", length(all_Ks), ", j = ", j, " of ", B, "\n", sep="")



#         this_theta_tilde = some_theta_tildes[j,]

#         this_b_Y = this_theta_tilde[1:5]
#         this_theta_Y = this_theta_tilde[6:8]
#         this_b_M = this_theta_tilde[9:12]
#         this_theta_M = this_theta_tilde[13:15]

#         this_ENC_tilde = all_ENCs(w, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs)

#         some_ENC_tildes = rbind(some_ENC_tildes, this_ENC_tilde)
#     }


#     list_all_ENC_tildes[[i]] = some_ENC_tildes
#     cov_ENC_tilde = cov(some_ENC_tildes)


#     list_cov_ENC_tildes[[i]] = cov_ENC_tilde
# }


# save(all_Ks, B, list_all_ENC_tildes, list_cov_ENC_tildes, file="ENC_MC_Delta.RData")



data_err_tilde = data.frame()
for(i in seq_along(all_Ks)){
    this_K = all_Ks[i]

    this_emp_cov = list_ENC_emp_covs[[i]]
    this_cov_tilde = list_cov_ENC_tildes[[i]]

    abs_err = norm(this_emp_cov - this_cov_tilde, type="2")
    scaled_err = abs_err * this_K

    this_info = c(this_K, abs_err, scaled_err)
    data_err_tilde = rbind(data_err_tilde, this_info)
}
colnames(data_err_tilde) = c("K", "Abs", "Scaled")
data_err_tilde

# Plot log abs error vs log K

ggplot(data = data_err_tilde, aes(x = log(K), y = log(Abs))) + geom_point() + geom_line()
# ggplot(data = data_err_tilde[c(-4),], aes(x = log(K), y = log(Abs))) + geom_point() + geom_line()

fit_log_tilde = lm(log(Abs) ~ log(K), data = data_err_tilde[c(-4),])
summary(fit_log_tilde)

plot(fit_log_tilde, 1)



# Trajectories of error as B increases
B_seq = seq(1000, B, by = 1000)
norms_by_MC_size = data.frame()

for(i in seq_along(all_Ks)){
    for(r in seq_along(B_seq)){
        this_B = B_seq[r]

        # this_emp_cov = cov(list_par_hats2[[i]][1:this_B,])
        # this_mean_cov = Reduce("+", list_par_cov_hats2[[i]][1:this_B]) / this_B

        # # Pool both simulations' results
        this_emp_cov = list_ENC_emp_covs[[i]]
        this_MC_cov = cov(list_all_ENC_tildes[[i]][1:this_B,])

        this_err = norm(this_emp_cov - this_MC_cov, type = "2")

        this_err_scaled = this_err * all_Ks[i]
        this_err_extra_scaled = this_err_scaled * sqrt(all_Ks[i])

        this_info = c(all_Ks[i], this_B, this_err, this_err_scaled, this_err_extra_scaled)


        norms_by_MC_size = rbind(norms_by_MC_size, this_info)
    }
}
colnames(norms_by_MC_size) = c("K", "B", "Abs", "Scaled", "Extra_Scaled")


# Plot abs error vs B for each K
ggplot(data = norms_by_MC_size, aes(x = B, y = Abs)) + geom_point() + geom_line() + facet_wrap(~K, scales = "free")

# Plot log abs error vs K for largest value of B
norms_by_MC_size %>% 
    filter(B == max(norms_by_MC_size$B)) %>% 
    ggplot(aes(x = K, y = Abs)) + geom_point() + geom_line()

norms_by_MC_size %>% 
    # filter(B == max(norms_by_MC_size$B), K %in% c(100, 400, 800)) %>% 
    filter(B == max(norms_by_MC_size$B)) %>% 
    lm(log(Abs) ~ log(K), data = .) %>% summary(.)

norms_by_MC_size %>% 
    lm(log(Abs) ~ log(K), data = .) %>% summary(.)
