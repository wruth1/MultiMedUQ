


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
source("R/Exact_Asymptotics/Imai Method.r")
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



set.seed(2222)




load("Par_Hat_MC-Large_K_Pooled.RData", verbose = TRUE)
load("ENC_Cov_Hats.RData", verbose = TRUE)



# ---------------------------------------------------------------------------- #
#        Simple MC study comparing MC delta method with my delta method        #
# ---------------------------------------------------------------------------- #


B = 10000
# B = 500

list_cov_ME_tildes = list()
list_all_ME_tildes = list()

for(i in seq_along(all_Ks)){


    print(paste0("K = ", all_Ks[i], "; number ", i, " of ", length(all_Ks)))


    this_Theta_hat = apply(list_par_hats[[i]], 2, mean)

    this_cov_hat = list_mean_covs[[i]]


    some_Theta_tildes = sim_Theta_tildes(B, this_Theta_hat, this_cov_hat)

    some_ME_tildes = Theta_tildes_2_MEs(scale, w, some_Theta_tildes, which_REs)


    list_all_ME_tildes[[i]] = some_ME_tildes

    cov_ME_tilde = cov(some_ME_tildes)
    list_cov_ME_tildes[[i]] = cov_ME_tilde
}




# ---------------- Direct delta method with same initial covs ---------------- #

list_cov_ME_deltas = list()

for(i in seq_along(all_Ks)){


    print(paste0("K = ", all_Ks[i], "; number ", i, " of ", length(all_Ks)))


    this_Theta_hat = apply(list_par_hats[[i]], 2, mean)

    this_cov_hat = list_mean_covs[[i]]


    # The number of parameters in theta_Y based on the number of REs for Y
    len_theta_Y = which_REs %>%
                    num_Y_REs() %>%
                    num_REs2theta_length()

                    
    # Extract parameters
    b_Y = this_Theta_hat[1:5]
    theta_Y = this_Theta_hat[6:(5 + len_theta_Y)]
    b_M = this_Theta_hat[(6 + len_theta_Y):(9 + len_theta_Y)]
    theta_M = this_Theta_hat[(10 + len_theta_Y):length(this_Theta_hat)]

    # Compute covariance matrix for MEs
    cov_ME_delta = all_covs_MEs_pars(scale, w, this_cov_hat, b_Y, theta_Y, b_M, theta_M, which_REs)


    list_cov_ME_deltas[[i]] = cov_ME_delta
}


# -------------------------- Compare both approaches ------------------------- #

MC_vs_delta_errs = c()
for(i in seq_along(all_Ks)){
    this_err = norm(list_cov_ME_deltas[[i]] - list_cov_ME_tildes[[i]], type="2") / norm(list_cov_ME_deltas[[i]], type="2")

    MC_vs_delta_errs = c(MC_vs_delta_errs, this_err)
}
MC_vs_delta_errs





# ---------------------------------------------------------------------------- #
#                     Run MC delta method on full MC study                     #
# ---------------------------------------------------------------------------- #

B_MC_delta = 1000

list_cov_tildes = list()

tic()

for(i in seq_along(list_par_hats)){

    some_Theta_hats = list_par_hats[[i]]
    some_cov_hats = list_par_cov_hats[[i]]

    some_cov_tildes = list()

    for(j in 1:nrow(some_Theta_hats)){

        if(j %% 100 == 0) print(paste0("i = ", i, " of ", length(list_par_hats), ", j = ", j, " of ", nrow(some_Theta_hats)))

        this_Theta_hat = some_Theta_hats[j,]

        this_cov_hat = some_cov_hats[[j]]


        #! Some estimated covariance matrices are not positive definite. In this case, skip the current iteration
        #! This isn't really a great option. All the covariance matrices should be positive definite.
        tryCatch({
            this_cov_tilde = MC_delta(B_MC_delta, this_Theta_hat, this_cov_hat, scale, w, which_REs)
            some_cov_tildes[[j]] = this_cov_tilde

            }, error = function(e){
                    print(paste0("Error at i = ", i, ", j = ", j, ": ", e))
            }
        )
    }

    list_cov_tildes[[i]] = some_cov_tildes
}

toc()

test = list()
test[[1]] = 1
test[[2]] = NA
test[[3]] = 3