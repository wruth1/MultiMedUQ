
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
source("R/Exact_Asymptotics/Imai Method.r")
devtools::load_all()





# Set parameters and fit models

# N = 100
# N = 20
# N = 40
# N = 60
n = 50

# all_Ks = c(50, 100, 200, 400, 800)
# all_Ks = c(50, 100, 200)
# all_Ks = 50 * (2:6)
K = 200

# num_reps = 30
# num_reps = 500
num_reps = 100


# which_REs = c("Y.Int", "Y.X", "M.All")
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")


x = 0
x_m = 1



w = c(2,3)


## Non-trivial values for the b's and theta's. Former based on output from another MC study. Latter chosen arbitrarily.
## Crucially, no parameters are equal to zero.
##? We choose the intercepts to that the mean of the linear predictor is zero. Doing this for M makes it easier to do so for Y.
# b_Y_int_old = 0.0376828219852018
b_Y_X = 0.966486302988689
b_Y_M = 1.99644760563721
b_Y_C1 = -0.00556557712859059
b_Y_C2 = 0.000826754128449799
b_Y_int = - sum(b_Y_X, b_Y_M, b_Y_C1, b_Y_C2) / 2       # ~ -1.48
b_Y = c(b_Y_int, b_Y_X, b_Y_M, b_Y_C1, b_Y_C2)
# b_Y = c(0.0376828219852018, 0.966486302988689, 1.99644760563721, -0.00556557712859059, 0.000826754128449799)

# b_M_int_old = -0.0990439890654785
b_M_X = 1.76353928991247
b_M_C1 = 0.0128566136999183
b_M_C2 = 0.00711746366915989
b_M_int = -sum(b_M_X, b_M_C1, b_M_C2) / 2       # ~ -0.89
b_M = c(b_M_int, b_M_X, b_M_C1, b_M_C2)
# b_M = c(-0.0990439890654785, 1.76353928991247, 0.0128566136999183, 0.00711746366915989)




# Choose theta_Y and theta_M based on the values of b_Y and b_M
sigma = 1 / 10
rho = 0.5
theta_Y = c(sigma, rho, rho, sigma, rho, sigma)
theta_M = c(sigma, rho, sigma)
# theta_Y = c(sqrt(0.5), 0.5, 0.5, 1, 0.5, sqrt(0.5)) / 10
# # theta_Y = c(sqrt(0.5), 0.5, 1)
# # theta_M = c(1, 0.5, 2)
# theta_M = c(sqrt(0.5), 0.5, 1) / 10


all_reg_pars = c(b_Y, theta_Y, b_M, theta_M)




# set.seed(1)

# all_datasets = pblapply(seq_len(num_reps), function(i){
#     make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
# })


set.seed(123456)

all_datasets = pblapply(seq_len(num_reps), function(i){
    make_validation_data(n, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
})
save(all_datasets, file = "R/Paper MC Study/all_datasets (fitted pars).RData")

for(i in seq_along(all_datasets)){
    if(i %% 10 == 0) print(paste0("i = ", i, " of ", length(all_datasets)))
    data = all_datasets[[i]]
    save(data, file = paste0("R/Paper MC Study/Datasets - Fitted Pars/", i, ".RData"))
}


w = c(2,3)
B = 500
scale = c("diff", "rat", "OR")
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")



# Setup cluster
# cl = makeCluster(detectCores() - 2)
# cl = makeCluster(15)
cl = makeCluster(10)
# clusterExport(cl, c("N", "b_Y", "theta_Y", "b_M", "theta_M", "which_REs"))
clusterExport(cl, c("w", "B", "scale", "which_REs")) #, "n", "K"))
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
    source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
    source("R/Exact_Asymptotics/Imai Method.r")
    devtools::load_all()
})
clusterSetRNGStream(cl = cl, 123)



# all_ME_hats = list()
# all_cov_hats_delta = list()
# all_cov_hats_MC_delta = list()

# total_runtime_delta = 0
# total_runtime_MC_delta = 0

# #? Note: Extracting the SE matrix for the reg pars is slow (using merDeriv::vcov.glmerMod()). I don't want to duplicate this step (or fitting the models), so I separated out model fitting/SE extraction from my method and the MC delta.


MC_results_delta_MC_delta = pblapply(1:num_reps, function(i) {
    load(paste0("R/Paper MC Study/Datasets - Fitted Pars/", i, ".RData"))

    # #! Remove half the groups
    # groups_keep = paste0("G", 1:(K / 2))
    # data %<>% dplyr::filter(group %in% groups_keep)

    # #! Remove half the samples in each group
    # data %<>% dplyr::group_by(group) %>%
    #     dplyr::slice_sample(n = n / 2) %>% dplyr::ungroup()


    tryCatch({

        ## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
        (fit_Y = suppressMessages(lme4::glmer(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))
        (fit_M = suppressMessages(lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))

        pars_Y = get_model_pars(fit_Y, format = "vector")
        pars_M = get_model_pars(fit_M, format = "vector")

        this_pars = c(pars_Y, pars_M)
        save(this_pars, file = paste0("R/Paper MC Study/Results - Fitted Pars/", i, ".RData"))


        }, error = function(e){
            this_pars = NULL
            save(this_pars, file = paste0("R/Paper MC Study/Results - Fitted Pars/", i, ".RData"))

    })

    # stop("Error: Loop should never reach this point.")

}, cl = cl)

stopCluster(cl)


#* Compile results into a single data frame

all_par_hats = t(pbsapply(1:100, function(i) {
    load(paste0("R/Paper MC Study/Results - Fitted Pars/", i, ".RData"))
    return(this_pars)
}))

all_par_hat_means = colSums(all_par_hats) / nrow(all_par_hats)
all_par_hat_SDs = apply(all_par_hats, 2, sd)
all_par_hat_SEs = all_par_hat_SDs / sqrt(nrow(all_par_hats))

(all_par_hat_means - all_reg_pars) / all_reg_pars

K_500_n_50 = data.frame(mean = all_par_hat_means, true = all_reg_pars, SE = all_par_hat_SEs)
rel_500 = K_500_n_50 %>% mutate(rel_err = (mean - true) / true) %>% select(rel_err)

K_1000_n_50 = data.frame(mean = all_par_hat_means, true = all_reg_pars, SE = all_par_hat_SEs)
K_1000_n_50 %>% mutate(rel_err = (mean - true) / true)
rel_1000 = K_1000_n_50 %>% mutate(rel_err = (mean - true) / true) %>% select(rel_err)

data_rel_errs = abs(cbind(rel_500, rel_1000))
apply(data_rel_errs, 2, formatC, format = "e", digits = 3)


apply(all_par_hats, 2, function(x) sum(is.na(x)))
nrow(na.omit(all_par_hats))
