
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
N = 100
n = N

# all_Ks = c(50, 100, 200, 400, 800)
# all_Ks = c(50, 100, 200)
# all_Ks = 50 * (2:6)
K = 200

# num_reps = 30
# num_reps = 500
num_reps = 1000


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
theta_Y = c(sqrt(0.5), 0.5, 0.5, 1, 0.5, sqrt(0.5)) / 2
# theta_Y = c(sqrt(0.5), 0.5, 1)
# theta_M = c(1, 0.5, 2)
theta_M = c(sqrt(0.5), 0.5, 1) / 2


all_reg_pars = c(b_Y, theta_Y, b_M, theta_M)




# set.seed(1)

# all_datasets = pblapply(seq_len(num_reps), function(i){
#     make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
# })


set.seed(123456)

more_datasets = pblapply(seq_len(num_reps), function(i){
    make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
})
load("R/Paper MC Study/all_datasets.RData", verbose = TRUE)
all_datasets = c(all_datasets, more_datasets)
save(all_datasets, file = "R/Paper MC Study/all_datasets.RData")

q = c(all_datasets, more_datasets)

# all_datasets = list()

# for (i in 1:num_reps) {
#     print(paste0("i = ", i, " of ", num_reps))
#     all_datasets[[i]] = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs)
# }

# save(all_datasets, file = "R/Paper MC Study/all_datasets.RData")
load("R/Paper MC Study/all_datasets.RData", verbose = TRUE)

for(i in seq_along(all_datasets)){
    if(i %% 100 == 0) print(paste0("i = ", i, " of ", length(all_datasets)))
    data = all_datasets[[i]]
    save(data, file = paste0("R/Paper MC Study/Datasets/", i, ".RData"))
}





# Compute true values of mediation effects
scale = c("diff", "rat", "OR")
true_MEs = all_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M, which_REs = which_REs)

save(true_MEs, file = "R/Paper MC Study/true_MEs (new).RData")
