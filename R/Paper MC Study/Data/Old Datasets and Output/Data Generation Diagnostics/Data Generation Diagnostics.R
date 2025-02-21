
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
N = 1000
n = N

# all_Ks = c(50, 100, 200, 400, 800)
# all_Ks = c(50, 100, 200)
# all_Ks = 50 * (2:6)
# K = 200
K = 1000

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



data_and_REs = make_validation_data(N, K, b_Y, theta_Y, b_M, theta_M, output_list = F, which_REs = which_REs, return_REs = T)

data = data_and_REs[["data"]]
all_REs = data_and_REs[["all_REs"]]




lin_preds_from_vecs = function(data, b_fix, inds_fix, b_ran, inds_ran){
    p = ncol(data)
    b = rep(0, times = p)

    for(i in inds_fix){
        b[i] = b[i] + b_fix[i]
    }

    for(i in inds_ran){
        b[i] = b[i] + b_ran[i]
    }

    lin_preds = lin_pred_contrib(data, b, add_intercept = F)

    return(lin_preds)
    
}




#! WARNING: The order of groups in data is not the same as the order of groups in all_REs

all_E_Ms = c()
all_M_bars = c()

for(i in seq_along(all_REs)){
    if (i %% 50 == 0) print(paste0("i = ", i, " of ", length(all_REs)))

    # i=1
    this_group = all_REs[[i]][["data"]]
    # this_group = filter(data, group == paste0("G", i))
    this_REs = all_REs[[i]][["M"]]



    this_data_M = this_group %>%
                    select(-Y, -M) %>% 
                    cbind(1, .)

    lin_preds_M = lin_preds_from_vecs(this_data_M, b_M, seq_len(ncol(this_data_M)), this_REs, c(1, 2))




    #* Compute conditional expectation of M given U

    data_all_combs = cbind(1, expand.grid(c(0,1), c(0,1), c(0,1)))
    lin_pred_all_combs = lin_preds_from_vecs(data_all_combs, b_M, seq_len(ncol(data_all_combs)), this_REs, c(1,2))
    probs_all_combs = expit(lin_pred_all_combs)

    this_E_M = mean(probs_all_combs)
    all_E_Ms[i] = this_E_M


    #* Estimate conditional expectation of M given U using MC sample
    this_M_bar = mean(this_group$M)
    all_M_bars[i] = this_M_bar
    # this_M_SD = sd(this_group$M)
}



all_E_Ms
all_M_bars

rel_errs = abs(all_M_bars - all_E_Ms)/ all_E_Ms
mean(rel_errs)
sd(rel_errs)
max(rel_errs)

cor(all_E_Ms, all_M_bars)


#* Verify that conditional variances are correct
#* I.e. V(M|U) = E(M|U) * (1 - E(M|U))


all_emp_vars = sapply(all_REs, function(x) {
    this_data = x[["data"]]
    return(var(this_data$M))
})

all_cond_vars = all_E_Ms * (1 - all_E_Ms)


mean(abs(all_emp_vars - all_cond_vars) / all_cond_vars)
cor(all_emp_vars, all_cond_vars)




