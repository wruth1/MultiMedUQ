


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

B = 500

load("R/Paper MC Study/all_datasets.RData", verbose = TRUE)


num_reps = length(all_datasets)

all_ME_hats_MC_delta = list()
all_cov_hats_MC_delta = list()

tic()
for(i in seq_len(num_reps)){

    print(paste0("i = ", i, " of ", num_reps))

    data = all_datasets[[i]]


    tryCatch({
            
        ## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
        (fit_Y = suppressMessages(lme4::glmer(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))
        (fit_M = suppressMessages(lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))

        # Estimate mediation effects
        this_MEs = all_MEs_models(w = w, fit_Y = fit_Y, fit_M = fit_M)

        
        # Estimate SE of mediation effects

        ## Extract model parameter estimates
        info_Y = get_model_pars(fit_Y)
        info_M = get_model_pars(fit_M)
        this_Theta_hat = c(unlist(info_Y), unlist(info_M))
        this_cov_hat = all_pars_cov_mat(fit_Y, fit_M)

        some_Theta_tildes = sim_Theta_tildes(B, this_Theta_hat, this_cov_hat)
        some_ME_tildes = Theta_tildes_2_MEs(scale, w, some_Theta_tildes, which_REs)
        cov_ME_tilde = cov(some_ME_tildes)


        }, error = function(e){
            this_MEs = NULL
            this_ME_cov = NULL
    })


    all_ME_hats_MC_delta[[i]] = this_MEs
    all_cov_hats_MC_delta[[i]] = cov_ME_tilde


}

time_info = toc()
runtime_MC_delta = time_info$toc - time_info$tic

save(all_ME_hats_MC_delta, all_cov_hats_MC_delta, runtime_MC_delta, file = "R/Paper MC Study/Paper - Quasi Bayesian Data.RData")