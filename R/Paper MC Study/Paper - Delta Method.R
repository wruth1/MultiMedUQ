


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



load("R/Paper MC Study/all_datasets.RData", verbose = TRUE)


num_reps = length(all_datasets)

all_ME_hats_delta = list()
all_cov_hats_delta = list()

runtime_delta_specific = 0

tic()
for(i in seq_len(num_reps)){

    print(paste0("i = ", i, " of ", num_reps))

    data = all_datasets[[i]]


    tryCatch({
            
        ## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
        (fit_Y = suppressMessages(lme4::glmer(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))
        (fit_M = suppressMessages(lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))

        
        # Compute parameter estimates
        this_MEs = all_MEs_models(w = w, fit_Y = fit_Y, fit_M = fit_M)

        ## Extract model parameter estimates
        info_Y = get_model_pars(fit_Y)
        info_M = get_model_pars(fit_M)
        this_Theta_hat = c(unlist(info_Y), unlist(info_M))
        this_cov_hat = all_pars_cov_mat(fit_Y, fit_M)

        # Compute covariance matrix
        #? Note: A previous implementation used all_cov_MEs here, which computes the reg par covariance matrix and applies the delta method. This was slow because merDeriv::vcov.glmerMod() is slow (all the code that I wrote runs quickly). I've now separated out these two steps so I can time only my method.
        time = Sys.time()
        this_ME_cov = all_covs_MEs_models(w = w, Sigma = this_cov_hat, fit_Y = fit_Y, fit_M = fit_M)
        this_runtime_specific = Sys.time() - time
        runtime_delta_specific = runtime_delta_specific + this_runtime_specific

        }, error = function(e){
            this_MEs = NULL
            this_ME_cov = NULL
    })


    all_ME_hats_delta[[i]] = this_MEs
    all_cov_hats_delta[[i]] = this_ME_cov


}

time_info = toc()
runtime_delta = time_info$toc - time_info$tic

save(all_ME_hats_delta, all_cov_hats_delta, runtime_delta, file = "R/Paper MC Study/Paper - Delta Method Data.RData")


tic()

toc()
