




# TODO: Figure out why there is a discrepancy between the current implementation and the one from the Exact_Asymptotics project

#* To run the following:
library(lme4)
library(merDeriv)
source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
devtools::load_all()


# Set parameters and fit models
set.seed(1)


x = 0
x_m = 1

b_Y = c(0,0,1,0,0)
# theta_Y = c(sqrt(0.5), 0.5, 0, 1, 0.5, sqrt(0.5))
theta_Y = c(sqrt(0.5), 0, 1)

b_M = c(0,0,0,0)
theta_M = c(1, 0.5, 2)


data = make_validation_data(20, 100, b_Y, theta_Y, b_M, theta_M, output_list = F)

w = c(0,0)

## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
(fit_Y = lme4::glmer(Y ~ X + M + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))
(fit_M = lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))



med_hats = all_MEs_models(scale = "OR", w, fit_Y, fit_M, which_REs = c("Y.Int", "Y.X", "M.All")) 

cov_hat = all_cov_MEs(scale = "OR", w, fit_Y, fit_M, which_REs = c("Y.Int", "Y.X", "M.All"))
