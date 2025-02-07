


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
library(broom.mixed)
library(glmmTMB)
library(stringr)
source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.r")
source("R/Exact_Asymptotics/Imai Method.r")
devtools::load_all()


#? I denote regression coefficients for Y as b_Y and for M as b_M
#? Random effect parameters for Y are theta_Y and for M as theta_M
    #? These RE parameters are SDs and correlations. They are organized as, e.g. with 3 REs, SD_1, corr_12, corr_13, SD_2, corr_23, SD_3
#? We denote the single vector containing all parameters by Theta = (c(b_Y, theta_Y, b_M, theta_M))
    #? I.e. Capital Theta contains both small theta vectors

# Set parameters

B = 500
scale = c("diff", "rat", "OR")
which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")



load("R/Trust_Study/Trust_Study_Data.RData", verbose = T)


dat.ma %<>% filter(q8.pcis %in% c("low", "sure"))


#* Fit models
fit_Y = glmmTMB(q5.fc ~ q8.pcis + q4.src + q1.cc  + q2.dc + q3.pc + age_group + gender + q7.la  + q9.edu + (q8.pcis + q4.src |country), data = dat.ma, family = "binomial")
fit_M = glmmTMB(q4.src ~ q8.pcis + q1.cc + q2.dc + q3.pc + age_group + gender + q7.la + q9.edu + (q8.pcis |country), data = dat.ma, family = "binomial")







#* Set value for confounders
#? Later, I will need a vector of confounder values. For now, I'm going to just use the first row of the design matrix. 
#ToDo: Check how sensitive the mediation effects are to this choice. Try other rows.

confounder_row_number = 2
w = broom.helpers::model_get_model_matrix(fit_M) %>% 
        as_tibble() %>%
        dplyr::slice(confounder_row_number) %>%
        select(-contains("Intercept"), -contains("q8"), -contains("country")) %>%
        as.matrix()









#* Extract parameters
theta_hat_Y = get_model_pars_TMB(fit_Y)
theta_hat_M = get_model_pars_TMB(fit_M)
Theta_hat = c(unlist(theta_hat_Y), unlist(theta_hat_M))

b_Y = theta_hat_Y[["b"]]
theta_Y = theta_hat_Y[["theta"]]
b_M = theta_hat_M[["b"]]
theta_M = theta_hat_M[["theta"]]
len_par_vecs = sapply(list(b_Y, theta_Y, b_M, theta_M), length)

## Fitted covariance matrices
all_cov_hats = all_pars_cov_mat_TMB(fit_Y, fit_M, return_separate_covs = TRUE)
cov_hat = all_cov_hats$joint_cov
cov_hat_Y = all_cov_hats$Y_cov
cov_hat_M = all_cov_hats$M_cov


#* Mediation Effects
MEs = all_MEs_pars(scale, w, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)
cov_MEs_delta = all_covs_MEs_pars(scale, w, cov_hat, b_Y, theta_Y, b_M, theta_M, which_REs =  which_REs)

#* Construct CIs - Our Method
CIs_raw = build_one_CI_many_pars(MEs, sqrt(diag(cov_MEs_delta)))
CIs_delta = sapply(CIs_raw, function(x) c(x$lcl, x$ucl)) %>% t()
rownames(CIs_delta) = names(MEs)
colnames(CIs_delta) = c("lcl", "ucl")
CIs_delta



# ---------------------------------------------------------------------------- #
#                                   MC Delta                                   #
# ---------------------------------------------------------------------------- #
set.seed(1)
B = 500

# Generate some sample parameter values, Theta_tilde_1, ..., Theta_tilde_B
some_Theta_tildes = sim_TMB_Theta_tildes(B, fit_Y, fit_M, cov_hat_Y = cov_hat_Y, cov_hat_M = cov_hat_M)
# Compute mediation effects for each simulated set of parameter values
some_ME_tildes = Theta_tildes_2_MEs(scale, w, some_Theta_tildes, which_REs, len_par_vecs = len_par_vecs)
# Estimate covariance matrix of sampling distribution for mediation effects
cov_MEs_MC_delta = cov(some_ME_tildes)

#* Construct CIs - MC Delta
CIs_raw = build_one_CI_many_pars(MEs, sqrt(diag(cov_MEs_MC_delta)))
CIs_MC_delta = sapply(CIs_raw, function(x) c(x$lcl, x$ucl)) %>% t()
rownames(CIs_MC_delta) = names(MEs)
colnames(CIs_MC_delta) = c("lcl", "ucl")
CIs_MC_delta




# ---------------------------------------------------------------------------- #
#                                Plot Intervals                                #
# ---------------------------------------------------------------------------- #

CIs_delta_plot = CIs_delta %>% as.data.frame() %>% 
                    mutate(method = "Delta", estimate = MEs, effect = names(MEs))
CIs_MC_delta_plot = CIs_MC_delta %>% as.data.frame() %>% 
                    mutate(method = "MC-Delta", estimate = MEs, effect = names(MEs))
CIs_plot = rbind(CIs_delta_plot, CIs_MC_delta_plot) %>% 
                mutate(scale = str_extract(effect, "diff|rat|OR"), 
                        effect = str_extract(effect, "total|direct|indirect"),
                        effect_method = paste(effect, method, sep = "_"))


this_scale = "diff"
# this_scale = "rat"
# this_scale = "OR"

ref_val = switch(this_scale, diff = 0, rat = 1, OR = 1)
scale_name = switch(this_scale, diff = "Difference", rat = "Ratio", OR = "Odds-Ratio")

pdf(paste0("R/Trust_Study/CI_Plot-", this_scale, ".pdf"), width=10, height=7)

filter(CIs_plot, scale == this_scale) %>%
ggplot(aes(y = effect_method, x = estimate, xmin = lcl, xmax = ucl, color = method)) +
  geom_pointrange() + ylab("Effect") + xlab("Estimate") +
  geom_vline(xintercept = ref_val, linetype = "dashed") +
  scale_y_discrete(labels = rep(c("Total", "Indirect", "Direct"), each = 2)) +
  theme(text = element_text(size = 30), plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste("95% CIs (", scale_name, " Scale)", sep = ""))



dev.off()