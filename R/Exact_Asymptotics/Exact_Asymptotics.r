
# Note: This code chunk is not run here, but it is called and run at the start of every subsequent section. Any changes made here will be duplicated throughout the rest of the document.

# set.seed(1)
set.seed(12345)

num_reps = 1000

all_Ns = c(100, 200, 500, 1000, 2000, 5000, 10000)
p_conf = 3  # Number of confounders


# Regression coefficients
## Intercepts are set separately for each model so that the mean of the linear predictor is approximately zero

## M model
a_1 = 1                       # Coefficient for X
A_2 = rep(1, times = p_conf)  # Coefficients for confounders

# Y model
b_1 = 1                       # Coefficient for M
b_2 = 1                       # Coefficient for X
B_3 = rep(1, times = p_conf)  # Coefficients for confounders


# Containers for output
all_a_hats = list()
all_b_hats = list()

all_a_SEs = list()
all_b_SEs = list()

all_med_hats = list()
all_med_SEs = list()








source("R/Exact_Asymptotics/Exact_Asymptotics_Helpers.R")

library(lme4)     # For fitting GLMMs using glmer()
library(MASS)     # For simulating multivariate normals using mvrnorm()
library(merDeriv) # For computing the full information matrix from models fit using lme4
library(tictoc)
library(ggplot2)
library(magrittr)
library(dplyr)

### Easier version of problem ###

all_Ns = c(50, 200, 500, 1000)  # Sample size for each cluster
all_Ks = c(5, 20, 50)           # Number of clusters
N_K_combs = expand.grid(N = all_Ns, K = all_Ks)

# Number of replicates for each sample size
num_reps = 200

# Values of X and W for which we compute the total effect
x_pred = 0
W_pred = c(1,1,1)

# Fixed effects for intercepts
a_0 = 0
a = c(a_0, a_1, A_2)
b_0 = -0.5
b = c(b_0, b_1, b_2, B_3)

# Random effects
sigma_a_0 = 0.5
sigma_a_1 = 0.5 * abs(a_1)
cor_a0_a1 = 0.5
cov_a0_a1 = sigma_a_0 * sigma_a_1 * cor_a0_a1
theta = c(sigma_a_0, sigma_a_1, cor_a0_a1)
Sigma_a = matrix(c(sigma_a_0^2, cov_a0_a1, cov_a0_a1, sigma_a_1^2), nrow = 2, ncol = 2)

sigma_b_0 = 0.5 * abs(b_0)
sigma_b_1 = 0.5 * abs(b_1)
cor_b0_b1 = 0.5
cov_b0_b1 = sigma_b_0 * sigma_b_1 * cor_b0_b1
gamma = c(sigma_b_0, sigma_b_1, cor_b0_b1)
Sigma_b = matrix(c(sigma_b_0^2, cov_b0_b1, cov_b0_b1, sigma_b_1^2), nrow = 2, ncol = 2)




### Harder version ###

# # Fewer large sample sizes
# # all_Ns = c(100, 200, 500, 1000, 2000)
# # all_Ns = c(200, 500, 1000, 2000)
# # all_Ns = c(200, 500, 1000)
# # all_Ns = c(200, 500, 1000)
# 
# Values of X and W for which we compute the total effect
# x_pred = 0
# W_pred = c(1,1,1)
# 
# # Fixed effects for intercepts
# a_0 = 0
# a = c(a_0, a_1, A_2)
# b_0 = -0.5
# b = c(b_0, b_1, b_2, B_3)
# 
# # # Random effects
# # sigma_a_0 = 0.2
# # sigma_a_1 = 0.2 * abs(a_1)
# # cor_a0_a1 = 0.2
# # cov_a0_a1 = sigma_a_0 * sigma_a_1 * cor_a0_a1
# # theta = c(sigma_a_0, sigma_a_1, cor_a0_a1)
# # Sigma_a = matrix(c(sigma_a_0^2, cov_a0_a1, cov_a0_a1, sigma_a_1^2), nrow = 2, ncol = 2)
# # 
# # sigma_b_0 = 0.2 * abs(b_0)
# # sigma_b_1 = 0.2 * abs(b_1)
# # cor_b0_b1 = 0.2
# # cov_b0_b1 = sigma_b_0 * sigma_b_1 * cor_b0_b1
# # gamma = c(sigma_b_0, sigma_b_1, cor_b0_b1)
# # Sigma_b = matrix(c(sigma_b_0^2, cov_b0_b1, cov_b0_b1, sigma_b_1^2), nrow = 2, ncol = 2)
# 
# # Number of clusters
# # all_Ks = c(5, 20, 50, 100)
# # K=5
# 
# # N_K_combs = expand.grid(N = all_Ns, K = all_Ks)
# 
# # Number of replicates for each sample size
# # We use fewer here because model fitting takes longer
# num_reps = 200





# Commented out to avoid lengthy runtime

library(doParallel)
library(pbapply)

# tic()
# 
# set.seed(1)
# 
# # Containers for output
# all_a_hats = list()
# all_theta_hats = list()
# all_b_hats = list()
# all_gamma_hats = list()
# 
# all_M_covs = list()
# all_Y_covs = list()
# 
# all_med_hats = list()
# all_med_SEs = list()
# 
# 
# # Initialize cluster on my machine ----
# # n_cores = 2
# n_cores = 10
# # n_cores = parallel::detectCores() - 1
# my_cluster = makeCluster(n_cores, type = "PSOCK")
# registerDoParallel(my_cluster)
# clusterEvalQ(my_cluster,{
#   library(MASS)
#   library(lme4)
#   library(merDeriv)
#   source("../src/Helpers.R")
# })
# clusterExport(my_cluster, c("N_K_combs", "num_reps", "a_0", "a_1", "A_2", "a", "b_0", "b_1", "b_2", "B_3", "b", "Sigma_a", "theta", "Sigma_b", "gamma", "p_conf", "x_pred", "W_pred"))
# clusterSetRNGStream(my_cluster, iseed = 1)
# 
# # output_bin_bin_ran = pblapply(seq_len(nrow(N_K_combs)), function(j){
# # output_bin_bin_ran = pblapply(nrow(N_K_combs), function(j){
# for(j in seq_len(nrow(N_K_combs))){
#   n = N_K_combs[j, "N"]
#   K = N_K_combs[j, "K"]
# 
#   clusterExport(my_cluster, c("n", "K"))
# 
#   print(paste0("j = ", j, " of ", nrow(N_K_combs),  " (n=", n, ", K=", K, ")"))
# 
# 
#   # # Containers for output with this value of n
#   # some_a_hats = list()
#   # some_theta_hats = list()
#   # some_b_hats = list()
#   # some_gamma_hats = list()
#   #
#   # some_M_covs = list()
#   # some_Y_covs = list()
#   #
#   # some_med_hats = list()
#   # some_med_SEs = list()
# 
#   # for(i in 1:num_reps){
#   some_output = pblapply(seq_len(num_reps), function(i){
#     # print(paste0("n = ", n, ", K = ", K, ", i = ", i, " of ", num_reps))
# 
#     this_output = NULL
# 
#     # Use tryCatch to omit any datasets leading to convergence issues
#     tryCatch({                                                                    ############## Un-quote
#       all_Xs = list()
#       all_Ws = list()
# 
#       for(k in 1:K){
#         X = rnorm(n, mean=0, sd=1)
#         W = matrix(rnorm(n*p_conf, mean=0, sd=1), nrow = n, ncol = p_conf)
# 
#         all_Xs[[k]] = X
#         all_Ws[[k]] = W
#       }
# 
#       # Generate M
#       all_Ms = list()
#       for(k in 1:K){
#         eta_vec_fixed = a_0 + a_1*all_Xs[[k]] + all_Ws[[k]]%*%A_2
# 
#         ## Add random effects
#         a_ran = mvrnorm(1, mu = rep(0, 2), Sigma = Sigma_a)
#         eta_vec = eta_vec_fixed + a_ran[1] + a_ran[2]*all_Xs[[k]]
# 
#         ## Generate M
#         p_M_vec = expit(eta_vec)
#         M = rbinom(n, size = 1, prob = p_M_vec)
#         all_Ms[[k]] = M
#       }
# 
#       # Generate Y
#       all_Ys = list()
#       for(k in 1:K){
#         zeta_vec_fixed = b_0 + b_1*all_Ms[[k]] + b_2 * all_Xs[[k]] + all_Ws[[k]]%*%B_3
# 
#         ## Add random effects
#         b_ran = mvrnorm(1, mu = rep(0, 2), Sigma = Sigma_b)
#         zeta_vec = zeta_vec_fixed + b_ran[1] + b_ran[2]*all_Xs[[k]]
# 
#         ## Generate Y
#         p_Y_vec = expit(zeta_vec)
#         Y = rbinom(n, size = 1, prob = p_Y_vec)
#         all_Ys[[k]] = Y
#       }
# 
# 
#       # Consolidate groups
#       X = do.call(c, all_Xs)
#       W = do.call(rbind, all_Ws)
#       M = do.call(c, all_Ms)
#       Y = do.call(c, all_Ys)
#       group = rep(1:K, each = n)
# 
# 
#       # Fit models
# 
#       ## M
#       M_data = data.frame(M, X, W1 = W[, 1], W2 = W[, 2], W3 = W[,3], group = group)
#       M_model = glmer(M ~ X + W1 + W2 + W3 + (1 + X | group), data = M_data, family = binomial(link = "logit"))
#       M_model_info = attributes(VarCorr(M_model)$group)
# 
#       ### Fitted parameters
#       a_hat = fixef(M_model)
#       a_RE_sds = M_model_info$stddev
#       a_RE_cor = M_model_info$correlation[2,1]
#       theta_hat = c(a_RE_sds, a_RE_cor)
#       if(any(is.nan(theta_hat))) stop("NaNs in theta_hat")  # Skip rest of current analysis if correlation is 0/0
# 
#       ### Estimated SE covariance matrix
#       M_cov = vcov(M_model, full=TRUE, ranpar="sd")
# 
#       # some_a_hats[[i]] = a_hat
#       # some_theta_hats[[i]] = theta_hat
#       # some_M_covs[[i]] = M_cov
# 
# 
#       ## Y
#       Y_data = data.frame(Y, M, X, W1 = W[, 1], W2 = W[, 2], W3 = W[,3], group = group)
#       Y_model = glmer(Y ~ M + X + W1 + W2 + W3 + (1 + X | group), data = Y_data, family = binomial(link = "logit"))
#       Y_model_info = attributes(VarCorr(Y_model)$group)
# 
#       ### Fitted parameters
#       b_hat = fixef(Y_model)
#       b_RE_sds = Y_model_info$stddev
#       b_RE_cor = Y_model_info$correlation[2,1]
#       gamma_hat = c(b_RE_sds, b_RE_cor)
#       if(any(is.nan(gamma_hat))) stop("NaNs in gamma_hat")  # Skip rest of current analysis if correlation is 0/0
# 
#       ### Estimated SE covariance matrix
#       Y_cov = vcov(Y_model, full=TRUE, ranpar="sd")
# 
#       # some_b_hats[[i]] = b_hat
#       # some_gamma_hats[[i]] = gamma_hat
#       # some_Y_covs[[i]] = Y_cov
# 
# 
# 
#       # Estimate mediation effect
# 
#       ## Fixed-effects
# 
#       a_0_hat = a_hat[1]
#       a_x_hat = a_hat[2]
#       A_2_hat = a_hat[3:5]
# 
#       b_0_hat = b_hat[1]
#       b_m_hat = b_hat[2]
#       b_x_hat = b_hat[3]
#       B_3_hat = b_hat[4:6]
# 
# 
#       ## Linear predictors
#       eta_hat = as.numeric(a_0_hat + a_x_hat*x_pred + W_pred%*%A_2_hat)
#       zeta_hat = as.numeric(b_0_hat + b_x_hat * x_pred + W_pred %*% B_3_hat)
# 
# 
#       ## Random effects covariances
#       s_M_0 = a_RE_sds[1]
#       s_M_x = a_RE_sds[2]
#       rho_M = a_RE_cor
# 
#       s_Y_0 = b_RE_sds[1]
#       s_Y_x = b_RE_sds[2]
#       rho_Y = b_RE_cor
# 
# 
#       ## Sigma functions
#       sigma_M1 = sigma_fun(x_pred, s_M_0, s_M_x, rho_M)
#       sigma_M2 = sigma_fun(x_pred + 1, s_M_x, s_M_0, rho_M)
# 
#       sigma_Y1 = sigma_fun(x_pred, s_Y_0, s_Y_x, rho_Y)
#       sigma_Y2 = sigma_fun(x_pred + 1, s_Y_x, s_Y_0, rho_Y)
# 
#       ## Mediation effect
#       ### See Helpers.R for the function Phi, which computes the mediation effect on odds-ratio scale
#       med_hat = Phi(eta_hat, zeta_hat, a_x_hat, b_m_hat, b_x_hat, sigma_M2, sigma_Y2, sigma_M1, sigma_Y1)
#       # some_med_hats[[i]] = med_hat
# 
# 
# 
# 
# 
#       # Estimate SE
# 
#       ## Gradient of OR wrt regression coefficients
#       grad_Phi_obs = grad_Phi(eta_hat, zeta_hat, a_x_hat, b_m_hat, b_x_hat, sigma_M2, sigma_Y2, sigma_M1, sigma_Y1)
#       grad_xi_obs = grad_xi(a_hat, theta_hat, b_hat, gamma_hat, x_pred, W_pred)
#       d_Phi_d_GLMM_pars = grad_Phi_obs %*% grad_xi_obs
# 
#       ## Get asymptotic SE using delta method
# 
#       ### Build joint covariance matrix of regression coefficients
#       M_length = nrow(M_cov)
#       Y_length = nrow(Y_cov)
#       joint_cov = matrix(0, nrow = M_length + Y_length, ncol = M_length +
#       Y_length)
#       joint_cov[1:M_length, 1:M_length] = M_cov
#       joint_cov[(M_length + 1):(M_length + Y_length), (M_length + 1):(M_length +
#       Y_length)] = Y_cov
# 
#       ### Convert to asymptotic covariance matrix
#       asymp_reg_cov = n * joint_cov
# 
#       ### Pre- and post-multiply asymptotic covariance by gradient of Phi wrt GLMM parameters
#       med_asymp_var = d_Phi_d_GLMM_pars %*% asymp_reg_cov %*% t(d_Phi_d_GLMM_pars)
# 
#       ### Get small-sample standard error
#       med_asymp_SE = sqrt(med_asymp_var)
#       med_SE = med_asymp_SE / sqrt(n)
# 
#       # some_med_SEs[[i]] = med_SE
# 
#       this_output = list(a_hat = a_hat, theta_hat = theta_hat, b_hat = b_hat, gamma_hat = gamma_hat, M_cov = M_cov, Y_cov = Y_cov, med_hat = med_hat, med_SE = med_SE)
#     # }, error = function(e){grad_Phi_obs <<- e}) # For troubleshooting parallel errors
#     }, error = function(e){})                                                       ############ Un-quote
# 
#     return(this_output)
# 
# 
#   }, cl = my_cluster)
#   # })
# 
#   # Store output for current values of n and K
# 
#   all_a_hats[[j]] = lapply(some_output, function(x) x$a_hat)
#   all_theta_hats[[j]] = lapply(some_output, function(x) x$theta_hat)
#   all_b_hats[[j]] = lapply(some_output, function(x) x$b_hat)
#   all_gamma_hats[[j]] = lapply(some_output, function(x) x$gamma_hat)
#   all_M_covs[[j]] = lapply(some_output, function(x) x$M_cov)
#   all_Y_covs[[j]] = lapply(some_output, function(x) x$Y_cov)
#   all_med_hats[[j]] = lapply(some_output, function(x) x$med_hat)
#   all_med_SEs[[j]] = lapply(some_output, function(x) x$med_SE)
# }
# 
# 
# runtime = toc()
# 
# parallel::stopCluster(my_cluster)


# Save output
# save(N_K_combs, all_a_hats, all_theta_hats, all_b_hats, all_gamma_hats, all_M_covs, all_Y_covs, all_med_hats, all_med_SEs, runtime, file = "William-Analysis_cache/MC-bin-bin-mix.RData")

load("R/Exact_Asymptotics/MC-bin-bin-mix.RData", verbose = TRUE)



# Re-format output to be more useful. 
data_med_hats = lapply(all_med_hats, unlist)
data_med_SEs = all_med_SEs %>% lapply(unlist) %>% lapply(na.omit)

# SD of estimates across Monte Carlo samples
emp_SEs = sapply(data_med_hats, sd)

# Summaries of estimated SEs (mean and median)
mean_delta_SEs = sapply(data_med_SEs, mean)
median_delta_SEs = sapply(data_med_SEs, median)

# Combine into a table
results = data.frame(n = N_K_combs$N, K = N_K_combs$K, Empirical = emp_SEs,
                     Mean = mean_delta_SEs, Median = median_delta_SEs,
                     "Percent_Err_Mean" = 100*(mean_delta_SEs - emp_SEs)/emp_SEs,
                     "Percent_Err_Median" = 100*(median_delta_SEs - emp_SEs)/emp_SEs)


# Add extra values of n and K to list
# all_N_combs = c(N_K_combs$N)
# all_K_combs = c(N_K_combs$K)
# results = data.frame(n = all_N_combs, K = all_K_combs, Empirical = emp_SEs,
#                      Mean = mean_delta_SEs, Median = median_delta_SEs,
#                      "Percent_Err_Mean" = 100*(mean_delta_SEs - emp_SEs)/emp_SEs,
#                      "Percent_Err_Median" = 100*(median_delta_SEs - emp_SEs)/emp_SEs)

```



```{r results="asis", echo=FALSE}
## Print table
kbl(results, format = "latex", booktabs = T, caption = "Summary of SE estimates for mediation effect under continuous response, continuous mediator, fixed-effects model \\label{tab:results_bin_bin_mix}")
```



There's lots that I don't understand about the above results. One think I decided to explore is the normality of the estimators. To this end, I make some density plots, and also perform Shapiro-Wilk tests for normality. See Table \ref{tab:results_bin_bin_mix-Shapiro} for the test results.

```{r}
med_hats_to_plot = unlist(data_med_hats)
Ns_to_plot = rep(N_K_combs$N, times = sapply(data_med_hats, length))
Ks_to_plot = rep(N_K_combs$K, times = sapply(data_med_hats, length))
data_med_hat_to_plot = data.frame(med_hat = med_hats_to_plot, N = Ns_to_plot, K = Ks_to_plot)

data_med_hat_to_plot_OLD = data_med_hat_to_plot
data_med_hat_to_plot %<>% dplyr::filter(med_hat < 4)
```


```{r}
ggplot(data_med_hat_to_plot, aes(x = med_hat)) +
  geom_density() +
  facet_grid(rows = vars(N), cols = vars(K), scales = "free") +
  ggtitle("Density of mediation effect estimates")
  # theme_minimal()
```


```{r}
# QQ plot for each sample size
ggplot(data_med_hat_to_plot, aes(sample = med_hat)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(rows = vars(N), cols = vars(K), scales = "free") +
  ggtitle("QQ plot of mediation effect estimates")
  # theme_minimal()
```


```{r}
# Test normality in each sample size
all_norm_p_vals = c()
for(i in 1:nrow(N_K_combs)){
  p_val = shapiro.test(data_med_hats[[i]])$p.value
  all_norm_p_vals = c(all_norm_p_vals, p_val)
}

results$Shapiro_p = all_norm_p_vals
# results

# # Check how p-value varies with N and K
# p_vals_by_N = results %>% group_by(n) %>% summarize(mean_p = mean(Shapiro_p), sd_p = sd(Shapiro_p))
# p_vals_by_K = results %>% group_by(K) %>% summarize(mean_p = mean(Shapiro_p), sd_p = sd(Shapiro_p))

```

```{r results="asis", echo=FALSE}
## Print table
kbl(results, format = "latex", booktabs = T, caption = "Summary of SE estimates for mediation effect under continuous response, continuous mediator, fixed-effects model \\label{tab:results_bin_bin_mix-Shapiro}")
```
