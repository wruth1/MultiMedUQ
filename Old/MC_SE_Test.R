
# Simulated validation data ####

## Initialize variables

set.seed(1)

n = 50
K = 500


x = 0
x_m = 1

b_Y = c(0,0,1,0,0)
theta_Y = c(sqrt(0.5), 0.5, 0, 1, 0.5, sqrt(0.5))

b_M = c(0,0,0,0)
theta_M = c(1, 0.5, 2)



w = c(0,0)

scale = c("diff", "rat", "OR")



# Run simulation

num_reps = 10

all_MEs = c()
all_ME_covs = list()

tic()

for(i in 1:num_reps){
  print(paste0(i, " of ", num_reps))



## Generate data
data = make_validation_data(n, K, b_Y, theta_Y, b_M, theta_M, output_list = F)


  ## Note: glmer wasn't converging with default values. I chose one of the default optimizers, and increased the number of function evaluations. Both bobyqa and the other default use this limiter instead of the number of iterations.
  (fit_Y = lme4::glmer(Y ~ X + M + C1 + C2 + (X + M | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))
  (fit_M = lme4::glmer(M ~ X + C1 + C2 + (X | group), data = data, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))



  this_MEs = all_MEs_models(scale, w, fit_Y, fit_M)


  this_ME_covs = all_cov_MEs(scale, w, fit_Y, fit_M)


  all_MEs = rbind( all_MEs, this_MEs)
  all_ME_covs[[i]] = this_ME_covs

}

toc()

emp_ME_covs = cov(all_MEs)

sum_ME_covs = matrix(0,nrow(all_ME_covs[[1]]), ncol(all_ME_covs[[1]]))
for(i in 1:length(all_ME_covs)){
  sum_ME_covs = sum_ME_covs + all_ME_covs[[i]]
}
mean_ME_covs = sum_ME_covs / length(all_ME_covs)

norm(emp_ME_covs - mean_ME_covs) / norm(emp_ME_covs)
norm(emp_ME_covs - sqrt(K)*mean_ME_covs) / norm(emp_ME_covs)


plot(unlist(emp_ME_covs), unlist(mean_ME_covs))
abline(a = 0, b = 1)
