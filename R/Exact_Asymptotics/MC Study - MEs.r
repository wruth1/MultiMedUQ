
# Covariances of mediation effects


list_ME_hats = list()
list_ME_cov_hats = list()

this_scale = c("diff", "rat", "OR")

for(i in seq_along(all_Ks)){
    print(paste0("K = ", all_Ks[i], "; number ", i, " of ", length(all_Ks)))

    some_ME_hats = data.frame()
    some_ME_cov_hats = list()

    for(j in seq_len(num_reps)){

        this_par_hat = list_par_hats[[i]][j,]

        this_b_Y = this_par_hat[1:5]
        this_theta_Y = this_par_hat[6:8]
        this_b_M = this_par_hat[9:12]
        this_theta_M = this_par_hat[13:15]


        this_ME_hat = all_MEs_pars(this_scale, w, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
        some_ME_hats = rbind(some_ME_hats, this_ME_hat)


        this_Sigma = list_par_cov_hats[[i]][[j]]
        this_ME_cov_hat = all_covs_MEs_pars(this_scale, w, this_Sigma, this_b_Y, this_theta_Y, this_b_M, this_theta_M, which_REs=which_REs)
        some_ME_cov_hats[[j]] = this_ME_cov_hat

    }

    colnames(some_ME_hats) = c("11", "10", "01", "00")
    list_ME_hats[[i]] = some_ME_hats
    list_ME_cov_hats[[i]] = some_ME_cov_hats
}


## Summarize covariance estimates

### Empirical

list_ME_emp_covs = list()

for(i in seq_along(all_Ks)){
    list_ME_emp_covs[[i]] = cov(list_ME_hats[[i]])
}


### Estimated

list_ME_mean_covs = list()

for(i in seq_along(all_Ks)){
    some_ME_cov_hats = list_ME_cov_hats[[i]]
    this_mean_ME_cov_hat = Reduce("+", some_ME_cov_hats) / num_reps
    list_ME_mean_covs[[i]] = this_mean_ME_cov_hat
}



## Compare covariance estimates

### Easier: Variances (i.e. diagonal elements)
info_var_diffs = data.frame()
for(i in seq_along(all_Ks)){
    some_var_diffs = diag(list_ME_emp_covs[[i]] - list_ME_mean_covs[[i]])
    some_rel_var_diffs = some_var_diffs / diag(list_ME_emp_covs[[i]])

    this_info = c(all_Ks[i], mean(some_var_diffs), sd(some_var_diffs), mean(some_rel_var_diffs), sd(some_rel_var_diffs))

    info_var_diffs = rbind(info_var_diffs, this_info)
}
colnames(info_var_diffs) = c("K", "Mean_Diff", "SD_Diff", "Mean_Rel_Diff", "SD_Rel_Diff")

info_var_diffs

info_var_diffs %>% mutate(scale_abs_diff = Mean_Diff * K) %>% select(-"Mean_Rel_Diff", -"SD_Rel_Diff")


### Harder: Matrix norms
info_ME_norms = data.frame()
for(i in seq_along(all_Ks)){
    diff_mat = list_ME_emp_covs[[i]] - list_ME_mean_covs[[i]]
    diff_mat_norm = norm(diff_mat, type = "2")

    # rel_norm = diff_mat_norm / norm(list_ME_emp_covs[[i]], type = "2")
    rel_norm = diff_mat_norm / norm(list_ME_emp_covs[[i]], type = "F")

    this_info = c(all_Ks[i], diff_mat_norm, rel_norm)
    info_ME_norms = rbind(info_ME_norms, this_info)
}
colnames(info_ME_norms) = c("K", "Norm_Diff", "Rel_Norm_Diff")

info_ME_norms

scaled_ME_norms = info_ME_norms %>% mutate(scale_abs_diff = Norm_Diff * K) %>% select(-"Rel_Norm_Diff")



# Summarize errors for all covariance estimates

theta_mat = scaled_rel_norms %>% select(-Norm_Rel, -scale_rel_norm)
ENC_mat = scaled_ENC_norms
ME_mat = scaled_ME_norms

colnames(theta_mat) = c("K", "theta-abs", "theta-scaled")
colnames(ENC_mat) = c("K", "ENC-abs", "ENC-scaled")
colnames(ME_mat) = c("K", "ME-abs", "ME-scaled")

summary_mat = cbind(theta_mat, ENC_mat[,-1], ME_mat[,-1])

summary_mat %>% 
    dplyr::mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
    kbl(., format = "latex", booktabs = F, 
        caption = "Raw and scaled absolute errors in estimating empirical covariance matrices of estimators for the parameters, ENCs and mediation effects.",
        label = "tab:err_summary",
        col.names = c("K", "Abs $\\theta$", "Scaled $\\theta$", "Abs ENC", "Scaled ENC", "Abs ME", "Scaled ME"),
        align = "c", escape = F
    )

